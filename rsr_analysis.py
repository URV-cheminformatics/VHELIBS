#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2011 - 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, sys, urllib2
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
else:
    import multiprocessing

import PDBfiles, EDS_parser
from cofactors import ligand_blacklist

RSR_upper = 0.4
RSR_lower = 0.24
outer_distance = 10**2
inner_distance = 4.5**2

###Create the argument parser###
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--pdbids', nargs='+', default=[], type=str, metavar='PDBID', help='list of PDB ids')
parser.add_argument('-s','--swissprot', nargs='+', default=[], type=str, metavar='SP_AN [or SP_NAME]', help='list of Swiss-Prot protein names or accession numbers')
parser.add_argument('-u','--rsr_upper', type=float, default=RSR_upper, metavar='FLOAT', help='set maximum RSR value for each residue (residues with a higher RSR will be discarded)')
parser.add_argument('-l','--rsr_lower', type=float, default=RSR_lower, metavar='FLOAT', help='set minimum RSR value for each residue (residues with a lower RSR value will be directly considered right)')
parser.add_argument('-d','--distance', type=float, default=4.5, metavar='Å', help='consider part of the binding sites all the residues nearer than this to the ligand (in Å)')
parser.add_argument('-f','--pdbidfile', metavar='PATH', type=unicode, default=None, required=False, help='text file containing a list of PDB ids, one per line')
parser.add_argument('-o','--outputfile', metavar='PATH', type=unicode, default='rsr_analysis.csv', required=False, help='output file name')
#######################

def get_sptopdb_dict():
    """
    Returns a dictionary containing the pdb entries for each Swiss-Prot entry
    """
    url = "http://www.uniprot.org/docs/pdbtosp.txt"
    sptopdb_dict = {}
    temppdbdict = {}
    reader = urllib2.urlopen(url)
    pdbid = None
    for line in reader:
        if not len(line) > 28:
            continue
        if line[:28].strip():
            line_stripped = line.strip()
            if line_stripped:
                pdbid_candidate = line_stripped.split()[0]
                if len(pdbid_candidate) == 4 and pdbid_candidate.isupper():
                    pdbid = pdbid_candidate
                    if pdbid not in temppdbdict:
                        temppdbdict[pdbid] = set()
        if pdbid:
            spinfo = line[28:].strip()
            temppdbdict[pdbid].update([item.strip() for item in spinfo.split(',') if item])
    reader.close()
    for pdbid in temppdbdict:
        for sp_id in temppdbdict[pdbid]:
            if sp_id not in sptopdb_dict:
                sptopdb_dict[sp_id] = set()
            sptopdb_dict[sp_id].add(pdbid)
    return sptopdb_dict

def parse_binding_site(argtuple):
    """
    argtuple = (pdbid, rsr_upper, rsr_lower)
    """
    #print argtuple
    pdbid, rsr_upper, rsr_lower = argtuple
    try:
        hetids_list = PDBfiles.hetdict[pdbid.lower()]
    except KeyError:
        hetids_list =[]
        future_hetids_list = set()
    residues_to_exam = set()
    ligand_residues = set()
    binding_sites = set()
    good_rsr = set()
    dubious_rsr = set()
    bad_rsr = set()
    protein_atoms = set()
    ligand_all_atoms_dict = {}
    for hetid in hetids_list:
        #print hetid
        ligand_all_atoms_dict[hetid] = set()
    protein_ca_atoms = set()
    pdbfilepath = os.path.join(PDBfiles.PREFIX, pdbid.lower(), pdbid.upper() + ".pdb.gz")
    pdbdict, rsrdict = EDS_parser.get_EDS(pdbid)
    if pdbdict['IN_EDS'] != 'TRUE':
        print "No EDS data available for %s, it will be discarded" % pdbid
        return  (None, None, None, None)
    if not os.path.isfile(pdbfilepath):
        PDBfiles.get_pdb_file(pdbid.upper(), pdbfilepath)
    pdbfile = gzip.GzipFile(pdbfilepath)
    #print pdbfilepath
    try:
        error = False
        for line in pdbfile:
            line = line.strip()
            if line[:6] in ("ATOM  ", "HETATM"):
                atom = PdbAtom(line)
                if atom.residue in rsrdict:
                    atom.rsr = float(rsrdict[atom.residue])
                if line[:6] == "ATOM  ":
                    protein_atoms.add(atom)
                    if atom.name[1:3] == 'CA':  #Is alpha-carbon
                        protein_ca_atoms.add(atom)
                elif atom.hetid in hetids_list or (not hetids_list and atom.hetid not in ligand_blacklist):
                    ligand_residues.add(atom.residue)
                    if not hetids_list and atom.hetid not in ligand_all_atoms_dict:
                        future_hetids_list.add(atom.hetid)
                        ligand_all_atoms_dict[atom.hetid] = set()
                    ligand_all_atoms_dict[atom.hetid].add(atom)
    except IOError, error:
        print pdbfilepath
        print error
    finally:
        pdbfile.close()
        if error:
            return  (None, None, None, None)
    #Now let's find binding site distance
    if not hetids_list:
        hetids_list = future_hetids_list
    for hetid in hetids_list:
        #print hetid
        #Generate outer distance set
        outer_binding_site = set()
        for atom in protein_ca_atoms:
            for ligandatom in ligand_all_atoms_dict[hetid]:
                distance = atom | ligandatom
                if distance <= outer_distance:
                    outer_binding_site.add(atom.residue)
                    break
        #Generate inner distance set
        inner_binding_site = set()
        for atom in protein_atoms:
            if atom.residue not in outer_binding_site:
                continue
            for ligandatom in ligand_all_atoms_dict[hetid]:
                distance = atom | ligandatom
                if distance <= inner_distance:
                    inner_binding_site.add(atom.residue)
                    classificate_residue(atom, good_rsr, dubious_rsr, bad_rsr, rsr_upper=rsr_upper, rsr_lower = rsr_lower)
                    break
        if not inner_binding_site <= good_rsr:
            #Not all  the residues from here have good rsr
            residues_to_exam.update(dubious_rsr | bad_rsr)
        binding_sites.update(inner_binding_site)
    #print ligand_residues
    return (pdbid, ligand_residues, residues_to_exam, binding_sites)



def classificate_residue(atom, good_rsr, dubious_rsr, bad_rsr, rsr_upper=RSR_upper, rsr_lower = RSR_lower):
    if atom.rsr != None and atom.rsr <= rsr_upper:
        if atom.rsr <= rsr_lower:
            good_rsr.add(atom.residue)
        else:
            dubious_rsr.add(atom.residue)
    else:
        bad_rsr.add(atom.residue)

class PdbAtom(object):
    """
    Represents an atom from a PDB file
    """
    def __init__(self, record):
        """
        Needs an ATOM or HETATM record
        """
        self.name = record[12:16]
        self.residue = record[17:27]
        self.hetid = self.residue[:3].strip()
        self.xyz = (float(record[30:38]), float(record[38:46]), float(record[46:54]))
        self.rsr = None
    def __or__(self, other):
        """
        Return squared distance
        """
        if not hasattr(other, 'xyz') or len(other.xyz) != 3:
            raise ValueError('invalid coordinates')
        else:
            return (self.xyz[0] - other.xyz[0])**2 + (self.xyz[1] - other.xyz[1])**2 + (self.xyz[2] - other.xyz[2])**2


def main(filepath, pdbidslist=[], swissprotlist = [], rsr_upper=RSR_upper, rsr_lower = RSR_lower, distance=None, outputfile='rsr_analysis.csv'):
    import csv, itertools
    if not rsr_upper > rsr_lower:
        print '%s is higher than %s!' % (rsr_lower, rsr_upper)
        raise ValueError
    if distance != None:
        global outer_distance
        global inner_distance
        distfactor = outer_distance/inner_distance
        inner_distance = distance**2
        if not outer_distance > inner_distance:
           outer_distance =  inner_distance*distfactor
    #print sys.argv[-1]
    PDBfiles.setglobaldicts()
    pdblist = pdbidslist
    if swissprotlist:
        sptopdb_dict = get_sptopdb_dict()
        for swissprot_id in swissprotlist:
            for key in sptopdb_dict:
                if swissprot_id in key:
                    pdblist = itertools.chain(pdblist, sptopdb_dict[key])
    if filepath:
        pdblistfile = open(filepath, 'rb')
        pdblist = itertools.chain(pdblist, (line.strip() for line in pdblistfile if line.strip()))
    if not rsr_upper > rsr_lower:
        print '%s is higher than %s!' % (rsr_lower, rsr_upper)
        raise ValueError
    resultdict = {}
    outfile = open(outputfile, 'wb')
    csvfile = csv.writer(outfile)
    csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
    basename = os.path.splitext(os.path.basename(outputfile))[0]
    goodfile = None
    argsarray = ((pdbid.upper(), rsr_upper, rsr_lower) for pdbid in pdblist if pdbid)
    if filepath:
        pdblistfile.close()
    #results = (parse_binding_site(argstuple) for argstuple in argsarray)
    #chunksize = int(math.sqrt(len(argsarray)))
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    results = pool.imap_unordered(parse_binding_site, argsarray)
    print 'Calculating...'
    for pdbid, ligandresidues, residues_to_exam, binding_site in results:
        if pdbid == None:
            continue
        if not residues_to_exam:
            if not ligandresidues:
                print '%s has no actual ligands, it will be discarded' % pdbid
            else:
                if not goodfile:
                    goodfilename = basename + '_good.csv'
                    goodfile = open(goodfilename,'ab')
                    goodwriter = csv.writer(goodfile)
                    goodwriter.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
                goodwriter.writerow([pdbid, '', ';'.join(ligandresidues),';'.join(binding_site)])
        else:
            csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
    outfile.close()
    pool.terminate()
    pool.join()

if __name__ == '__main__':
    values = parser.parse_args()
    print values
    main(values.pdbidfile, pdbidslist = values.pdbids, swissprotlist =values.swissprot , rsr_upper=values.rsr_upper, rsr_lower = values.rsr_lower, distance=values.distance, outputfile = values.outputfile)
