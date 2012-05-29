#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2011 - 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, sys, urllib2, csv, itertools
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
#    print inner_distance
#    print outer_distance
    #print argtuple
    pdbid, rsr_upper, rsr_lower = argtuple
    future_hetids_list = set()
    try:
        hetids_list = PDBfiles.hetdict[pdbid.lower()]
    except KeyError:
        hetids_list =[]
    residues_to_exam = set()
    ligand_residues = set()
    binding_sites = set()
    good_rsr = set()
    dubious_rsr = set()
    bad_rsr = set()
    protein_atoms = set()
    ligand_all_atoms_dict = {}
    ligand_res_atom_dict = {}
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
        seqres = set()
        links = []
        error = False
        for line in pdbfile:
            line = line.strip()
            label = line[:6].strip()
            if label in ("ATOM", "HETATM"):
                atom = PdbAtom(line)
                if atom.residue in rsrdict:
                    atom.rsr = float(rsrdict[atom.residue])
                if label == "ATOM" and inner_distance: #Don't care about protein when distance = 0
                    protein_atoms.add(atom)
                    if atom.name[1:3] == 'CA':  #Is alpha-carbon
                        protein_ca_atoms.add(atom)
                    seqres.add(atom.residue)
                elif label == 'HETATM':
                    if atom.hetid in hetids_list or (not hetids_list and (atom.hetid not in ligand_blacklist)):
                        ligand_residues.add(atom.residue)
                        if not hetids_list and atom.hetid not in ligand_all_atoms_dict:
                            future_hetids_list.add(atom.hetid)
                            ligand_all_atoms_dict[atom.hetid] = set()
                        ligand_all_atoms_dict[atom.hetid].add(atom)
                        if not atom.residue in ligand_res_atom_dict:
                            ligand_res_atom_dict[atom.residue] = set()
                        ligand_res_atom_dict[atom.residue].add(atom)
            elif label == 'LINK':
                links.append((line[17:27].strip(),  line[47:57].strip(), float(line[73:78]))) #distancia
#            elif label == 'SEQRES':
#                seqres.update(set(line[19:].split()))
    except IOError, error:
        print pdbfilepath
        print error
    finally:
        pdbfile.close()
        if error:
            return  (None, None, None, None)
    #Now let's prune covalently bound ligands
    notligands = set()
    alllinksparsed = True
    while alllinksparsed:
        alllinksparsed = False
        for res1,  res2,  blen in links:
            inseqres = 0
            ligres = sres = None
            if res1[:3] in seqres or res1 in seqres:
                inseqres +=1
                sres,  ligres = res1, res2
            if res2[:3] in seqres or res2 in seqres:
                inseqres +=1
                sres,  ligres = res2, res1
            if res1[:3] in ligand_blacklist or res2[:3] in ligand_blacklist:
                print 'Binding to a blacklisted ligand: %s - %s' % (res1, res2)
                notligands.add(res1)
                seqres.add(res1)
                ligres = res2
                inseqres = 1
            if inseqres == 1:
                if not blen or blen >= 2.1: #Disulfide bonds are about 2.05;
                    print 'Bond distance big enough (%s) between %s and %s' % (blen, res1,  res2)
                    continue
                notligands.add(ligres)
                seqres.add(ligres)
                links.remove((res1,  res2,  blen))
                alllinksparsed = True
                break
    for nonligand in notligands:
        print nonligand + 'is not a ligand!'
        for hetlist in (future_hetids_list, hetids_list):
            if nonligand[:3] in hetlist:
                hetlist.remove(nonligand[:3])
        if nonligand in ligand_residues:
            ligand_residues.remove(nonligand)
        if nonligand[:3] in ligand_all_atoms_dict:
            ligand_all_atoms_dict.pop(nonligand[:3])

    def classificate_residue(atom):
        #print 'comparing %s with upper %s' % ( atom.rsr ,  rsr_upper)
        if atom.rsr != None and atom.rsr <= rsr_upper:
            #print 'comparing %s with lower %s' % ( atom.rsr ,  rsr_lower)
            if atom.rsr <= rsr_lower:
                good_rsr.add(atom.residue)
                #print '%s is lower!' % atom.rsr
            else:
                #print 'added to dubious'
                dubious_rsr.add(atom.residue)
        else:
            #print 'added to bad'
            bad_rsr.add(atom.residue)

    for res in ligand_all_atoms_dict.itervalues():
        for atom in res:
            classificate_residue(atom)
            break

    def group_ligands(ligand_residues):
        """
        Split all the ligand residues into different molecules
        """
        print links
        ligands = []
        for lres in ligand_residues:
            ligand = set()
            ligand.add(lres)
            for res1,  res2,  blen in links:
                added = False
                if res1 == lres:
                    otherres = res2
                elif res2 == lres:
                    otherres = res1
                else:
                    otherres = None
                if otherres:
                    print lres, res1, res2
                    links.remove((res1,  res2,  blen))
                    if otherres in ligand_residues:
                        ligand.add(otherres)
                        for kligand in ligands:
                            if otherres in kligand:
                                kligand.update(ligand)
                                added = True
                                break
            else:
                if not added:
                    ligands.append(ligand)
                    added = True
        return ligands

    ligands = group_ligands(ligand_residues)

    def get_binding_site(ligand):
        """
        Get the binding site residues for the provided ligand and return them in a tuple
        """
        #Generate outer distance set
        outer_binding_site = set()
        for ligandres in ligand:
            for atom in protein_ca_atoms:
                for ligandatom in ligand_res_atom_dict[ligandres]:
                    distance = atom | ligandatom
                    if distance <= outer_distance:
                        outer_binding_site.add(atom.residue)
                        break
        #Generate inner distance set
        inner_binding_site = set()
        for atom in protein_atoms:
            if atom.residue not in outer_binding_site:
                continue
            for ligandatom in ligand_res_atom_dict[ligandres]:
                distance = atom | ligandatom
                if distance <= inner_distance:
                    inner_binding_site.add(atom.residue)
                    classificate_residue(atom)
                    break
        if not (inner_binding_site <= good_rsr and ligand_residues <= good_rsr):
            #Not all  the residues from here have good rsr
            residues_to_exam.update(dubious_rsr | bad_rsr)
        return inner_binding_site

    ligand_bs_list = []
    for ligand in ligands:
        bs = get_binding_site(ligand)
        rte = set()
        for resex in residues_to_exam:
            if resex in ligand or resex in bs:
                rte.add(resex)
        ligand_bs_list.append((ligand, bs, rte))
    print ligand_bs_list

    return (pdbid, ligand_bs_list)

class PdbAtom(object):
    """
    Represents an atom from a PDB file
    """
    def __init__(self, record):
        """
        Needs an ATOM or HETATM record
        """
        self.name = record[12:16]
        self.residue = record[17:27].strip()
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

def results_to_csv(results, outputfile):
    """
    Writes the output of parse_binding_site's to a csv file
    """
    outfile = open(outputfile, 'wb')
    csvfile = csv.writer(outfile)
    csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
    basename = os.path.splitext(os.path.basename(outputfile))[0]
    goodfile = None
    print 'Calculating...'
    datawritten = False
    for pdbid, ligand_bs_list in results:
        if pdbid == None:
            continue
        for ligandresidues, binding_site, residues_to_exam in ligand_bs_list:
            id = pdbid
            if not residues_to_exam:
                if not ligandresidues:
                    print '%s has no actual ligands, it will be discarded' % pdbid
                else:
                    if not goodfile:
                        goodfilename = basename + '_good.csv'
                        goodfile = open(goodfilename,'ab')
                        goodwriter = csv.writer(goodfile)
                        goodwriter.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
                    goodwriter.writerow([id, '', ';'.join(ligandresidues),';'.join(binding_site)])
                    goodfile.flush()
            else:
                csvfile.writerow([id, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
                outfile.flush()
                datawritten = True
    outfile.close()
    if not datawritten:
        os.remove(outputfile)
    if goodfile:
        print'Results for structures with a RSR value below the specified minimum were saved to %s' % goodfilename
    return datawritten

def main(filepath = None, pdbidslist=[], swissprotlist = [], rsr_upper=RSR_upper, rsr_lower = RSR_lower, distance=None, outputfile='rsr_analysis.csv'):
    if not rsr_upper > rsr_lower:
        print '%s is higher than %s!' % (rsr_lower, rsr_upper)
        raise ValueError
    if distance != None:
        global outer_distance
        global inner_distance
        distfactor = outer_distance/inner_distance
        inner_distance = distance**2
        outer_distance =  inner_distance*distfactor
    #print sys.argv[-1]
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
    argsarray = ((pdbid.upper(), rsr_upper, rsr_lower) for pdbid in pdblist if pdbid)
    if filepath:
        pdblistfile.close()
    #results = (parse_binding_site(argstuple) for argstuple in argsarray)
    #chunksize = int(math.sqrt(len(argsarray)))
    PDBfiles.setglobaldicts()
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    results = pool.imap(parse_binding_site, argsarray)
    datawritten = results_to_csv(results, outputfile)
    pool.terminate()
    pool.join()
    return datawritten

if __name__ == '__main__':
    values = parser.parse_args()
    #print values
    main(values.pdbidfile, pdbidslist = values.pdbids, swissprotlist =values.swissprot , rsr_upper=values.rsr_upper, rsr_lower = values.rsr_lower, distance=values.distance, outputfile = values.outputfile)
