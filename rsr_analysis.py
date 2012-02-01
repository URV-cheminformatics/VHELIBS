#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2011 - 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, sys
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
else:
    import multiprocessing

import PDBfiles, EDS_parser

RSR_upper = 0.4
RSR_lower = 0.24
outer_distance = 10**2
inner_distance = 4.5**2

def get_pdbs_with_good_rsr(pdblist, rsr_upper=RSR_upper, rsr_lower = RSR_lower):
    if not rsr_upper > rsr_lower:
        print '%s is higher than %s!' % (rsr_lower, rsr_upper)
        raise ValueError
    argsarray = []
    resultdict = {}
    for pdbid in pdblist:
        argsarray.append((pdbid.upper(), rsr_upper, rsr_lower))
    PDBfiles.setglobaldicts()
    #pool = multiprocessing.Pool(multiprocessing.cpu_count())
    #results = pool.map(parse_binding_site, argsarray)
    results = (parse_binding_site(argstuple) for argstuple in argsarray)
    for pdbid, ligandresidues, residues_to_exam, binding_site in results:
        resultdict[pdbid] = ligandresidues, residues_to_exam, binding_site
    return resultdict

def parse_binding_site(argtuple):
    """
    argtuple = (pdbid, rsr_upper, rsr_lower)
    """
    #print argtuple
    pdbid, rsr_upper, rsr_lower = argtuple
    hetids_list = PDBfiles.hetdict[pdbid.lower()]
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
                elif atom.hetid in hetids_list:
                    ligand_residues.add(atom.residue)
                    ligand_all_atoms_dict[atom.hetid].add(atom)
    except IOError, error:
        print pdbfilepath
        print error
    finally:
        pdbfile.close()
        if error:
            return  (None, None, None, None)
    #Now let's find binding site distance
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
            classificate_residue(atom, good_rsr, dubious_rsr, bad_rsr, rsr_upper=rsr_upper, rsr_lower = rsr_lower)
            for ligandatom in ligand_all_atoms_dict[hetid]:
                distance = atom | ligandatom
                if distance <= outer_distance:
                    inner_binding_site.add(atom.residue)
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


def main(argv=[]):
    import csv, math
    #print sys.argv[-1]
    PDBfiles.setglobaldicts()
    rsr_upper=RSR_upper
    rsr_lower = RSR_lower
    pdblistfile = open(argv[-1], 'rb')
    pdblist = (line.strip() for line in pdblistfile)
    if not rsr_upper > rsr_lower:
        print '%s is higher than %s!' % (rsr_lower, rsr_upper)
        raise ValueError
    argsarray = []
    resultdict = {}
    outfile = open('rsr_analysis.csv', 'wb')
    csvfile = csv.writer(outfile)
    csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
    for pdbid in pdblist:
        argsarray.append((pdbid.upper(), rsr_upper, rsr_lower))
    pdblistfile.close()
    #results = (parse_binding_site(argstuple) for argstuple in argsarray)
    chunksize = int(math.sqrt(len(argsarray)))
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    results = pool.imap_unordered(parse_binding_site, argsarray, chunksize)
    print 'Calculating...'
    for pdbid, ligandresidues, residues_to_exam, binding_site in results:
        if pdbid == None:
            continue
        if not residues_to_exam:
            continue
        csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
    outfile.close()
    pool.terminate()
    pool.join()

if __name__ == '__main__':
    main(sys.argv)
