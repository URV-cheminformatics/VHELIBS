#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2010, 2011 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, multiprocessing
import db_updater, PDBfiles, EDS_parser

from config import host,  user,  password,  database


RSR_upper = 0.4
RSR_lower = 0.24
outer_distance = 10**2
inner_distance = 4.5**2
sql_query = """SELECT DISTINCT pdb_ligands.PDB_ID, pdb_ligands.HET_ID, pdb.UNIPROT_ACCN FROM pdb_ligands inner join pdb
    ON pdb_ligands.PDB_ID = pdb.PDB_ID
WHERE RSR < %s """

def get_pdbs_with_good_rsr(rsr_upper=RSR_upper):
    """
    """
    pdbid_listdict = {}
    db = db_updater.myDB(host,  user,  password,  database)
    queryresults = db.execute(sql_query % rsr_upper)
    for queryrow in queryresults:
        pdbid, hetid, accn = queryrow
        if not pdbid in pdbid_listdict:
            pdbid_listdict[pdbid] = (pdbid, [hetid, ], accn)
        else:
            pdbid_listdict[pdbid][1].append(hetid)

    pool = multiprocessing.Pool()
    results = pool.map(parse_binding_site, pdbid_listdict.values())
    #results = (parse_binding_site(value) for value in pdbid_listdict.itervalues())
    resultdict = {}
    for pdbid, ligandresidues, residues_to_exam, binding_site in results:
        resultdict[pdbid] = ligandresidues, residues_to_exam, binding_site
    return resultdict

def parse_binding_site(queryrow):
    """
    queryrow = (PDB_ID, [HET_ID,], UNIPROT_ACCN)
    """
    pdbid, hetids_list, accn = queryrow
    residues_to_exam = set()
    ligand_residues = set()
    binding_sites = set()
    good_rsr = set()
    dubious_rsr = set()
    bad_rsr = set()
    protein_atoms = set()
    ligand_all_atoms_dict = {hetid:set() for hetid in hetids_list}
    protein_ca_atoms = set()
    pdbfilepath = os.path.join(PDBfiles.PREFIX, accn, pdbid + ".pdb.gz")
    pdbdict, rsrdict = EDS_parser.get_EDS(pdbid)
    if not os.path.isfile(pdbfilepath):
        print('%s is not a file!' % pdbfilepath)
        return None
    pdbfile = gzip.GzipFile(pdbfilepath)
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
    pdbfile.close()
    #Now let's find binding site distance
    for hetid in hetids_list:
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
            classificate_residue(atom, good_rsr, dubious_rsr, bad_rsr)
            for ligandatom in ligand_all_atoms_dict[hetid]:
                distance = atom | ligandatom
                if distance <= outer_distance:
                    inner_binding_site.add(atom.residue)
                    break
        if not inner_binding_site <= good_rsr:
            #Not all  the residues from here have good rsr
            residues_to_exam.update(dubious_rsr | bad_rsr)
        binding_sites.update(inner_binding_site)
    return (pdbid, ligand_residues, residues_to_exam, binding_sites)



def classificate_residue(atom, good_rsr, dubious_rsr, bad_rsr):
    if atom.rsr != None and atom.rsr <= RSR_upper:
        if atom.rsr <= RSR_lower:
            good_rsr.add(atom.residue)
        else:
            dubious_rsr.add(atom.residue)
    else:
        bad_rsr.add(atom.residue)



#                if atom.name[1:3] == 'CA': #Is alpha-carbon
#                    pass


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


if __name__ == '__main__':
    import csv
    resultdict = get_pdbs_with_good_rsr()
    outfile = open('rsr_analysis.csv', 'wb')
    csvfile = csv.writer(outfile)
    csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
    for pdbid in resultdict:
        ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
        if not residues_to_exam:
            continue
        csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
    outfile.close()





