#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import cython
from cpython cimport bool
cimport PDBfiles
cimport cofactors
cimport EDS_parser

cython.declare(
    RSR_upper=cython.float
    , RSR_lower=cython.float
    , inner_distance=cython.float
    , titles=list
    )

cdef cython.int dbg(str string)

@cython.locals(sptopdb_dict=dict,temppdbdict=dict)
cdef dict get_sptopdb_dict()

@cython.locals(
pdbid=str
,rsr_upper=cython.float
,rsr_lower=cython.float
,future_hetids_list=set
,hetids_list=list
,ligand_residues = set
,binding_sites = set
,good_rsr = set
,dubious_rsr = set
,bad_rsr = set
,protein_atoms = set
,ligand_all_atoms_dict = dict
,ligand_res_atom_dict = dict
,notligands = dict
,seqres = set
,links = list
,line = str
,label = str
,alllinksparsed = cython.bint
,ligdiff=set
,checklink = cython.int)
cpdef tuple parse_binding_site(tuple argtuple)

@cython.locals(rsr=cython.float, Natom=cython.float, S_occ=cython.float)
cdef cython.int classificate_residue(
    residue
    , dict edd_dict
    , set good_rsr
    , set dubious_rsr
    , set bad_rsr
    , cython.float rsr_upper
    ,  cython.float rsr_lower
    )

@cython.locals(
ligands = list
,ligand_links = list
,linked_ligand_res = set
,ligand=set)
cdef list group_ligands(set ligand_residues, list links)

@cython.locals(inner_binding_site = set, rte=set)
cdef tuple get_binding_site(set ligand, set good_rsr, set bad_rsr, set dubious_rsr, str pdbid, set protein_atoms, ligands, dict ligand_res_atom_dict, cython.float rsr_upper, cython.float rsr_lower, dict edd_dict)

cdef object validate(set residues, set good_rsr, set bad_rsr, set dubious_rsr, str pdbid)

@cython.locals(datawritten=cython.bint, restuplelen=cython.int, pdbid=str)
cdef bool results_to_csv(results, outputfile)

cpdef main(values)
