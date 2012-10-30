import cython

cython.declare(PDBbase=str)

@cython.locals(pdbcode=str, url=str, tries=cython.int, downloaded=cython.bint)
cdef str get_pdb_file(str pdbcode, str filename=*)

@cython.locals( outdict=dict, destfile=str, dicturl=str, line=str, ligand=str,  pdb_codes=str, pdb_code=str)
cdef dict get_ligand_pdb_dict(blacklist=*)

cpdef cython.int setglobaldicts()
