import cython

cython.declare(edsurl=str, residuelist=str)

@cython.locals(pdbdict=dict, rsrdict=dict, downloaddir=str, statfilepath=str, statfilelines=list, tries=cython.int, line=str, rsr=str, residue=str)
cpdef tuple get_EDS(str pdbid)

cpdef tuple parse_EDS(pdblist, ligandlist=*)
