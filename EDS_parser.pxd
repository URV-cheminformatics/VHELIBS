import cython

cython.declare(edsurl=str)

@cython.locals(
    pdbdict=dict
    , edd_dict=dict
    , downloaddir=str
    , statfilepath=str
    , statfilelines=list
    , tries=cython.int
    , line=str
    , rscc=str
    , rsr=str
    , owab=str
    , natom=str
    , s_occ=str
    , residue=str
    )
cdef tuple get_EDS(str pdbid)
