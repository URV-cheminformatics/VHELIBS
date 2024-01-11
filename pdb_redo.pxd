import cython

cython.declare(PDBREDObase=str,CACHE_EXPIRED=cython.bint)

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
    , preresidue=list
    , residue=str
    , header=list
    , row=list
    , tries=cython.int
    )
cdef dict get_ED_data(str pdbid)

@cython.locals(
    alldatapath=str
    ,tries = cython.int
    ,data_started = cython.int
    ,download = cython.bint
    ,n = cython.bint
    ,firstlines = str
    ,row = list
    ,header = list
    ,datadict = dict
    ,parseddata = dict
)
cdef dict get_pdbredo_data(list pdbids=*)
