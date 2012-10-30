import cython

cdef update_lists(dict new_m=* , dict new_lb=*)

@cython.locals(metalstring=str, metalhetids=list, hetid=str, lblstring=str, lblhetids=list)
cdef cython.int dump_lists(str fname=*)

@cython.locals(new_m=dict, new_lb=dict, lines=list, line=str, dirty_hetid=str, dirty_desc=str)
cdef cython.int load_lists(str fname)
