# -*- coding: utf-8 -*-
#
#   Copyright 2012-2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import cython

cython.declare(PDBbase=str, PDBREDObase=str)

@cython.locals(pdbcode=str, url=str, tries=cython.int)
cdef str get_pdb_file(str pdbcode, str filename=*)
