#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import cython
from cpython cimport bool
import argparse
import csv

cython.declare(
    RSR_upper=cython.float
    , RSR_lower=cython.float
    , inner_distance=cython.float
    , titles=list
    )

cpdef int dbg(str string)

@cython.locals(sptopdb_dict=dict,temppdbdict=dict)
cpdef dict get_sptopdb_dict()

@cython.locals(datawritten=bool, restuplelen=cython.int, pdbid=str)
cpdef bool results_to_csv(results, outputfile)
