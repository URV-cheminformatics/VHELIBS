#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2012-2020 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import sys; sys.path.append("one-jar/src/pylib/")
import rsr_analysis
if __name__ == '__main__':
    values = rsr_analysis.parser.parse_args()
    if not (values.pdbids or values.swissprot or values.pdbidfile):
        values = rsr_analysis.parser.parse_args(['--help'])
    rsr_analysis.main(values)
