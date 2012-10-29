#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
from rsr_analysis import parser, main
if __name__ == '__main__':
    values = parser.parse_args()
    if not (values.pdbids or values.swissprot or values.pdbidfile):
        values = parser.parse_args(['--help'])
    main(values)
