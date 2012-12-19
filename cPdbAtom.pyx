# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
cdef class PdbAtom(object):
    """
    Represents an atom from a PDB file
    """
    cdef public str residue
    cdef public str hetid
    cdef float[3] xyz
    def __init__(PdbAtom self, str record):
        """
        Needs an ATOM or HETATM record
        """
        self.residue = record[17:27]
        self.hetid = self.residue[:3].strip()
        self.xyz[0] = float(record[30:38])
        self.xyz[1] = float(record[38:46])
        self.xyz[2] = float(record[46:54])
    def __or__(PdbAtom self, PdbAtom other):
        """
        Return squared distance
        """
        cdef float d
        d = distance(self, other)
        return d

cdef float distance(PdbAtom self, PdbAtom other):
    return (self.xyz[0] - other.xyz[0])**2 + (self.xyz[1] - other.xyz[1])**2 + (self.xyz[2] - other.xyz[2])**2
