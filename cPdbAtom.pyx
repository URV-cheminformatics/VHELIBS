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
    cdef public float occupancy
    cdef public str variant
    cdef float[3] xyz
    def __init__(PdbAtom self, dict atom_dict):
        cdef str pos
        pos = atom_dict["auth_seq_id"]
        while len(pos) < 4:
            pos = " " + pos
        self.residue = "{} {}{}".format(atom_dict["auth_comp_id"],  atom_dict["auth_asym_id"], pos)
        self.hetid = atom_dict["auth_comp_id"]
        self.xyz[0] = float(atom_dict["Cartn_x"])
        self.xyz[1] = float(atom_dict["Cartn_y"])
        self.xyz[2] = float(atom_dict["Cartn_z"])
        self.occupancy = float(atom_dict["occupancy"])
        self.variant = atom_dict["label_alt_id"]

    def __or__(PdbAtom self, PdbAtom other):
        """
        Return squared distance
        """
        cdef float d
        d = distance(self, other)
        return d

cdef float distance(PdbAtom self, PdbAtom other):
    return (self.xyz[0] - other.xyz[0])**2 + (self.xyz[1] - other.xyz[1])**2 + (self.xyz[2] - other.xyz[2])**2
