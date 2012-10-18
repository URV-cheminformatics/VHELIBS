#-*- coding: utf-8 -*-
class PdbAtom(object):
    """
    Represents an atom from a PDB file
    """
    def __init__(self, record):
        """
        Needs an ATOM or HETATM record
        """
        self.name = record[12:16]
        self.residue = record[17:27]
        self.hetid = self.residue[:3].strip()
        self.xyz = (float(record[30:38]), float(record[38:46]), float(record[46:54]))
    def __or__(self, other):
        """
        Return squared distance
        """
        sx, sy, sz = self.xyz
        ox, oy, oz = other.xyz
        return (sx - ox)**2 + (sy - oy)**2 + (sz - oz)**2
