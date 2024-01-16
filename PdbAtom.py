#-*- coding: utf-8 -*-
class PdbAtom(object):
    """
    Represents an atom from a PDB file
    """
    def __init__(self, atom_dict):
        """
        Needs an ATOM or HETATM record
        """
        self.name = atom_dict["auth_atom_id"]
        pos = atom_dict["auth_seq_id"]
        while len(pos) < 4:
            pos = " " + pos
        self.residue = "{} {}{}".format(atom_dict["auth_comp_id"],  atom_dict["auth_asym_id"], pos)
        self.hetid = atom_dict["auth_comp_id"]
        self.xyz = (float(atom_dict["Cartn_x"]), float(atom_dict["Cartn_y"]), float(atom_dict["Cartn_z"]))
        self.occupancy = float(atom_dict["occupancy"])
        self.variant = atom_dict["label_alt_id"]
    def __or__(self, other):
        """
        Return squared distance
        """
        sx, sy, sz = self.xyz
        ox, oy, oz = other.xyz
        return (sx - ox)**2 + (sy - oy)**2 + (sz - oz)**2
