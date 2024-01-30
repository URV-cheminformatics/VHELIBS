# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2024 Adrià Cereto Massagué <adria.cereto@fundacio.urv.cat>
#
"""
Module dealing with getting EDM data and validation stats
"""
import os
import sys
import time
import xml.etree.ElementTree as ET
import PDBfiles
import requests

edmapsurl = "https://edmaps.rcsb.org/maps/{}_2fofc.dsn6"
edstatsurl = "https://www.ebi.ac.uk/pdbe/entry-files/download/{}_validation.xml"


def get_EDM(pdbid):
    """
    Downloads the EDM of the given pdb code
    """
    edmurl = edmapsurl.format(pdbid.lower())
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    omapfile = os.path.join(downloaddir,  os.path.basename(edmurl))
    tries = 0
    while tries <= 3:
        try:
            if os.path.isfile(omapfile):
                sigma = 1
                return omapfile, sigma
            req = urlopen(edmurl)
            omap = req.read()
            req.close()
            outfile = open(omapfile, 'wb')
            outfile.write(omap)
            outfile.close()
        except Exception as e:
            print(e)
            tries += 1
    return None, None


def get_EDS(pdbid):
    """
    Extract data from EDS site for a given PDB code
    """
    pdbid = pdbid.lower()
    pdbdict = {pdbid: None}
    edd_dict = {}
    url = edstatsurl.format(pdbid.lower())
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    statfilepath = os.path.join(downloaddir,  os.path.basename(url))
    try:
        if not os.path.isfile(statfilepath):
            tries = 0

            print('Downloading %s' % url)
            statfilelines = []
            while tries <= 3:
                tries += 1
                try:
                    req = requests.get(url)
                    statfilelines = req.content
                    if not statfilelines:
                        print('could not read stat file')
                        pdbdict[pdbid] = False
                        return pdbdict, edd_dict
                    statfile = open(statfilepath, 'wb')
                    statfile.write(statfilelines)
                    statfile.close()
                    tries = 999
                except Exception as e:
                    if tries > 3:
                        print(e)
                        raise e
                    time.sleep(1)
        if not os.path.isfile(statfilepath):
            print('could not read stat file')
            pdbdict[pdbid] = False
            return pdbdict, edd_dict

        tree = ET.parse(statfilepath)
        for res in tree.findall("ModelledSubgroup"):

            resname = list("   ")
            for i, c in enumerate(res.get("resname")[::-1]):
                resname[2-i] = c
            resname = "".join(resname)

            resnum = list("    ")
            for i, c in enumerate(res.get("resnum")[::-1]):
                resnum[3-i] = c
            resnum = "".join(resnum)

            residue = "{} {}{}{}".format(resname, res.get(
                "chain"), resnum, res.get("icode")).strip()

            resdict = {"RSR": float(res.get("rsr", 100)), "RSCC": float(res.get("rscc", 0)), "OWAB": float(res.get("owab", 1000)), "RSRZ": float(res.get("rsrz", 9999)), "occupancy": float(res.get("avgoccu", 0))
                       }

            edd_dict[residue] = resdict

        pdbdict[pdbid] = True

    except Exception as e:
        if hasattr(e, 'reason'):
            pdbdict[pdbid] = e.reason
        elif hasattr(e, 'code'):
            if e.code == 404:
                pdbdict[pdbid] = False
            else:
                pdbdict[pdbid] = e.code
        else:
            pdbdict[pdbid] = str(e)
    return pdbdict, edd_dict
