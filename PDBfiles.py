# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2024 Adrià Cereto Massagué <adria.cereto@fundacio.urv.cat>
#
"""
Handle and download PDB files
"""
import os
import time
import tempfile
import json
import sys
import requests

PDBbase = "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/dz/{}.cif.gz"
PDBREDObase_full = "https://pdb-redo.eu/db/PDBID/PDBID_final.cif"
# On guardarem els fitxers:
CACHEDIR = tempfile.mkdtemp()

# SERVICELOCATION="http://www.rcsb.org/pdb/rest/customReport"

# QUERY_TPL = "?pdbids=%s&customReportColumns={}&service=wsfile&format=csv".format(",".join(columns))
QUERY_TPL = "https://data.rcsb.org/rest/v1/core/entry/{}"


def get_custom_report(pdbid):
    urlstring = QUERY_TPL.format(pdbid)
    # print(urlstring)
    cachedir = os.path.join(CACHEDIR, pdbid)
    try:
        os.makedirs(cachedir)
    except:
        pass
    alldatapath = os.path.join(cachedir, "pdb_stats.json")
    try:
        if os.path.isfile(alldatapath) and os.path.getsize(alldatapath) > 0:
            print("loading {}".format(alldatapath))
            rawdict = json.load(open(alldatapath, 'rt'))
        else:
            print("Downloading {}".format(urlstring))
            rawdict = requests.get(urlstring).json()
            if rawdict:
                with open(alldatapath, "wt") as cache_file:
                    cache_file.write(json.dumps(rawdict))
            else:
                print("Empty stats for {}?".format(pdbid))
                return {}
        rawdict["rcsb_entry_info"]
        rawdict["refine"]
        rawdict["cell"]
    except Exception as e:
        print(e)
        return {}
    rowdict = {}
    rowdict["experimentalTechnique"] = rawdict["rcsb_entry_info"]["experimental_method"]
    rowdict["rFree"] = rawdict["refine"][0].get("ls_rfactor_rfree", 9999)
    rowdict["rWork"] = rawdict["refine"][0].get("ls_rfactor_rwork", 9999)
    rowdict["refinementResolution"] = rawdict["refine"][0]["ls_dres_high"]
    rowdict["nreflections"] = rawdict["refine"][0].get(
        "ls_number_reflns_rfree", 0)
    rowdict["unitCellAngleAlpha"] = rawdict["cell"]["angle_alpha"]
    rowdict["unitCellAngleBeta"] = rawdict["cell"]["angle_beta"]
    rowdict["unitCellAngleGamma"] = rawdict["cell"]["angle_gamma"]
    rowdict["lengthOfUnitCellLatticeA"] = rawdict["cell"]["length_a"]
    rowdict["lengthOfUnitCellLatticeB"] = rawdict["cell"]["length_b"]
    rowdict["lengthOfUnitCellLatticeC"] = rawdict["cell"]["length_c"]
    return {pdbid.upper(): rowdict}


def get_pdb_file(pdbcode, pdb_redo=False):
    """
    Downloads a PDB file and stores it with the specified filename
    """
    global CACHEDIR
    if not os.path.isdir(CACHEDIR):
        os.makedirs(CACHEDIR)
    if not pdb_redo:
        url = PDBbase.format(pdbcode.lower())
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".cif.gz")
    else:
        pdbcode = pdbcode.lower()
        url = PDBREDObase_full.replace('PDBID', pdbcode)
        filename = os.path.join(CACHEDIR, os.path.basename(url))
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        print("loading {}".format(filename))
        return filename
    print("Downloading {} to {}".format(url, filename))
    tries = 3
    while tries > 0:
        tries -= 1
        try:
            filehandle = open(filename, "wb")
            filehandle.write(requests.get(url).content)
            filehandle.close()
            return filename
        except Exception as e:
            print("Could not download",  url)
            print(e)
            print("Retrying...",  tries)
            time.sleep(1)
    print("%s could not be downloaded!!!!" % url)
    return ''

########################################################
