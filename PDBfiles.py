# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2024 Adrià Cereto Massagué <adria.cereto@fundacio.urv.cat>
#
"""
Handle and download PDB files
"""
import os, time, tempfile, json

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
    import ssl
    ssl._create_default_https_context = ssl._create_unverified_context
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

PDBbase = "http://www.rcsb.org/pdb/files/{}.pdb.gz"
PDBREDObase_full = "https://pdb-redo.eu/db/PDBID/PDBID_final.pdb"
#On guardarem els fitxers:
CACHEDIR = tempfile.mkdtemp()

#SERVICELOCATION="http://www.rcsb.org/pdb/rest/customReport"

#QUERY_TPL = "?pdbids=%s&customReportColumns={}&service=wsfile&format=csv".format(",".join(columns))
QUERY_TPL = "https://data.rcsb.org/rest/v1/core/entry/{}"

def get_custom_report(pdbid):
    urlstring = QUERY_TPL.format(pdbid)
    print(urlstring)
    rawdict = json.load(urlopen(urlstring))
    rowdict = {}
    rowdict["experimentalTechnique"] = rawdict["rcsb_entry_info"]["experimental_method"]
    rowdict["rFree"] = rawdict["refine"][0]["ls_rfactor_rfree"]
    rowdict["rWork"] = rawdict["refine"][0]["ls_rfactor_rwork"]
    rowdict["refinementResolution"] = rawdict["refine"][0]["ls_dres_high"]
    rowdict["nreflections"] = rawdict["refine"][0]["ls_number_reflns_rfree"]
    rowdict["unitCellAngleAlpha"] = rawdict["cell"]["angle_alpha"]
    rowdict["unitCellAngleBeta"] = rawdict["cell"]["angle_beta"]
    rowdict["unitCellAngleGamma"] = rawdict["cell"]["angle_gamma"]
    rowdict["lengthOfUnitCellLatticeA"] = rawdict["cell"]["length_a"]
    rowdict["lengthOfUnitCellLatticeB"] = rawdict["cell"]["length_b"]
    rowdict["lengthOfUnitCellLatticeC"] = rawdict["cell"]["length_c"]
    return {pdbid.upper():rowdict}

def get_pdb_file(pdbcode, pdb_redo = False):
    """
    Downloads a PDB file and stores it with the specified filename
    """
    global CACHEDIR
    if not os.path.isdir(CACHEDIR):
        os.makedirs(CACHEDIR)
    if not pdb_redo:
        url = PDBbase.format(pdbcode)
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".pdb.gz")
        print("Downloading {} to {}".format(url, filename))
    else:
         pdbcode = pdbcode.lower()
         url = PDBREDObase_full.replace('PDBID', pdbcode)
         filename = os.path.join(CACHEDIR, os.path.basename(url))
    if os.path.isfile(filename):
        return filename
    tries = 0
    while tries <= 3:
        tries += 1
        try:
            handler = urlopen(url)
            filehandle = open(filename, "wb")
            filehandle.write(handler.read())
            filehandle.close()
            handler.close()
            return filename
        except Exception as e:
            print("Could not download",  url)
            print(e)
            print("Retrying...",  tries)
            time.sleep(1)
    print("%s could not be downloaded!!!!" % url)
    return ''

########################################################

