# -*- coding: utf-8 -*-
#
#   Copyright 2013-2024 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
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
import  sys, os, time, csv, datetime
if sys.platform.startswith('java'):
    from java.util import Locale
import PDBfiles
import json

PDB_REDO_ed_data_url_tmpl = "https://pdb-redo.eu/db/PDBID/PDBID_final.json"
PDB_REDO_edm_url_tmpl = "https://pdb-redo.eu/db/PDBID/PDBID_final.mtz"
ALLDATA_URL = "https://pdb-redo.eu/db/PDBID/data.json"
CACHE_EXPIRED = True

def get_EDM(pdbid):
    """Downloads the EDM and returns its path"""
    pdbid = pdbid.lower()
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    url = PDB_REDO_edm_url_tmpl.replace('PDBID', pdbid)
    filename = os.path.join(downloaddir, pdbid + "_final.mtz")
    if not (os.path.isfile(filename) and os.path.getsize(filename) > 0): #Download
        print("Downloading %s" % url)
        tries = 3
        while tries > 0:
            tries -= 1
            try:
                rfh = urlopen(url)
                ed_file = open(filename, 'wb')
                ed_file.write((rfh.read()))
                ed_file.close()
                rfh.close()
                break
            except Exception as e:
                print("Could not download " +  url)
                print(e)
                print("Retrying... {}".format(tries))
                time.sleep(1)
        else:
            print("Unable to download %s" % url)
            return
    return filename

def get_ED_data(pdbid):
    """
    Extract data from PDB_REDO site for a given PDB code
    """
    pdbid = pdbid.lower()
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    url = PDB_REDO_ed_data_url_tmpl.replace('MIDDLE', pdbid[1:3]).replace('PDBID', pdbid)
    filename = os.path.join(downloaddir, os.path.basename(url))
    if not (os.path.isfile(filename) and os.path.getsize(filename) > 0): #Download
        print("Downloading %s" % url)
        tries = 3
        while tries > 0:
            tries -= 1
            try:
                rfh = urlopen(url)
                ed_file = open(filename, 'wb')
                ed_file.write(rfh.read())
                ed_file.close()
                rfh.close()
                break
            except Exception as e:
                print("Could not download " +  url)
                print(e)
                print("Retrying... {}".format(tries))
                time.sleep(1)
        else:
            print("Unable to download %s" % url)
            return
    #pdbdict={pdbid:None}
    edd_dict = {}
    edfile = open(filename, 'r')
    ed_data = json.load(edfile)
    for comp in ed_data:
        d = {"RSR":float(comp["RSR"] or 100)
                ,"RSCC": float(comp["RSCCS"] or 0)
    #            , "occupancy": 1.0 #FIXME
                             }
        hetid = comp['pdb']["compID"]
        while len(hetid) < 3:
            hetid = " " + hetid
        strand = comp['pdb']["strandID"]
        seqnum = str(comp['pdb']["seqNum"])
        while len(seqnum) < 4:
            seqnum = " " + seqnum
        residue = hetid + " " + strand + seqnum
        edd_dict[residue] = d
    return edd_dict

def get_pdbredo_data(pdbid):
    tries = 3
    url = ALLDATA_URL.replace("PDBID", pdbid)
    cachedir = os.path.join(PDBfiles.CACHEDIR, pdbid)
    try:
        os.makedirs(cachedir)
    except:
        pass
    alldatapath = os.path.join(cachedir, os.path.basename(url))
    #if not os.path.isfile(alldatapath): #Download
    while tries > 0:
        tries -= 1
        e = "unknown"
        if sys.platform.startswith('java'):
            oldlocale = Locale.getDefault()
            Locale.setDefault(Locale.ENGLISH)
        try:
            if os.path.isfile(alldatapath)  and os.path.getsize(alldatapath) > 0:
                print("loading {}".format(alldatapath))
                rawdict = json.load(open(alldatapath, 'rt'))
                break
            else:
                rawdict = json.load(urlopen(url))
                if not rawdict: return
                with open(alldatapath, "wt") as cache_file:
                        cache_file.write(json.dumps(rawdict))
                break
        except HTTPError as e:
            print("Could not download {}".format( url))
            print(e)
            print("Retrying... {}".format(tries))
            time.sleep(1)
        if sys.platform.startswith('java'):
            Locale.setDefault(oldlocale)
    else:
        print("Unable to download %s" % url)
        return
    #line 7 contains date
    #lines[11:111] are the columns explained
        #We need columns:
        #0   PDBID
        #14  RFIN
        #15  RFFIN
        #62  RESOLUTION
        # 65  NREFCNT    Number of reflections
        # 70  AAXIS      Length of cell axis a
        # 71  BAXIS      Length of cell axis b
        # 72  CAXIS      Length of cell axis c
        # 73  ALPHA      Cell angle alpha
        # 74  BETA       Cell angle beta
        # 75  GAMMA      Cell angle gamma
    #lines[112] == '#START DATA'
    #lines[113] are the column headers
    #here the data [114:-3]
    #lines[-3] == ''
    #lines[-2] == '#END DATA'
    rowdict = {}
    rowdict["experimentalTechnique"] = rawdict["properties"]["EXPTYP"]
    rowdict["rFree"] = rawdict["properties"]["RFFIN"]
    rowdict["rWork"] = rawdict["properties"]["RFIN"]
    rowdict["refinementResolution"] = rawdict["properties"]["RESOLUTION"]
    rowdict["unitCellAngleAlpha"] = rawdict["properties"]["ALPHA"]
    rowdict["unitCellAngleBeta"] = rawdict["properties"]["BETA"]
    rowdict["unitCellAngleGamma"] = rawdict["properties"]["GAMMA"]
    rowdict["lengthOfUnitCellLatticeA"] = rawdict["properties"]["AAXIS"]
    rowdict["lengthOfUnitCellLatticeB"] = rawdict["properties"]["BAXIS"]
    rowdict["lengthOfUnitCellLatticeC"] = rawdict["properties"]["CAXIS"]
    rowdict["nreflections"] = rawdict["properties"]["NREFCNT"]
    return rowdict

if __name__ == "__main__":
    print('start')
    d = get_pdbredo_data(["3dzu"])
    print(len(d))

