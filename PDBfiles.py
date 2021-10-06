# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2015 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Handle and download PDB files
"""
import os, time, tempfile

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
#On guardarem els fitxers:
CACHEDIR = tempfile.mkdtemp()

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
        raise Exception("PDB_REDO not supported anymore")
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

