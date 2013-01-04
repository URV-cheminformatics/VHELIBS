# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Handle and download PDB files
"""
import os, gzip, urllib2, time, sys, tempfile
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
else:
    import multiprocessing
from cofactors import ligand_blacklist

PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"
#On guardarem els fitxers:
CACHEDIR = tempfile.mkdtemp()

def get_pdb_file(pdbcode, filename = ''):
    """
    Downloads a PDB file and stores it with the specified filename
    """
    global CACHEDIR
    if not os.path.isdir(CACHEDIR):
        os.makedirs(CACHEDIR)
    url = PDBbase % pdbcode

    if not filename:
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".pdb.gz")
    tries = 0
    downloaded = False
    while tries <= 3 and not downloaded:
        tries += 1
        try:
            if os.path.isfile(filename):
                return filename
            handler = urllib2.urlopen(url)
            filehandle = open(filename, "wb")
            filehandle.write(handler.read())
            filehandle.close()
            handler.close()
            downloaded = True
            return filename
        except Exception, e:
            print "Could not download",  url
            print e
            print "Retrying...",  tries
            time.sleep(1)
    if tries > 3:
        print "%s could not be downloaded!!!!" % url
        return ''

########################################################

