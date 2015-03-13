# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2015 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Handle and download PDB files
"""
import os, urllib2, time, tempfile

PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"
#PDBREDObase_conservative = "http://www.cmbi.ru.nl/pdb_redo/MIDDLE/PDBID/PDBID_besttls.pdb"
PDBREDObase_full = "http://www.cmbi.ru.nl/pdb_redo/MIDDLE/PDBID/PDBID_final.pdb"
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
        url = PDBbase % pdbcode
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".pdb.gz")
    else:
        pdbcode = pdbcode.lower()
        url = PDBREDObase_full.replace('MIDDLE', pdbcode[1:3]).replace('PDBID', pdbcode)
        filename = os.path.join(CACHEDIR, os.path.basename(url))
    if os.path.isfile(filename):
        return filename
    tries = 0
    while tries <= 3:
        tries += 1
        try:
            handler = urllib2.urlopen(url)
            filehandle = open(filename, "wb")
            filehandle.write(handler.read())
            filehandle.close()
            handler.close()
            return filename
        except Exception, e:
            print "Could not download",  url
            print e
            print "Retrying...",  tries
            time.sleep(1)
    print "%s could not be downloaded!!!!" % url
    return ''

########################################################

