# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
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

hetdict = None

def get_pdb_file(pdbcode, filename = None):
    """
    Downloads a PDB file and stores it with the specified filename
    """
    url = PDBbase % pdbcode

    if not filename:
        filename = os.path.split(url)[1]
    tries = 0
    downloaded = False
    while tries <= 3 and not downloaded:
        tries += 1
        try:
            if os.path.isfile(filename):
                return
            handler = urllib2.urlopen(url)
            filehandle = open(filename, "wb")
            filehandle.write(handler.read())
            filehandle.close()
            handler.close()
            downloaded = True
            return
            tries = 999
        except Exception, e:
            print "Could not download",  url
            print e
            print "Retrying...",  tries
            time.sleep(1)
    if tries > 3:
        print "%s could not be downloaded!!!!" % url

########################################################

def get_ligand_pdb_dict(blacklist = False):
    """Returns a pdb code - ligands dictionary"""
    outdict={}
    destfile = os.path.join(CACHEDIR, 'cc-to-pdb.tdd')
    if os.path.isfile(destfile):
        handler = open(destfile, 'r')
        writer = None
    else:
        writer = open(destfile, 'w')
        print(u'Downloading PDBID-HETID dictionary...'),
        dicturl = "http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd"
        handler = urllib2.urlopen(dicturl)
    for line in handler:
        if writer:
            writer.write(line)
        ligand,  pdb_codes = line.strip().split("\t")
        if not blacklist or ligand not in ligand_blacklist:
            for pdb_code in pdb_codes.split():
                if pdb_code not in outdict:
                    outdict[pdb_code] = [ligand, ]
                else:
                    outdict[pdb_code].append(ligand)
    handler.close()
    return outdict

def setglobaldicts():
    global hetdict
    if not hetdict:
        hetdict = get_ligand_pdb_dict()
