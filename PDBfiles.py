# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2012 Adrián Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Mòdul encarregat de manipular i descarregar fitxers del PDB
"""
import os, shutil, gzip, urllib2, time, sys
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
else:
    import multiprocessing
from cofactors import ligand_blacklist

#Els fitxers .pdb es poden baixar comprimits des de la següent url, substituint %s pel codi del PDB
PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"
#On guardarem els fitxers:
PREFIX = "data"

hetdict = None
inchidict = None

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
                filesize = handler.info().get('Content-Length')
                if filesize:
                    filesize = int(filesize)
                localsize = os.path.getsize(filename)
                download_needed = localsize != filesize
                if not download_needed:
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
            print "Retrying...",  tries
            time.sleep(1)
    if tries > 3:
        print "%s could not be downloaded!!!!" % url

########################################################

def get_ligand_pdb_dict(blacklist = True):
    """Returns a pdb code - ligands dictionary"""
    outdict={}
    print(u'Descarregant diccionari PDBID-HETID...'),
    dicturl = "http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd"
    handler = urllib2.urlopen(dicturl)
    for line in handler:
        ligand,  pdb_codes = line.strip().split("\t")
        if not blacklist or ligand not in ligand_blacklist:
            for pdb_code in pdb_codes.split():
                if not outdict.has_key(pdb_code):
                    outdict[pdb_code] = [ligand, ]
                else:
                    outdict[pdb_code].append(ligand)
    print('fet')
    handler.close()
    return outdict

def get_inchi_list():
    """
    Return {HETID:(InChI, inchikey, NAME),}
    """
    print(u'Descarregant diccionari InChI...'),
    inchilist = urllib2.urlopen('http://ligand-expo.rcsb.org/dictionaries/Components-inchi.ich')
    inchikeylist = urllib2.urlopen('http://ligand-expo.rcsb.org/dictionaries/Components-inchikey.ich')
    inchidict ={}
    for line in inchilist:
        inchikeyline = inchikeylist.readline()
        try:
            inchi, het, name = line.strip().split('\t')
            inchikey, hetk, namek = inchikeyline.strip().split('\t')
            if het != het:
                print 'differenthet_order_error'
                print het
                print hetk
            if het not in ligand_blacklist:
                inchidict[het] = (inchi, inchikey,  name)
        except Exception, e:
#            print e
            continue
    print('fet')
    inchilist.close()
    inchikeylist.close()
    return inchidict

def setglobaldicts():
    """
    """
    global hetdict
    #global inchidict
    if not hetdict:
        hetdict = get_ligand_pdb_dict()
    #inchidict = get_inchi_list()

if __name__ == "__main__":
    print "No"

