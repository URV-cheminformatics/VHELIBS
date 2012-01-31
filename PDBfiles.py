# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2012 Adrián Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Mòdul encarregat de manipular i descarregar fitxers del PDB
"""
import os, shutil, gzip, multiprocessing, urllib2, time
from cofactors import ligand_blacklist

#Els fitxers .pdb es poden baixar comprimits des de la següent url, substituint %s pel codi del PDB
PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"
#On guardarem els fitxers:
PREFIX = "data"

hetdict = None
inchidict = None

try:
    from config import rapid
except:
    rapid = True


def get_pdb_file(pdbcode, filename = None):
    """
    Es descarrega un fitxer del PDB
    """
    url = PDBbase % pdbcode

    if not filename:
        filename = os.path.split(url)[1]
    tries = 0
    while tries <= 3:
        tries += 1
        try:
            if os.path.isfile(filename):
                if rapid:
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
            tries = 999
        except Exception, e:
            print "No s'ha pogut descarregar",  url
            print "Tornant-ho a intentar...",  tries
            time.sleep(1)
    if tries > 3:
        print "No s'ha descarregat %s !!!!" % url

def get_pdbs_from_SP(pdbrefsdict, download_all_pdb = False):
    """
    Descarrega els fitxers del PDB els codis dels quals estàn continguts en una llista:
    [id, swissprot accn, pdb id] i els ordena per carpetes segons el swissprot accn
    """
    #####################################################

    #####################################################
    print 'Obtenint fitxers del pdb... (això pot trigar una estona)'
    pool = multiprocessing.Pool(multiprocessing.cpu_count()*20 + 1)
    pool.map(download_pdb_to_accn_dir, pdbrefsdict.values())
    pool.close()
    pool.join()


def download_pdb_to_accn_dir(pdbdict):
    """
    Desa un fitxer del pdb sota un directori amb el nom de la seva fitxa corresponent del SwissProt
    """
    saveto = os.path.join(PREFIX,pdbdict['UNIPROT_ACCN'])
    basename = pdbdict['PDB_ID'] + ".pdb.gz"
    if not os.path.isdir(saveto):
        try:
            os.makedirs(saveto)
        except OSError,  e:
            print e
            #Potser s'ha fet des d'un altre fil
    filename = os.path.join(saveto, basename)
    get_pdb_file(pdbdict['PDB_ID'],filename)

########################################################

def parse_query_row(row):
    """
    Funció dedicada a buscar llligands dins dels fitxers pdb a partir de la consulta a la base de dades
    """
    pdb_id, accn = row
    pdbdictlist = []
    liganddictdict = {}
    if hetdict.has_key(pdb_id.lower()):
        for ligand in hetdict[pdb_id.lower()]:
            pdbdict = {'PDB_ID':pdb_id, 'HET_ID': ligand}
            pdbgz = os.path.join(PREFIX, accn, pdb_id + ".pdb.gz")
            if  os.path.isfile(pdbgz):
                pdb_contents = parse_pdb_file(pdbgz)
                if  pdbdict['HET_ID'] not in pdb_contents['SEQRES']:
                    pdbdict['COVALENTLY_BOUND'] = 'FALSE'
                    if pdb_contents.has_key('LINK'):
                        for res1,  res2,  blen in pdb_contents['LINK']:
                            if blen and float(blen) <= 1.9\
                            and (res1 not in ligand_blacklist and res2 not in ligand_blacklist)\
                            and  ((res1 not in pdb_contents['SEQRES'] and res2 in pdb_contents['SEQRES'])\
                            or (res1 in pdb_contents['SEQRES'] and res2 not in pdb_contents['SEQRES'])):
                                pdbdict['COVALENTLY_BOUND'] = 'TRUE'
                    if 'INHIB' in pdb_contents['TITLE'].upper():
                        pdbdict['IS_INHIBITOR'] = 'TRUE'
                    else:
                        pdbdict['IS_INHIBITOR'] = 'FALSE'
                    pdbdictlist.append(pdbdict)
                    if inchidict.has_key(ligand):
                        inchikey = inchidict[ligand][1]
                        inchi = inchidict[ligand][0]
                        hetnam = inchidict[ligand][2]
                        if liganddictdict.has_key(inchikey):
                            liganddict = liganddictdict[inchikey]
                        else:
                            liganddict = {}
                        liganddict['InChI'] = inchi
                        liganddict['HET_ID'] = ligand
                        liganddict['InChIkey'] = inchikey
                        liganddict['HETNAM'] = hetnam
                        if pdb_contents.has_key('HETSYN'):
                            for hetsyn in pdb_contents['HETSYN']:
                                if hetsyn[0] == ligand:
                                    liganddict['HETSYN'] = hetsyn[1]
                        liganddictdict[inchikey] = liganddict
    return pdbdictlist, liganddictdict


def parse_pdb_file(pdbgz):
    """Ha de ser en .gz"""
    pdb_contents = {'TITLE':'', 'SEQRES':set()}
    pdbfile = gzip.GzipFile(pdbgz)
    for line in pdbfile:
        line = line.strip()
        content = None
        if len(line) < 6:
            continue
        label = line[:6].strip()
        if label == 'LINK':
            content = (line[17:20].strip(),  line[47:50].strip(),  line[73:78].strip()) #distancia
        elif label == 'SEQRES':
            content = set(line[19:].split())
        elif label == 'TITLE':
            content = line[6:].strip()
        elif label in ('HETNAM', 'HETSYN'):
            if not line[6:11].strip():
                content = [line[11:14].strip() ,  line[14:].strip()]
            else:
                content = None
                if line[11:14].strip() == pdb_contents[label][-1][0]:
                    pdb_contents[label][-1][1] += line[14:]
                    pdb_contents[label][-1][1] = pdb_contents[label][-1][1].strip()
                    while '  ' in pdb_contents[label][-1][1]:
                        pdb_contents[label][-1][1] = pdb_contents[label][-1][1].replace('  ', ' ')
                    while '- ' in pdb_contents[label][-1][1]:
                        pdb_contents[label][-1][1] = pdb_contents[label][-1][1].replace('- ', '-')
                    while ' -' in pdb_contents[label][-1][1]:
                        pdb_contents[label][-1][1] = pdb_contents[label][-1][1].replace(' -', '-')
                else:
                    print('Error: format inesperat')
        if content:
            if not pdb_contents.has_key(label):#Primer registre LINK
                if label != 'SEQRES':
                    content = [content]
                pdb_contents[label] = content
            else:
                if label == 'SEQRES':
                    pdb_contents[label].update(content)
                elif label == 'TITLE':
                    pdb_contents[label] += ' ' + content
                elif label in ('HETNAM', 'HETSYN', 'LINK'):
                    pdb_contents[label].append(content)
    pdbfile.close()
    return pdb_contents


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
    global inchidict
    hetdict = get_ligand_pdb_dict()
    inchidict = get_inchi_list()

if __name__ == "__main__":
    print "No"

