# -*- coding: utf-8 -*-
#
#   Copyright 2010, 2011 Adrián Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Mòdul que extreu informació la EDS database
"""
import urllib2, multiprocessing, sys, os, time, PDBfiles
from decimal import Decimal
from Common import logfile

edsurl = "http://eds.bmc.uu.se/eds/dfs/PDB2/PDB1/PDB1_stat.lis"
EDSdir = 'EDS_data'
logfile = logfile % "EDS_error_log"
residuelist = None

if not os.path.isdir(EDSdir):
    os.makedirs(EDSdir)

def get_EDS(pdbid):
    """
    Funció encarregada d'extreure informació d'una pàgina de l'EDS
    """
    pdbdict={'PDB_ID':pdbid.upper(),'IN_EDS':None}
    rsrdict = {}
    statfilepath = os.path.join(EDSdir, '%s_stat.lis' % pdbid.lower())
    try:
        try:
            statfile = open(statfilepath, 'rb')
            statfilelines = statfile.readlines()
        except:
            tries = 0
            url = edsurl.replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower())
            while tries <=3:
                tries += 1
                try:
                    req = urllib2.urlopen(url)
                    statfilelines = req.readlines()
                    statfile = open(statfilepath, 'wb')
                    statfile.writelines(statfilelines)
                    tries =999
                except Exception, e:
                    if tries >3:
                        raise e
                    time.sleep(1)
        for line in statfilelines:
            if not line.startswith('!'):
                rsr = line.split(']')[1].strip().split()[1]
                residue = line.strip().split('[')[1].split(']')[0]
                residue_only = residue.split()[0]
                if not residuelist or residue_only in residuelist:
                    rsrdict[residue] = rsr
        pdbdict['IN_EDS'] = 'TRUE'
        statfile.close()
    except urllib2.URLError, e:
        if hasattr(e, 'reason'):
            pdbdict['IN_EDS'] = e.reason
        elif hasattr(e, 'code'):
            if e.code == 404:
                pdbdict['IN_EDS'] = 'FALSE'
            else:
                pdbdict['IN_EDS'] = e.code
        else:
            pdbdict['IN_EDS'] = unicode(e)
    return pdbdict, rsrdict

def parse_EDS(pdblist, ligandlist = []):
    """
    Obté informació del EDS d'una llista de codis del PDB
    """
    global residuelist
    residuelist = ligandlist
    pdboutdict = {}
    macro_rsrdict = {}
    print "Obtenint informació de l'EDS, això pot trigar una estona..."
    pool = multiprocessing.Pool(multiprocessing.cpu_count()*10 + 1)
    uhlog = open(logfile, "wb")
    uhlog.write('')
    for pdbdict,  rsrdict in pool.map(get_EDS, pdblist):
        if pdbdict != {} and (pdbdict['IN_EDS'] == 'TRUE' or pdbdict['IN_EDS'] == 'FALSE'):
            pdboutdict[pdbdict['PDB_ID']] = pdbdict
        else:
            uhlog.write(pdbdict['PDB_ID'] + " " + str(pdbdict['IN_EDS']) + "\n")
        for ligand in rsrdict:
            macro_rsrdict[pdbdict['PDB_ID'].upper() + ligand.split()[0]] = {
            'PDB_ID':pdbdict['PDB_ID'].upper()
            , 'HET_ID': ligand.split()[0]
            , 'RSR':rsrdict[ligand]
            }
    uhlog.close()
    pool.close()
    pool.join()
    return pdboutdict,  macro_rsrdict

if __name__ == "__main__":
    if len(sys.argv[1:]) != 3:
       exit("""Us:
       %s fitxer_llista_pdbs minim maxim"""  % sys.argv[0])
    file = sys.argv[1]
    maxval = Decimal(sys.argv[3])
    minval = Decimal(sys.argv[2])
    vurl = 'http://eds.bmc.uu.se/cgi-bin/eds/eds_astex.pl?infile=%s'
    resultfile  = open('valors_rsr<%s_%s.csv' % (sys.argv[2] + '~' +sys.argv[3], os.path.basename(file)), 'wb')

    pdblist = open(file, 'rb')
    liganddict = PDBfiles.get_ligand_pdb_dict(blacklist = False)
    print "Llegint fitxers d'entrada"
    for pdbid in pdblist.readlines():
        pdbid = pdbid.lower().strip()
        print 'Mirant informaci de %s' % pdbid
        pdbdict,  rsrdict = get_EDS(pdbid)
        if pdbdict['IN_EDS'] == 'TRUE':
            if liganddict.has_key(pdbid):
                for residue in rsrdict:
                    hetid = residue.split()[0]
                    if hetid in liganddict[pdbid]:
                        if Decimal(rsrdict[residue]) <= maxval and Decimal(rsrdict[residue]) >= minval:
                            line = '"' + '"\t"'.join([pdbid, residue, rsrdict[residue],  vurl % pdbid]) + '"' + '\n'
                            resultfile.write(line)
            else:
                print pdbid, 'no te lligands'
        else:
            print pdbid, "no a l'EDS"
    resultfile.close()
    print 'Fet!'


