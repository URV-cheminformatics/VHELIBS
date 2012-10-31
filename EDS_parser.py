# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
"""
Mòdul que extreu informació la EDS database
"""
import urllib2, sys, os, time
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
else:
    import multiprocessing
import PDBfiles

edsurl = "http://eds.bmc.uu.se/eds/dfs/PDB2/PDB1/PDB1_stat.lis"
residuelist = ''

def get_EDM_sigma(pdbid):
    """
    Downloads the file containing EDM sigma
    """
    edmurl = edsurl.replace('_stat.lis', '.sfdat').replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower())
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    sfdatfile = os.path.join(downloaddir,  '%s.sfdat' % pdbid.lower())
    tries = 0
    while tries <=3:
        try:
            if os.path.isfile(sfdatfile):
                sfdatfilestream =  open(sfdatfile)
                sigma = float([line.split()[1] for line in sfdatfilestream if 'MAP_SIGMA' in line][0])
                sfdatfilestream.close()
                return sigma
            req = urllib2.urlopen(edmurl)
            sfdat = req.read()
            req.close()
            outfile = open(sfdatfile, 'wb')
            outfile.write(sfdat)
            outfile.close()
        except urllib2.URLError, e:
            print e
            tries +=1
    print "Unable to find the map sigma"
    return 1


def get_EDM(pdbid):
    """
    Downloads the EDM of the given pdb code
    """
    edmurl = edsurl.replace('_stat.lis', '.omap').replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower())
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    omapfile = os.path.join(downloaddir,  '%s.omap' % pdbid.lower())
    tries = 0
    while tries <=3:
        try:
            if os.path.isfile(omapfile):
                sigma = get_EDM_sigma(pdbid)
                return omapfile, sigma
            req = urllib2.urlopen(edmurl)
            omap = req.read()
            req.close()
            outfile = open(omapfile, 'wb')
            outfile.write(omap)
            outfile.close()
        except urllib2.URLError, e:
            print e
            tries +=1

def get_EDS(pdbid):
    """
    Extract data from EDS site for a given PDB code
    """
    pdbdict={'PDB_ID':pdbid.upper(),'IN_EDS':None}
    rsrdict = {}
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    statfilepath = os.path.join(downloaddir,  '%s_stat.lis' % pdbid.lower())
    try:
        if os.path.isfile(statfilepath):
            statfile = open(statfilepath, 'rb')
            statfilelines = statfile.readlines()
        else:
            tries = 0
            url = edsurl.replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower())
            print 'Downloading %s' % url
            statfilelines = []
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
                        print e
                        raise e
                    time.sleep(1)
        if not statfilelines:
            print 'could not read stat file'
            pdbdict['IN_EDS'] = 'FALSE'
            return pdbdict, rsrdict
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
    Get EDS data from a list of PDB codes
    """
    global residuelist
    residuelist = ligandlist
    pdboutdict = {}
    macro_rsrdict = {}
    print "Obtenint informació de l'EDS, això pot trigar una estona..."
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    for pdbdict,  rsrdict in pool.map(get_EDS, pdblist):
        if pdbdict != {} and (pdbdict['IN_EDS'] == 'TRUE' or pdbdict['IN_EDS'] == 'FALSE'):
            pdboutdict[pdbdict['PDB_ID']] = pdbdict
        for ligand in rsrdict:
            macro_rsrdict[pdbdict['PDB_ID'].upper() + ligand.split()[0]] = {
            'PDB_ID':pdbdict['PDB_ID'].upper()
            , 'HET_ID': ligand.split()[0]
            , 'RSR':rsrdict[ligand]
            }
    pool.close()
    pool.join()
    return pdboutdict,  macro_rsrdict
