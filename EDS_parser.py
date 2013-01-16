# -*- coding: utf-8 -*-
#
#   Copyright 2010 - 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
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
    pdbid = pdbid.lower()
    pdbdict={pdbid:None}
    edd_dict = {}
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
                    tries = 999
                except Exception, e:
                    if tries >3:
                        print e
                        raise e
                    time.sleep(1)
        if not statfilelines:
            print 'could not read stat file'
            pdbdict[pdbid] = False
            return pdbdict, edd_dict
        for line in statfilelines:
            if not line.startswith('!'):
                residue = line[8:18]
                rscc = line[21:26]
                rsr = line[27:33]
                owab = line[35:40]
#                natom = line[41:46]
#                s_occ = line[47:52]
                edd_dict[residue] = {"RSR":float(rsr) if rsr.strip() else 100
                                     ,"RSCC": float(rscc) if rscc.strip() else 0
                                     ,"OWAB": float(owab) if owab.strip() else 1000
#                                     ,"Natom":float(natom) if natom.strip() else None
#                                     ,"S_occ":float(s_occ) if s_occ.strip() else 0
                                     }
        pdbdict[pdbid] = True
        statfile.close()
    except urllib2.URLError, e:
        if hasattr(e, 'reason'):
            pdbdict[pdbid] = e.reason
        elif hasattr(e, 'code'):
            if e.code == 404:
                pdbdict[pdbid] = False
            else:
                pdbdict[pdbid] = e.code
        else:
            pdbdict[pdbid] = unicode(e)
    return pdbdict, edd_dict
