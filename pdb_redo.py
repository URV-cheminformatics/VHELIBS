# -*- coding: utf-8 -*-
#
#   Copyright 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import urllib2, sys, os, time, csv, datetime
if sys.platform.startswith('java'):
    import multithreading as multiprocessing
    from java.util import Locale
else:
    import multiprocessing
import PDBfiles

PDB_REDO_ed_data_url_tmpl = "http://www.cmbi.ru.nl/pdb_redo/MIDDLE/PDBID/PDBID_final.eds"
ALLDATA_URL = "http://www.cmbi.ru.nl/pdb_redo/others/alldata.txt"
CACHE_EXPIRED = True

def get_ED_data(pdbid):
    """
    Extract data from EDS site for a given PDB code
    """
    pdbid = pdbid.lower()
    downloaddir = os.path.join(PDBfiles.CACHEDIR, pdbid.lower())
    if not os.path.isdir(downloaddir):
        os.makedirs(downloaddir)
    url = PDB_REDO_ed_data_url_tmpl.replace('MIDDLE', pdbid[1:3]).replace('PDBID', pdbid)
    filename = os.path.join(downloaddir, os.path.basename(url))
    if CACHE_EXPIRED or not os.path.isfile(filename): #Download
        print "Downloading %s" % url
        tries = 3
        while tries > 0:
            tries -= 1
            try:
                rfh = urllib2.urlopen(url)
                ed_file = open(filename, 'w')
                ed_file.write(rfh.read())
                ed_file.close()
                rfh.close()
                break
            except Exception, e:
                print "Could not download",  url
                print e
                print "Retrying...",  tries
                time.sleep(1)
        else:
            print "Unable to download %s" % url
            print e
            return
    pdbdict={pdbid:None}
    edd_dict = {}
    edfile = open(filename, 'r')
    header = []
    for row in csv.reader(edfile, delimiter="\t"):
        if not header:
            header = row
            continue
        preresidue = row[0].split('_')
        if len(preresidue[0]) < 3:
            preresidue[0] = " "*(3-len(preresidue[0]))+preresidue[0]
        preresidue[0] += " "
        if len(preresidue[2]) < 4:
            preresidue[2] = " "*(4-len(preresidue[2]))+preresidue[2]
        residue = ''.join(preresidue)
        if len(residue) < 10:
            residue += " "
        rsr = row[1]
        rscc = row[3]
        #ngrid = row[4]
        edd_dict[residue] = {"RSR":float(rsr) if rsr.strip() else 100
                             ,"RSCC": float(rscc) if rscc.strip() else 0
                             }
    edfile.close()
    return edd_dict

def get_pdbredo_data(pdbids=[]):
    global CACHE_EXPIRED
    alldatapath = os.path.join(PDBfiles.CACHEDIR, os.path.basename(ALLDATA_URL))
    #if not os.path.isfile(alldatapath): #Download
    tries = 3
    while tries > 0:
        tries -= 1
        if sys.platform.startswith('java'):
            oldlocale = Locale.getDefault()
            Locale.setDefault(Locale.ENGLISH)
        try:
            download = False
            olddate = datetime.datetime(1,1,1)
            if os.path.isfile(alldatapath):
                alldata = open(alldatapath, 'r')
                n = 0
                for line in alldata:
                    n +=1
                    if n == 8:
                        olddate = datetime.datetime.strptime(line.split('+')[0], "Created on %a, %d %b %Y %H:%M:%S ")
                alldata.close()
            rfh = urllib2.urlopen(ALLDATA_URL)
            n = 0
            firstlines = ''
            for line in rfh:
                if download:
                    alldata.write(line)
                    continue
                firstlines += line
                n +=1
                if n == 8:
                    newdate = datetime.datetime.strptime(line.split('+')[0], "Created on %a, %d %b %Y %H:%M:%S ")
                    if newdate <= olddate:
                        print "%s is up to date (%s)" % (alldatapath, olddate.__str__())
                        CACHE_EXPIRED = False
                        break
                    else:
                        CACHE_EXPIRED = True
                        download = True
                        alldata = open(alldatapath, 'w')
                        alldata.write(firstlines)
                        print "Downloading %s" % ALLDATA_URL
            #alldata.write(rfh.read())
            alldata.close()
            rfh.close()
            break
        except Exception, e:
            print "Could not download",  ALLDATA_URL
            print e
            print "Retrying...",  tries
            time.sleep(1)
        if sys.platform.startswith('java'):
            Locale.setDefault(oldlocale)
    else:
        print "Unable to download %s" % ALLDATA_URL
        print e
        return
    print "reading alldata.txt"
    alldata = open(alldatapath, 'r')
    #line 7 contains date
    #lines[11:111] are the columns explained
        #We need columns:
        #0   PDBID
        #15  RFFIN
        #62  RESOLUTION
    #lines[112] == '#START DATA'
    #lines[113] are the column headers
    #here the data [114:-3]
    #lines[-3] == ''
    #lines[-2] == '#END DATA'
    parseddata = {}
    data_started = 0
    #pdbids = ('100d', '3dzu')
    header = []
    for row in csv.reader(alldata, delimiter=" ", quotechar="'", quoting = csv.QUOTE_MINIMAL, skipinitialspace = True):
        if row:
            if data_started == 2:
                if row[0] == '#END':
                    break
                elif pdbids and row[0] not in pdbids:
                    continue
                datadict = {}
                #for i in xrange(1, len(row)):
                for i in (15, 62):
                    datadict[header[i]] = row[i]
                parseddata[row[0].upper()] = datadict
            elif data_started == 1:
                data_started = 2
                header = row
            else:
                if row[0] == '#START':
                    data_started = 1
        pass # do stuff
    return parseddata

if __name__ == "__main__":
    print 'start'
    d = get_pdbredo_data()
    print len(d)

