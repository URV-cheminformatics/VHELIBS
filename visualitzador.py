# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import sys
import os
import csv
import urllib2
import tarfile
from java.awt import Frame, Panel, BorderLayout, FlowLayout, Label, TextField
from java.lang import System
import astex.MoleculeViewer as MoleculeViewer
#print 'astex importat'

import EDS_parser
edsurl = EDS_parser.edsurl

import PDBfiles
PDBfiles.PREFIX = EDS_parser.EDSdir
datadir = PDBfiles.PREFIX

PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"

if not len(sys.argv):
    sys.argv.append('-h')
import rsr_analysis
parser = rsr_analysis.parser
parser.add_argument('-c','--csvfile', metavar='CSVFILE', type=unicode, default=None, required=False, help='CSV file containing results from a previous RSR analysis')
parser.add_argument('--no-view', required=False, action='store_true', help="Do not visualize the generated csv file")

if not os.path.isdir(datadir):
    os.mkdir(datadir)

def reslist_to_sel(reslist):
    sellist = []
    for res in reslist:
        if res.strip():
            try:
                resname = res[:3]
                chain = res[4]
                resnum = res[5:].strip()
                if not resnum.isdigit():
                    for char in resnum:
                        if not char.isdigit():
                            resnum = resnum.replace(char, '')
                #print 'Transformant %s en %s' % (res, "(chain '%s' and residue %s)" % (chain, resnum))
                sellist.append("(chain '%s' and residue %s)" % (chain, resnum))
            except IndexError:
                print "Malformed residue string:"
                print res
    return sellist

def main():
    values = parser.parse_args(sys.argv)
    if not (values.csvfile or values.pdbidfile or values.pdbids or values.swissprot) :
        print "Use the -h or --help options to see how to use this program"
        return
    csvfilename = values.csvfile
    if not csvfilename:
        rsr_analysis.main(values.pdbidfile, pdbidslist = values.pdbids, swissprotlist =values.swissprot , rsr_upper=values.rsr_upper, rsr_lower = values.rsr_lower, distance=values.distance, outputfile = values.outputfile)
        if values.no_view:
            return 0
        else:
            csvfilename = values.outputfile
    print 'Loading data from %s...' % csvfilename,
    resultdict = {}
    if not os.path.isfile(csvfilename):
        print 'File %s does not exist' % csvfilename
        return(1)
    basename = os.path.splitext(os.path.basename(csvfilename))[0]
    if csvfilename.endswith('_wip.csv'):
        basename = basename.replace('_wip', '')
    elif os.path.isfile(basename + '_wip.csv'):
        print "S'ha trobat un fitxer amb les estructures ja mirades desades."
        print "Vols carregar-lo per no haver de començar de nou? (Sí/No)"
        ans = raw_input()
        answered = False
        while not answered:
            if ans.lower().strip() in ('sí', 'si', 'yes', 'ja', 'da', 'bai', 'oui', 'oc', 'òc','jes', 'yeah', 'sim', 'ok', 'Oook', 'y', 's'):
                csvfilename = basename + '_wip.csv'
                answered = True
            elif ans.lower().strip() in ('no', 'non', 'nein', 'nope', 'ez', 'ne', 'pas', 'não', 'nao', 'Eeek', 'n', 'niet'):
                answered = True
            else:
                ans = raw_input()
    csvfile = open(csvfilename, 'rb')
    try:
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
    except csv.Error,  e:
        print e
        dialect = 'excel'
        print 'using default dialect: %s' % dialect
    csvfile.seek(0)
    reader = csv.reader(csvfile, dialect)
    for pdbid, residues_to_exam_string, ligandresidues_string, binding_site_string in reader:
        if len(pdbid) != 4:
            continue
        residues_to_exam = residues_to_exam_string.split(';')
        ligandresidues = ligandresidues_string.split(';')
        binding_site = binding_site_string.split(';')
        resultdict[pdbid] = ligandresidues, residues_to_exam, binding_site
    else:
        print 'Dades carregades'
    if not resultdict:
        print 'Fitxer sense dades!'
        return
    csvfile.close()
    goodfilename = basename + '_good.csv'
    goodfile = open(goodfilename,'ab')
    goodwriter = csv.writer(goodfile, dialect)
    badfilename = basename + '_bad.csv'
    badfile = open(badfilename,'ab')
    badwriter = csv.writer(badfile, dialect)
    dubiousfilename = basename + '_dubious.csv'
    dubiousfile = open(dubiousfilename,'ab')
    dubiouswriter = csv.writer(dubiousfile, dialect)
    writerdict = {'good':goodwriter, 'bad':badwriter, 'dubious':dubiouswriter}
    filesdict = {'good':goodfile, 'bad':badfile, 'dubious':dubiousfile}
    ###
    #Carrega l'astex
    ###
    frame = Frame()
    frame.setLayout(BorderLayout())
    moleculeViewer = MoleculeViewer()
    moleculeViewer.setApplication(False)
    frame.addWindowListener(moleculeViewer)
    frame.setMenuBar(moleculeViewer.createMenuBar())

    # // if the system property arraycopy is set to true
    # // then the renderer should use arraycopy to clear
    # // its buffers rather than brute force
    arraycopy = System.getProperty("arraycopy")
    if arraycopy:
        moleculeViewer.setArrayCopy(True)

    moleculeRenderer =  moleculeViewer.getMoleculeRenderer()

    mvp = Panel()
    mvp.setLayout(BorderLayout())
    mvp.add(moleculeViewer, BorderLayout.CENTER)
    frame.add(mvp, BorderLayout.CENTER)
    frame.pack()
    frame.show()
    moleculeViewer.moleculeRenderer.execute('set selectcount off; set symmetry on;')

    #moleculeViewer.moleculeRenderer.execute('molecule load pdb_files/P00491/2A0W.pdb.gz pdb_files/P00491/2A0W.pdb.gz;')


    ###
    helpmsg = ''
    helpmsg += """> Special commands:
help : print this message
good : Considera l'estructura com a bona, la desa a %s
bad : Considera l'estructura com a incorrecta, la desa a %s
dubious : Considera l'estructura com a dubtosa, la desa a %s
<pdbid> : load this structure from the queue
list : shows the queue of structures
selectbs : select binding site residues
selectligands : select ligand residues
selectresex : select residues to exam from the binding site
reclipmap : re-clips the map to the ligands and residues to exam
Com fer servir l'OpenAstexViewer:
http://openastexviewer.net/web/interface.html """ % (goodfilename, badfilename , dubiousfilename, )
    print helpmsg
    if not os.path.isdir('PDB'):
        os.mkdir('PDB')
    pdbid = None
    needreload = False
    while resultdict:
        inp = ';'
        if not pdbid:
            pdbid = resultdict.keys()[0]
            needreload = True
        if needreload:
            #Netejar
            moleculeViewer.moleculeRenderer.execute('molecule remove *;')
            moleculeViewer.moleculeRenderer.execute('map remove *;')
            ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
            if ligandresidues == ['']:
                print 'Structure without ligands!'
                print 'Skipping it'
                inp = 'dubious'
            else:
                #Descarregar pdb
                localpdb = os.path.join(datadir, pdbid.lower(), pdbid + '.pdb.gz')
                if not os.path.isfile(localpdb):
                    if not os.path.isdir(os.path.join(datadir, pdbid.lower())):
                        os.mkdir(os.path.join(datadir, pdbid.lower()))
                    PDBfiles.get_pdb_file(pdbid.lower(), localpdb)
                relmap = os.path.join(pdbid.lower(), pdbid.lower() + '.omap')
                localmap = os.path.join(datadir, relmap)
                #Visualitzar a l'astex
                moleculeViewer.moleculeRenderer.execute('molecule load %s %s;' % (pdbid, localpdb))
                moleculeViewer.moleculeRenderer.repaint()
                exam_residues_selection = []
                binding_site_selection = []
                ligands_selection = reslist_to_sel(ligandresidues)
                binding_site_selection = reslist_to_sel(binding_site)
                exam_residues_selection = reslist_to_sel(residues_to_exam)
                moleculeViewer.moleculeRenderer.execute('select all;')
                moleculeViewer.moleculeRenderer.execute('display lines off all;')
                moleculeViewer.moleculeRenderer.execute('select none;')
                for bsres in binding_site_selection:
                    moleculeViewer.moleculeRenderer.execute('append %s;' % bsres)
                    moleculeViewer.moleculeRenderer.execute('display lines on %s;' % bsres)
                    moleculeViewer.moleculeRenderer.execute('color white %s;' % bsres)
                    moleculeViewer.moleculeRenderer.repaint()
                for ligandres in ligands_selection:
                    moleculeViewer.moleculeRenderer.execute('append %s;' % ligandres)
                    moleculeViewer.moleculeRenderer.execute('display sticks on %s;' % ligandres)
                    moleculeViewer.moleculeRenderer.execute('color magenta %s;' % ligandres)
                    moleculeViewer.moleculeRenderer.repaint()
                selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
                moleculeRenderer.setCenter(selectedAtoms)
                moleculeViewer.moleculeRenderer.execute('select none;')
                for examres in exam_residues_selection:
                    moleculeViewer.moleculeRenderer.execute('append %s;' % examres)
                    moleculeViewer.moleculeRenderer.execute('display sticks on %s;' % examres)
                    moleculeViewer.moleculeRenderer.execute('color_by_atom;')
                    moleculeViewer.moleculeRenderer.repaint()
                for ligandres in ligands_selection:
                    moleculeViewer.moleculeRenderer.execute('append %s;' % ligandres)
                selectiondict = {'bs':binding_site_selection, 'ligands':ligands_selection, 'resex':exam_residues_selection}

                try:
                    if not os.path.isfile(localmap):
                        genurl = 'http://eds.bmc.uu.se/cgi-bin/eds/gen_zip.pl?pdbCode=' + pdbid.lower()
                        generate = urllib2.urlopen(genurl)
                        generate.read()
                        generate.close()
                        url = edsurl.replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower()).replace('_stat.lis', '.tar.gz')
                        print 'Descarregant %s' % url
                        remotearchive = urllib2.urlopen(url)
                        localarchive = os.path.join(datadir, pdbid.lower(), os.path.basename(url))
                        archivefile = open(localarchive, 'wb')
                        archivefile.write(remotearchive.read())
                        remotearchive.close()
                        archivefile.close()
                        archivefile = tarfile.open(localarchive)
                        archivefile.extract(relmap, datadir)
                        archivefile.close()
                    moleculeViewer.moleculeRenderer.execute('map load %s %s;' % (pdbid, localmap))
                    moleculeViewer.moleculeRenderer.execute('map %s contour 0 yellow;' % (pdbid))
                    #moleculeViewer.moleculeRenderer.execute('select %s or %s or %s;' % exam_residues_selection, ligands_selection, binding_site_selection)
                    selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
                    moleculeRenderer.clipMaps(None, selectedAtoms, True)
                    moleculeViewer.moleculeRenderer.execute('select none;')
                    moleculeViewer.moleculeRenderer.repaint()
                except Exception,  e:
                    print e
                    print 'Impossible carregar el mapa de densitat electronica per %s' % pdbid
                needreload = False
            print "\n####################################"
            print "Viewing structure %s" % pdbid
            print "####################################\n"
            #Descarregar estructura
            #Carregar estructura
            #Fer visible només els lligands i els residus d'interès
        while not inp or inp == ';':
            print ">>>",
            inp = raw_input().strip()
        #inp = 'bad'
        print
        if inp.upper() in resultdict:
            pdbid = inp.upper()
            needreload = True
        elif inp.lower() == 'list':
            print resultdict.keys()
        elif inp.lower() == 'help':
            print helpmsg
        elif inp.lower() in ('good', 'bad', 'dubious'):
            ligandresidues, residues_to_exam, binding_site = resultdict.pop(pdbid)
            writerdict[inp.lower()].writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
            filesdict[inp.lower()].flush()
            ###
            outfile = open(basename + '_wip.csv', 'wb')
            csvfile = csv.writer(outfile)
            csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
            for pdbid in resultdict:
                ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
                csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
            outfile.close()
            pdbid = None
        elif inp.lower() in ('selectbs', 'selectligands', 'selectresex'):
            target = inp.lower().strip().replace('select', '')
            for res in selectiondict[target]:
                moleculeViewer.moleculeRenderer.execute('append %s;' % res)
            moleculeViewer.moleculeRenderer.repaint()
        elif inp.lower().strip() == 'reclipmap':
            moleculeViewer.moleculeRenderer.execute('select none;')
            for res in (selectiondict['ligands'] + selectiondict['resex']):
                moleculeViewer.moleculeRenderer.execute('append %s;' % res)
            selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
            moleculeRenderer.clipMaps(None, selectedAtoms, True)
            moleculeViewer.moleculeRenderer.execute('select none;')
            moleculeViewer.moleculeRenderer.repaint()
        else:
            try:
                if inp[-1] != ';':
                    inp += ';'
                moleculeViewer.moleculeRenderer.execute(inp)
                moleculeViewer.moleculeRenderer.repaint()
            except Exception, e:
                print e
                print "### Ordre o codi PDB no reconeguts: '%s'. Torna-ho a provar ###" % inp
       #Netejar-ho tot
    for file in filesdict.values():
        file.close()
    print "Enhorabona! S'han acabat les estructures!"
