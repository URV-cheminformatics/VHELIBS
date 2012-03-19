# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#

#Python stuff
import sys
import os
import csv
import urllib2
import tarfile

#Java stuff
from java.lang.System import exit

from java.awt import BorderLayout, Container, Dimension, Graphics
from java.awt.event import WindowAdapter, WindowEvent

from javax.swing import JFrame, JPanel

#Jmol stuff
from org.jmol.adapter.smarter import SmarterJmolAdapter
from org.jmol.api import JmolViewer
from org.jmol.util import Logger
from org.openscience.jmol.app.jmolpanel import AppConsole

#Own stuff
import EDS_parser
edsurl = EDS_parser.edsurl #Deprecate

import PDBfiles
PDBfiles.PREFIX = EDS_parser.EDSdir
datadir = PDBfiles.PREFIX

PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz" #Deprecate

if not len(sys.argv):
    sys.argv.append('-h')
import rsr_analysis

###Define useful classes###
class ApplicationCloser(WindowAdapter):
    def windowClosing(self, e):
        exit(0)

class JmolPanel(JPanel):
    def __init__(self):
        self.currentSize = Dimension()
        self.viewer = JmolViewer.allocateViewer(self, SmarterJmolAdapter(), None, None, None, None, None)

    def paint(self, g):
        self.getSize(self.currentSize)
        self.viewer.renderScreenImage(g, self.currentSize.width, self.currentSize.height)

if not os.path.isdir(datadir): #Deprecate
    os.mkdir(datadir)

yes = ('sí', 'si', 'yes', 'ja', 'da', 'bai', 'oui', 'oc', 'òc','jes', 'yeah', 'sim', 'ok', 'oook', 'y', 's')
no = ('no', 'non', 'nein', 'nope', 'ez', 'ne', 'pas', 'não', 'nao', 'eeek', 'n', 'niet')

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
                sellist.append("[%s]%s:%s" % (resname, resnum, chain))
            except IndexError:
                print "Malformed residue string:"
                print res
    if sellist:
        return ' or '.join(sellist)
    else:
        return 'none'

def main():
    """
    """
    ### build the parser###
    print sys.argv
    parser = rsr_analysis.parser
    parser.add_argument('-c','--csvfile', metavar='CSVFILE', type=unicode, default=None, required=False, help='CSV file containing results from a previous RSR analysis')
    parser.add_argument('--no-view', required=False, action='store_true', help="Do not visualize the generated csv file")
    values = parser.parse_args(sys.argv)
    if not (values.csvfile or values.pdbidfile or values.pdbids or values.swissprot) :
        print "Use the -h or --help options to see how to use this program"
        return
    csvfilename = values.csvfile
    if not csvfilename:
        datawritten = rsr_analysis.main(values.pdbidfile, pdbidslist = values.pdbids, swissprotlist =values.swissprot , rsr_upper=values.rsr_upper, rsr_lower = values.rsr_lower, distance=values.distance, outputfile = values.outputfile)
        if values.no_view:
            exit(0)
        if datawritten:
            csvfilename = values.outputfile
        else:
            print 'No structures to be viewed.'
            exit(0)
    print 'Loading data from %s...' % csvfilename,
    resultdict = {}
    if not os.path.isfile(csvfilename):
        print 'File %s does not exist' % csvfilename
        exit(1)
    outdir = os.path.dirname(csvfilename)
    basename = os.path.splitext(os.path.basename(csvfilename))[0]
    if csvfilename.endswith('_wip.csv'):
        basename = basename.replace('_wip', '')
    elif os.path.isfile(os.path.join(outdir, basename + '_wip.csv')):
        print "S'ha trobat un fitxer amb les estructures ja mirades desades."
        print "Vols carregar-lo per no haver de començar de nou? (Sí/No)"
        ans = raw_input()
        answered = False
        while not answered:
            if ans.lower().strip() in yes:
                csvfilename = os.path.join(outdir, basename + '_wip.csv')
                answered = True
            elif ans.lower().strip() in no:
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
        exit(1)
    csvfile.close()
    goodfilename = os.path.join(outdir, basename + '_good.csv')
    goodfile = open(goodfilename,'ab')
    goodwriter = csv.writer(goodfile, dialect)
    badfilename = os.path.join(outdir, basename + '_bad.csv')
    badfile = open(badfilename,'ab')
    badwriter = csv.writer(badfile, dialect)
    dubiousfilename = os.path.join(outdir, basename + '_dubious.csv')
    dubiousfile = open(dubiousfilename,'ab')
    dubiouswriter = csv.writer(dubiousfile, dialect)
    writerdict = {'good':goodwriter, 'bad':badwriter, 'dubious':dubiouswriter}
    filesdict = {'good':goodfile, 'bad':badfile, 'dubious':dubiousfile}
    ###                    ###
    ###Load the GUI###
    ###                    ###
    frame = JFrame("Structure Validation Helper")
    frame.addWindowListener(ApplicationCloser())
    frame.setSize(410, 700);
    contentPane = frame.getContentPane()
    jmolPanel = JmolPanel()
    jmolPanel.setPreferredSize(Dimension(400, 400))
    # main panel -- Jmol panel on top
    panel = JPanel()
    panel.setLayout(BorderLayout())
    panel.add(jmolPanel)
    # main panel -- console panel on bottom
    panel2 = JPanel()
    panel2.setLayout(BorderLayout())
    panel2.setPreferredSize(Dimension(400, 300))
    console = AppConsole(jmolPanel.viewer, panel2,"History State Clear")
    jmolPanel.viewer.setJmolCallbackListener(console)
    panel.add("South", panel2)
    contentPane.add(panel)
    jmolPanel.viewer.evalString('wireframe off')
    jmolPanel.viewer.evalString('ribbon off')
    jmolPanel.viewer.evalString('meshribbon off')
    jmolPanel.viewer.evalString('meshribbon off')
    jmolPanel.viewer.evalString('cartoon off')
    jmolPanel.viewer.evalString('set bondMode OR')
    frame.setVisible(True)
#        Logger.error(strError)

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
exit: Exits the program
Com fer servir l'OpenAstexViewer:
http://openastexviewer.net/web/interface.html """ % (goodfilename, badfilename , dubiousfilename, )
    print helpmsg
    pdbid = None
    needreload = False
    wipfilename = os.path.join(outdir, basename + '_wip.csv')
    while resultdict:
        inp = ';'
        if not pdbid:
            pdbid = resultdict.keys()[0]
            needreload = True
        if needreload:
            #Netejar
            jmolPanel.viewer.evalString('delete')
            ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
            if ligandresidues == ['']:
                print 'Structure without ligands!'
                print 'Skipping it'
                inp = 'dubious'
            else:
                #load in Jmol
                jmolPanel.viewer.evalString('load "=%s"' % pdbid)
                ligands_selection = reslist_to_sel(ligandresidues)
                binding_site_selection = reslist_to_sel(binding_site)
                exam_residues_selection = reslist_to_sel(residues_to_exam)
                jmolPanel.viewer.evalString('select all')
                jmolPanel.viewer.evalString('wireframe off')
                jmolPanel.viewer.evalString('select none')

                jmolPanel.viewer.evalString('select(%s)' % binding_site_selection)
                jmolPanel.viewer.evalString('wireframe 0.01')
                jmolPanel.viewer.evalString('spacefill off')
                jmolPanel.viewer.evalString('color white')
                jmolPanel.viewer.evalString('select none')

                jmolPanel.viewer.evalString('select(%s)' % ligands_selection)
                jmolPanel.viewer.evalString('wireframe 0.1')
                jmolPanel.viewer.evalString('spacefill 0.2')
                jmolPanel.viewer.evalString('color magenta')
                jmolPanel.viewer.evalString('select none')

                jmolPanel.viewer.evalString('center {%s}' % ligands_selection)

                jmolPanel.viewer.evalString('select(%s)' % exam_residues_selection)
                jmolPanel.viewer.evalString('wireframe 0.1')
                jmolPanel.viewer.evalString('spacefill 0.2')
                jmolPanel.viewer.evalString('color cpk')
                jmolPanel.viewer.evalString('select none')

                selectiondict = {'bs':binding_site_selection, 'ligands':ligands_selection, 'resex':exam_residues_selection}

                try:
                    jmolPanel.viewer.evalString('isosurface s2 color yellow cutoff 0.33 within 2.0 {%s} "=%s" mesh nofill' %  (' or '.join([ligands_selection, exam_residues_selection]), pdbid ))
                except Exception,  e:
                    print e
                    print 'Impossible carregar el mapa de densitat electronica per a %s' % pdbid
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
            outfile = open(wipfilename, 'wb')
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
                jmolPanel.viewer.evalString('append %s;' % res)
        elif inp.lower().strip() == 'reclipmap':
            jmolPanel.viewer.evalString('select none;')
            for res in (selectiondict['ligands'] + selectiondict['resex']):
                jmolPanel.viewer.evalString('append %s;' % res)
#            selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
#            moleculeRenderer.clipMaps(None, selectedAtoms, True)
            jmolPanel.viewer.evalString('select none;')
        elif inp.lower().strip() == 'exit':
            exit(0)
        else:
            try:
                if inp[-1] != ';':
                    inp += ';'
                jmolPanel.viewer.evalString(inp)
            except Exception, e:
                print e
                print "### Ordre o codi PDB no reconeguts: '%s'. Torna-ho a provar ###" % inp
       #Netejar-ho tot
    for file in filesdict.values():
        file.close()
    try:
        os.remove(wipfilename)
    except:
        pass
    jmolPanel.viewer.evalString('delete;')
    print "No more structures to view!"
    print 'Do you want to exit now?'
    ans = raw_input().lower().strip()
    if ans in yes:
        exit(0)
    else:
        print 'you can still manually load more structures and do other operations through OpenAstexViewer scripting'

if __name__ == '__main__':
    sys.argv = [arg for arg in sys.argv if __file__ not in arg]
    main()
