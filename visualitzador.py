# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#

#Python stuff
import sys
import os
import csv

#Java stuff
from java.lang.System import exit

from java.awt import BorderLayout, Dimension, GridLayout

from javax.swing import JFrame, JPanel, JButton

#Jmol stuff
from org.jmol.adapter.smarter import SmarterJmolAdapter
from org.jmol.api import JmolViewer
from org.jmol.util import Logger
from org.openscience.jmol.app.jmolpanel import AppConsole

#Own stuff
if not len(sys.argv):
    sys.argv.append('-h')
import rsr_analysis

yes = ('sí', 'si', 'yes', 'ja', 'da', 'bai', 'oui', 'oc', 'òc','jes', 'yeah', 'sim', 'ok', 'oook', 'y', 's')
no = ('no', 'non', 'nein', 'nope', 'ez', 'ne', 'pas', 'não', 'nao', 'eeek', 'n', 'niet')

###Define useful classes###

class JmolPanel(JPanel):
    def __init__(self, preferredSize):
        self.currentSize = Dimension()
        self.preferredSize = preferredSize
        self.viewer = JmolViewer.allocateViewer(self, SmarterJmolAdapter(), None, None, None, None, None)

    def paint(self, g):
        self.getSize(self.currentSize)
        self.viewer.renderScreenImage(g, self.currentSize.width, self.currentSize.height)

class StruVa():
    def __init__(self, csvfilename=None):
        if csvfilename:
            self.loadCSV(csvfilename)
            self.setupUi()
            self.setVisible(True)
            self.console.sendConsoleMessage(self.helpmsg)
            self.start()

    def setupUi(self):
        frame = JFrame("Structure Validation Helper", defaultCloseOperation = JFrame.EXIT_ON_CLOSE, size = (700, 410))
        contentPane = frame.contentPane
        jmolPanel = JmolPanel(preferredSize = Dimension(400, 400))
        self.viewer = jmolPanel.viewer
        self.execute = self.viewer.evalStringQuiet
        # main panel -- Jmol panel on left
        panel = JPanel(layout = BorderLayout())
        panel.add(jmolPanel)
        # main panel -- console panel on right
        panel2 = JPanel(layout = BorderLayout(), preferredSize = Dimension(300, 400))
        self.console = AppConsole(jmolPanel.viewer, panel2,"History Variables State Clear Help")
        jmolPanel.viewer.jmolCallbackListener = self.console
        panel.add("East", panel2)
        contentPane.add(panel)
        self.execute('wireframe only')
        self.execute('wireframe off')
        self.execute('set bondMode OR')
        buttonPanel = JPanel(GridLayout(2, 2))
        self.good_button = JButton('Good', actionPerformed=self.nextStruct)
        self.bad_button = JButton('Bad', actionPerformed=self.nextStruct)
        self.dubious_button = JButton('Dubious', actionPerformed=self.nextStruct)
        self.skip_button = JButton('Skip', actionPerformed=self.nextStruct)
        buttonPanel.add(self.good_button)
        buttonPanel.add(self.bad_button)
        buttonPanel.add(self.dubious_button)
        buttonPanel.add(self.skip_button)
        panel.add(buttonPanel, BorderLayout.NORTH)
        self.setVisible = frame.setVisible

    def nextStruct(self, event = None, text = None):
        if event:
            print event
            print event.source == self.skip_button
            print event.source.text
            text = event.source.text
        if not text:
            return


    def reloadStruct(self, pdbid):
        #Netejar
        self.execute('delete')
        ligandresidues, residues_to_exam, binding_site = self.resultdict[pdbid]
        if ligandresidues == ['']:
            self.console.sendConsoleMessage( 'Structure without ligands!')
            self.console.sendConsoleMessage( 'Skipping it')
            self.nextStruct(text='Skip')
        else:
            #load in Jmol
            try:
                self.execute('load "=%s"' % pdbid)
                self.execute('select all')
                self.execute('wireframe only')
                self.execute('wireframe off')
                self.execute('select none')
                self.displayBindingSite(binding_site)
                self.displayResToExam(residues_to_exam)
                self.displayLigands(ligandresidues)
                self.loadEDM(pdbid)
            except Exception,  e:
                self.console.sendConsoleMessage("ERROR: " + e)
        self.console.sendConsoleEcho( "\n####################################")
        self.console.sendConsoleEcho( "Viewing structure %s" % pdbid)
        self.console.sendConsoleEcho( "####################################\n")

    def displayBindingSite(self, binding_site):
        binding_site_selection = reslist_to_sel(binding_site)
        self.execute('define binding_site (%s)' % binding_site_selection)
        self.execute('select(binding_site)')
        self.execute('wireframe 0.01')
        self.execute('spacefill off')
        self.execute('color white')
        self.execute('select none')

    def displayResToExam(self, residues_to_exam):
        exam_residues_selection = reslist_to_sel(residues_to_exam)
        self.execute('define res_to_exam (%s)' % exam_residues_selection)
        self.execute('select(res_to_exam)' )
        self.execute('wireframe 0.1')
        self.execute('spacefill 0.2')
        self.execute('color cpk')
        self.execute('select none')

    def displayLigands(self, ligandresidues):
        ligands_selection = reslist_to_sel(ligandresidues)
        self.execute('define ligands (%s)' % ligands_selection)
        self.execute('select(ligands)')
        self.execute('wireframe 0.1')
        self.execute('spacefill 0.2')
        self.execute('color magenta')
        self.execute('select none')
        self.execute('center ligands')

    def loadEDM(self, pdbid):
        self.execute('isosurface EDM color yellow sigma 1.0 within 2.0 {ligands or res_to_exam} "=%s" mesh nofill' %  pdbid)

    def saveWIP(self):
        outfile = open(self.wipfilename, 'wb')
        csvfile = csv.writer(outfile)
        csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
        for pdbid in self.resultdict:
            ligandresidues, residues_to_exam, binding_site = self.resultdict[pdbid]
            csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
        outfile.close()

    def updateOutFile(self, filetype, pdbid):
        filetype = filetype.lower()
        if filetype not in self.writerdict.keys():
            raise TypeError('Unknown destination file')
        ligandresidues, residues_to_exam, binding_site = self.resultdict.pop(pdbid)
        self.writerdict[filetype].writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
        self.filesdict[filetype].flush()

    def start(self):
        pdbid = None
        needreload = False
        while self.resultdict:
            inp = ';'
            if not pdbid:
                pdbid = self.resultdict.keys()[0]
                needreload = True
            if needreload:
                self.reloadStruct(pdbid)
                needreload = False
            while not inp or inp == ';':
                self.console.sendConsoleEcho( ">>>")
                inp = raw_input().strip()
            print
            if inp.upper() in self.resultdict:
                pdbid = inp.upper()
                needreload = True
            elif inp.lower() == 'list':
                self.console.sendConsoleEcho( self.resultdict.keys())
            elif inp.lower() == 'help':
                self.console.sendConsoleMessage( self.helpmsg)
            elif inp.lower() in ('good', 'bad', 'dubious'):
                self.updateOutFile(inp, pdbid)
                ###
                self.saveWIP()
                pdbid = None
            elif inp.lower().strip() == 'exit':
                exit(0)
            else:
                try:
                    if inp[-1] != ';':
                        inp += ';'
                    self.viewer.evalString(inp)
                except Exception, e:
                    self.console.sendConsoleMessage("ERROR: " + e)
                    self.console.sendConsoleMessage( "### Ordre o codi PDB no reconeguts: '%s'. Torna-ho a provar ###" % inp)
           #Netejar-ho tot
        for file in self.filesdict.values():
            file.close()
        try:
            os.remove(self.wipfilename)
        except:
            pass
        self.execute('delete;')
        self.console.sendConsoleMessage( "No more structures to view!")
        self.console.sendConsoleMessage( 'Do you want to exit now?')
        ans = raw_input().lower().strip()
        if ans in yes:
            exit(0)
        else:
            self.console.sendConsoleMessage( 'you can still manually load more structures and do other operations through OpenAstexViewer scripting')

    def loadCSV(self, csvfilename):
        self.resultdict = {}
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
            self.resultdict[pdbid] = ligandresidues, residues_to_exam, binding_site
        else:
            print 'Dades carregades'
        if not self.resultdict:
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
        self.writerdict = {'good':goodwriter, 'bad':badwriter, 'dubious':dubiouswriter}
        self.filesdict = {'good':goodfile, 'bad':badfile, 'dubious':dubiousfile}
        self.wipfilename = os.path.join(outdir, basename + '_wip.csv')
        self.helpmsg = """> Special commands:
help : print this message
good : Considera l'estructura com a bona, la desa a %s
bad : Considera l'estructura com a incorrecta, la desa a %s
dubious : Considera l'estructura com a dubtosa, la desa a %s
see <pdbid> : load this structure from the queue
list : shows the queue of structures
binding_site : select binding site residues
ligands : select ligand residues
res_to_exam : select residues to exam from the binding site
Jmol scripting manual:
http://chemapps.stolaf.edu/jmol/docs/?&fullmanual=1&ver=12.4 """ % (goodfilename, badfilename , dubiousfilename, )

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
    struva = StruVa(csvfilename)

if __name__ == '__main__':
    sys.argv = [arg for arg in sys.argv if __file__ not in arg]
    main()
