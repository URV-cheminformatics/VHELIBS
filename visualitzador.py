# -*- coding: utf-8 -*-
#
#   Copyright 2012, 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#

#Python stuff
######One-jar magic#######
import sys
sys.path.append('__pyclasspath__/pylib')
if not sys.prefix:
    sys.prefix='.'
######One-jar magic#######
import shutil, tempfile, os

import csv
import math
import time
from textwrap import dedent
from sys import exit
try:
    import cPickle as pickle
except ImportError:
    import pickle
#Java stuff
import java
from java.net import URL, URI
from java.lang import Runnable
from java.awt import BorderLayout, Dimension, GridLayout, GridBagLayout, GridBagConstraints, Insets, Cursor, Desktop
from java.awt.event import ItemEvent, ActionListener, WindowAdapter
from javax.swing import JFrame, JPanel, JButton, JOptionPane, JTextField, JCheckBox, JLabel, UIManager, JDialog, SwingUtilities, SwingWorker, JComboBox, ToolTipManager, ImageIcon

systemlaf = UIManager.getSystemLookAndFeelClassName()
_infoicon = UIManager.getIcon("OptionPane.informationIcon")
if systemlaf == u'javax.swing.plaf.metal.MetalLookAndFeel':
    try:
        import com.sun.java.swing.plaf.gtk.GTKLookAndFeel
        UIManager.setLookAndFeel(u'com.sun.java.swing.plaf.gtk.GTKLookAndFeel')
    except Exception,  e:
        print e
else:
    UIManager.setLookAndFeel(systemlaf)

r = java.lang.ClassLoader.getSystemClassLoader().getResource('icon.png')
if r:
    u = URL(str(r))
    vhelibsicon = ImageIcon(u).image
else:
    vhelibsicon=ImageIcon('icon.png').image

#Jython-specific stuff
from swingutils.preferences import getUserPrefs
from swingutils.dialogs.filechooser import showOpenDialog, SimpleFileFilter
from swingutils.dialogs.basic import showErrorDialog, showWarningDialog, showMessageDialog
prefs = getUserPrefs('struva')

#Jmol stuff
from org.jmol.adapter.smarter import SmarterJmolAdapter
from org.jmol.api import JmolViewer
from org.openscience.jmol.app.jmolpanel import AppConsole

VHELIBS_VERSION = "3.0.5"
TITLE =  "VHELIBS " + VHELIBS_VERSION

#Own stuff
sys.argv = [arg for arg in sys.argv if __file__ not in arg]
if not len(sys.argv):
    sys.argv.append('--no-args')
import rsr_analysis, PDBfiles,  EDS_parser, pdb_redo
### build the parser###
argparser = rsr_analysis.parser
argparser.add_argument('-c','--csvfile', metavar='CSVFILE', type=unicode, default=None, required=False, help='CSV file containing results from a previous RSR analysis')
argparser.add_argument('--no-view', required=False, action='store_true', help="Do not visualize the generated csv file")
argparser.add_argument('--no-args', required=False, action='store_true')

def prefbool(string):
    if type(string) == type(True):
        return string
    elif type(string) == type(None):
        return True
    if string.lower() == u'true':
        return True
    elif string.lower() == u'false':
        return False
    else:
        raise TypeError(string + ' cannot be booleaned!')

def load_model(pdbid, source):
    if source == 'PDB_REDO':
        return PDBfiles.get_pdb_file(pdbid, True)
    if source != 'PDB':
        print "WARNING: Unknown model source: %s"% source
        print "Loading from the PDB instead"
    return PDBfiles.get_pdb_file(pdbid, False)

###Define useful classes###
class Console(AppConsole):
    def __init__(self, viewer, panel):
        AppConsole.__init__(self, viewer, panel, '')
        self.pane.components[0].remove(1)

class JmolPanel(JPanel):
    def __init__(self, preferredSize):
        self.currentSize = Dimension()
        self.preferredSize = preferredSize
        self.viewer = JmolViewer.allocateViewer(self, SmarterJmolAdapter(), None, None, None, None, None)

    def paint(self, g):
        self.getSize(self.currentSize)
        self.viewer.renderScreenImage(g, self.currentSize.width, self.currentSize.height)

class make_listen(ActionListener):
    def __init__(self, callable):
        self.actionPerformed = callable

class StruVa(Runnable):
    actions = (u'toggle ligand', u'toggle binding site', u'toggle coordinates to exam',)
    frame = None
    viewer = None
    def __init__(self, values):
        self.values = values
        self.wd = WaitDialog()
        if self.values:
            try:
                restart = self.loadCSV()
                if restart:
                    self.restart()
                    return
                if not self.viewer:
                    SwingUtilities.invokeAndWait(self)
                else:
                    self.run()
                self.start()
            except Exception, e:
                self.e = e
                print e
                showErrorDialog('Unable to load RSR analysis results file:\n %s' % str(e))

    def run(self):
        if self.values:
            self.setupUi()
            self.frame.setVisible(True)

    def setupUi(self):
        self.actionsDict = {}
        self.frame = JFrame(TITLE, iconImage=vhelibsicon, defaultCloseOperation = JFrame.EXIT_ON_CLOSE, size = (700, 410))
        jmolPanel = JmolPanel(preferredSize = (500, 500))
        self.viewer = jmolPanel.viewer
        self.execute = self.viewer.evalStringQuiet
        # main panel -- Jmol panel on left
        panel = JPanel(layout = GridBagLayout())
        panelc = GridBagConstraints()
        panelc.gridwidth = 2
        panelc.gridx = 0
        panelc.gridy = 1
        panelc.weightx = 1
        panelc.weighty = 1
        panelc.fill = GridBagConstraints.BOTH
        panel.add(jmolPanel, panelc)
        # main panel -- console panel on right
        panel2 = JPanel(layout = BorderLayout(), preferredSize = (300, 500))
        self.console = Console(jmolPanel.viewer, panel2)
        jmolPanel.viewer.jmolCallbackListener = self.console
        self.viewer=jmolPanel.viewer
        panelc.gridwidth = 1
        panelc.gridheight = 2
        panelc.gridx = 2
        panelc.gridy = 0
        panelc.weightx = 0
        panelc.weighty = 0.5
        panelc.fill = GridBagConstraints.VERTICAL
        panel.add(panel2, panelc)
        self.frame.add(panel)
        self.execute(';'.join(['wireframe only'
        ,'wireframe off'
        ,'set frank off'
        ,'set bondMode OR'
        ,'set syncScript ON'
        ,'set antialiasDisplay ON'
        ,'set antialiasTranslucent ON'
        ]))
        constraints = GridBagConstraints()
        constraints.gridwidth = 1
        constraints.gridheight =1
        constraints.ipadx =0
        constraints.ipady =0
        constraints.weightx =0.5
        constraints.weighty =0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        cbpanel = JPanel(GridLayout(2, 3, 3, 3))
        later = []
        for action in self.actions:
            caction = action[0].upper() + action[1:]
            if 'toggle' in action:
                cbpanel.add(JCheckBox(caction, itemStateChanged=self.updateDisplay))
                if 'ligand' in action:
                    checked = prefs.get('ligand_edm', False)
                elif 'binding' in action:
                    checked = prefs.get('bindingsite_edm', False)
                elif 'exam' in action:
                    checked = prefs.get('coordstoexam_edm', True)
                later.append(JCheckBox(caction.replace('Toggle', 'EDM for'), prefbool(checked), itemStateChanged=self.updateDisplay))
        #Must first add the checkboxes to the panel, then refrence them
        #Otherwise their selection state is not correctly accessed
        for cb in later:
            cbpanel.add(cb)
        for comp in cbpanel.components:
            self.actionsDict[comp.text.lower()] = comp
        panelc.gridwidth = 2
        panelc.gridheight = 1
        panelc.gridx = 0
        panelc.gridy = 0
        panelc.weightx = 1
        panelc.weighty = 0
        panel.add(cbpanel, panelc)

        buttonPanel = JPanel(GridBagLayout())

        idlist = self.resultdict.keys()
        idlist.sort()
        self.structs_cbox = JComboBox(idlist)
        self.structs_cbox.addActionListener(make_listen(self.reloadStruct))

        about_button = JButton('About', actionPerformed=self.displayAbout)
        opt_button = JButton('Display settings', actionPerformed=self.showDisplaySettings)
        next_button = JButton('Save & Check Next Structure', toolTipText="Save the current structure with the selected flags to %s, then load another structure to check" % self.checkedfilename, actionPerformed=self.nextStruct)

        struct_lbl = JLabel('Current Structure:')
        lig_lbl = JLabel('Ligand')
        bs_lbl = JLabel('Binding Site')

        qualities = ('Good', 'Bad', 'Dubious')

        self.lig_cbox = JComboBox(qualities)
        self.bs_cbox = JComboBox(qualities)

        constraints = GridBagConstraints()
        constraints.gridwidth = 1
        constraints.gridheight =1
        constraints.ipadx =0
        constraints.ipady =0
        constraints.weightx =0.5
        constraints.weighty =0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        constraints.gridx = 0
        constraints.gridy = 0

        buttonPanel.add(about_button, constraints)
        constraints.gridx += 1
        buttonPanel.add(opt_button, constraints)
        constraints.gridx -= 1
        constraints.gridy += 1

        buttonPanel.add(struct_lbl, constraints)
        constraints.gridx += 1
        buttonPanel.add(self.structs_cbox, constraints)
        constraints.gridx -= 1
        constraints.gridy += 1

        buttonPanel.add(lig_lbl, constraints)
        constraints.gridx += 1
        buttonPanel.add(self.lig_cbox, constraints)
        constraints.gridx -= 1
        constraints.gridy += 1
        buttonPanel.add(bs_lbl, constraints)
        constraints.gridx += 1
        buttonPanel.add(self.bs_cbox, constraints)
        constraints.gridx -= 1
        constraints.gridy += 1
        constraints.gridwidth = 2
        buttonPanel.add(next_button, constraints)
        self.buttonPanel = buttonPanel

        panel2.add(self.buttonPanel, BorderLayout.NORTH)
        self.panel = panel
        self.frame.pack()
        self.optionsdiag = DisplaySettingsDialog(self)
        self.aboutdiag = AboutDialog(self)

    def updateDisplay(self, event):
        text = event.source.text.lower()
        ltext = text.lower()
        if 'toggle' in ltext:
            checked = event.getStateChange() == ItemEvent.SELECTED
            if 'ligand' in ltext:
                prefs['ligand'] = checked
                self.displayLigand(checked)
            elif 'binding' in ltext or 'site' in ltext:
                prefs['bindingsite'] = checked
                self.displayBindingSite(checked)
            elif 'exam' in ltext:
                prefs['coordstoexam'] = checked
                self.displayCoordsToExam(checked)
        elif 'edm for' in ltext:
            checked = event.getStateChange() == ItemEvent.SELECTED
            if 'ligand' in ltext:
                prefs['ligand_edm'] = checked
                self.displayLigand(self.actionsDict[u'toggle ligand'].selected)
            elif 'binding' in ltext or 'site' in ltext:
                prefs['bindingsite_edm'] = checked
                self.displayBindingSite(self.actionsDict[u'toggle binding site'].selected)
            elif 'exam' in ltext:
                prefs['coordstoexam_edm'] = checked
                self.displayCoordsToExam(self.actionsDict[u'toggle coordinates to exam'].selected)

    def nextStruct(self, event):
        if self.resultdict:
            bs_valid = self.bs_cbox.selectedItem.lower()
            ligand_valid = self.lig_cbox.selectedItem.lower()
            if bs_valid == 'good':
                self.resultdict[self.key][4] = True
            elif bs_valid == 'bad':
                self.resultdict[self.key][4] = False
            else:
                self.resultdict[self.key][4] = bs_valid

            if ligand_valid == 'good':
                self.resultdict[self.key][3] = True
            elif ligand_valid == 'bad':
                self.resultdict[self.key][3] = False
            else:
                self.resultdict[self.key][3] = ligand_valid
            self.updateOutFile()
            self.structs_cbox.removeItem(self.key)
            if not self.resultdict:
                self.key = None
                self.pdbid = None
                self.clean()
        else:
            self.clean()

    def showDisplaySettings(self, event):
        self.optionsdiag.show(not self.optionsdiag.diag.visible)
        self.optionsdiag.diag.pack()

    def displayAbout(self, event):
        self.aboutdiag.show()
        pass

    def clean(self):
        #Clean everything
        self.execute('delete;')
        game_over = JOptionPane.showConfirmDialog(self.frame,u'Continue working with other structures?',u'No more structures to view!',JOptionPane.YES_NO_OPTION)
        if game_over == JOptionPane.OK_OPTION:
            self.restart()
            return
        else:
            exit(0)

    def reloadStruct(self, event=None):
        if event and event.actionCommand != u'comboBoxChanged': return
        self.console.clearContent('')
        self.key = self.structs_cbox.selectedItem
        if not self.key: return
        self.pdbid = self.key.split('|')[0]
        self.wd.show(False)
        self.wd = WaitDialog(parent=self.frame,info='<html>Loading structure %s from PDB and EDS<br /> Please be patient</html>' % self.key)
        self.wd.dialog.pack()
        self.ds = DialogShower(self.wd, self.viewer)
        self.ds.execute()
        showbs = prefs.get('bindingsite', True)
        showre = prefs.get('coordstoexam', True)
        showlig = prefs.get('ligand', True)
        self.actionsDict[u'toggle binding site'].selected = False
        self.actionsDict[u'toggle coordinates to exam'].selected = False
        self.actionsDict[u'toggle ligand'].selected = False
        prefs['bindingsite'] = showbs
        prefs['coordstoexam'] = showre
        prefs['ligand'] = showlig
        #Clean up
        self.execute('delete')
        #Load relevant data
        self.ligandresidues, self.residues_to_exam, self.binding_site,  self.ligandgood, self.bsgood, self.source = self.resultdict[self.key]
        self.ligandresidues_IS = self.residues_to_exam_IS = self.binding_site_IS = None
        if self.ligandgood.lower() == 'true':
            self.lig_cbox.selectedItem = 'Good'
        elif self.ligandgood.lower() == 'false':
            self.lig_cbox.selectedItem = 'Bad'
        else:
            self.lig_cbox.selectedItem = 'Dubious'
        if self.bsgood.lower() == 'true':
            self.bs_cbox.selectedItem = 'Good'
        elif self.bsgood.lower() == 'false':
            self.bs_cbox.selectedItem = 'Bad'
        else:
            self.bs_cbox.selectedItem = 'Dubious'
        #load in Jmol
        try:
            self.execute(';'.join(['load "%s"' % load_model(self.pdbid, self.source).replace(os.sep, '/')
            ,'select all'
            ,'wireframe only'
            ,'wireframe off'
            ,'select none']))
            self.actionsDict[u'toggle binding site'].selected = prefbool(prefs['bindingsite'])
            self.actionsDict[u'toggle coordinates to exam'].selected = prefbool(prefs['coordstoexam'])
            self.actionsDict[u'toggle ligand'].selected = prefbool(prefs['ligand'])
            self.execute('zoom 0')
        except Exception,  e:
            self.console.sendConsoleMessage("ERROR: " + str(e))
            showErrorDialog(e)

    def displayEDM(self, atoms, name, color):
        if self.source == 'PDB':
            omapfile, sigma_c = EDS_parser.get_EDM(self.pdbid)
            if not omapfile:
                self.console.sendConsoleMessage("EDM unavailable at EDS")
                return
            self.execute('isosurface %s color %s cutoff %s within %s %s "%s" mesh dots fill translucent 0.3' %\
                        (name, color, sigma_c*float(prefs.get('sigma', '1.0')), prefs.get('edmdistance', '2.1'), atoms, omapfile.replace(os.sep, '/')))
        elif self.source == 'PDB_REDO':
            ccp4file = pdb_redo.get_EDM(self.pdbid)
            if ccp4file:
                self.execute('isosurface %s color %s sigma %s within %s %s insideout "%s" mesh dots fill translucent 0.3'%\
                        (name, color, float(prefs.get('sigma', '1.0')), prefs.get('edmdistance', '2.1'), atoms, ccp4file.replace(os.sep, '/')))
            else:
                self.console.sendConsoleMessage("EDM unavailable at %s" % self.source)

        else:
            self.console.sendConsoleMessage("Unable to load EDM from %s" % self.source)

    def displayBindingSite(self, visible=True):
        if not self.binding_site:
            return
        if not visible:
            self.execute(';'.join(['select(binding_site)'
            ,'wireframe only'
            ,'wireframe off'
            ,'select none']))
            if self.binding_site_IS:
                self.execute('isosurface BINDINGSITE off')
                self.binding_site_IS = 0
            return
        binding_site_selection = reslist_to_sel(self.binding_site)
        self.execute(';'.join(['define binding_site (%s)' % binding_site_selection
        ,'select(binding_site)'
        ,'wireframe %s' % prefs.get('bswfv','0.01')
        ,'spacefill %s' % prefs.get('bssfv', 'off')
        ,'color %s' % prefs.get('bscolor', 'white')
        ,'select none']))
        if prefs.get('bindingsite_edm', False):
            if self.binding_site_IS == 0:
                self.execute('isosurface BINDINGSITE on')
            elif not self.binding_site_IS:
                self.displayEDM("{binding_site}", "BINDINGSITE", prefs.get('bsedmcolor', 'cyan'))
            self.binding_site_IS = 1
        elif self.binding_site_IS:
            self.execute('isosurface BINDINGSITE off')
            self.binding_site_IS = 0

    def displayCoordsToExam(self, visible=True):
        if not self.residues_to_exam:
            return
        if not visible:
            self.execute(';'.join(['select(coords_to_exam)'
            ,'wireframe only'
            ,'wireframe off'
            ,'select none']))
            if self.residues_to_exam_IS:
                self.execute('isosurface COORDS_TO_EXAM off')
                self.residues_to_exam_IS = 0
            return
        exam_residues_selection = reslist_to_sel(self.residues_to_exam)
        self.execute(';'.join(['define coords_to_exam (%s)' % exam_residues_selection
        ,'select(coords_to_exam)'
        ,'wireframe %s' % prefs.get('rewfv', '0.1')
        ,'spacefill %s' % prefs.get('resfv', '0.2')
        ,'color %s' % prefs.get('recolor', 'temperature')
        ,'select none']))
        if prefs.get('coordstoexam_edm', True):
            if self.residues_to_exam_IS == 0:
                self.execute('isosurface COORDS_TO_EXAM on')
            elif not self.residues_to_exam_IS:
                self.displayEDM("{coords_to_exam}", "COORDS_TO_EXAM", prefs.get('reedmcolor', 'yellow'))
            self.residues_to_exam_IS = 1
        elif self.residues_to_exam_IS:
            self.execute('isosurface COORDS_TO_EXAM off')
            self.residues_to_exam_IS = 0

    def displayLigand(self, visible=True):
        if not self.ligandresidues:
            return
        if not visible:
            self.execute(';'.join(['select(svligand)'
            ,'wireframe only'
            ,'wireframe off'
            ,'select none']))
            if self.ligandresidues_IS:
                self.execute('isosurface LIGAND off')
                self.ligandresidues_IS = 0
            return
        ligands_selection = reslist_to_sel(self.ligandresidues)
        self.execute(';'.join(['define svligand (%s)' % ligands_selection
        ,'select(svligand)'
        ,'wireframe %s' % prefs.get('ligwfv', '0.1')
        ,'spacefill %s' % prefs.get('ligsfv', '0.2')
        ,'color %s' % prefs.get('ligcolor', 'magenta')
        ,'select none'
        ,'center svligand']))
        if prefs.get('ligand_edm', False):
            if self.ligandresidues_IS == 0:
                self.execute('isosurface LIGAND on')
            elif not self.ligandresidues_IS:
                self.displayEDM("{svligand}", "LIGAND", prefs.get('ligedmcolor', 'red'))
            self.ligandresidues_IS = 1
        elif self.ligandresidues_IS:
            self.execute('isosurface LIGAND off')
            self.ligandresidues_IS = 0

    def updateOutFile(self):
        outfile = open(self.checkedfilename, 'wb')
        csvfile = csv.writer(outfile)
        csvfile.writerow(rsr_analysis.titles)
        self.savedkeys[self.key] = self.resultdict.pop(self.key)
        for d in self.savedkeys, self.resultdict:
            for key in d:
                ligandresidues, residues_to_exam, binding_site, ligandgood, bsgood, source = d[key]
                csvfile.writerow([key.split('|')[0], ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site), ligandgood, bsgood, source])
        outfile.close()

    def start(self):
        self.savedkeys = {}
        self.key = self.structs_cbox.selectedItem
        self.structs_cbox.selectedItem = self.key

    def loadCSV(self):
        restart = False
        values = self.values
        self.resultdict = {}
        csvfilename = values.csvfile
        outdir = None
        if csvfilename:
            if not os.path.isfile(csvfilename):
                print 'File %s does not exist' % csvfilename
                showWarningDialog('File %s does not exist' % csvfilename)
                restart = True
                return restart
            outdir = os.path.dirname(csvfilename)
            basename = os.path.splitext(os.path.basename(csvfilename))[0]
            #Set some global variables
            rsr_analysis.main(values)
        else:
            self.wd.show()
            datawritten = rsr_analysis.main(values)
            self.wd.show(False)
            if values.no_view:
                exit(0)
            csvfilename = values.outputfile
            outdir = os.path.dirname(csvfilename)
            basename = os.path.splitext(os.path.basename(csvfilename))[0]
            rejfile = os.path.join(outdir, basename + '_rejected.txt')
            if os.path.isfile(rejfile):
                showMessageDialog('The reasons can be found in %s' % rejfile, "Some ligands or structures rejected")
            if datawritten:
                showMessageDialog('Analysis data saved to %s' % values.outputfile, 'Analysis completed')
            else:
                showWarningDialog('No structures to be viewed.')
                restart = True
                return restart
        self.checkedfilename = os.path.join(outdir, basename + '_checked.csv')

        print('Loading data from %s...' % csvfilename)
        csvfile = open(csvfilename, 'rb')
        reader = csv.reader(csvfile, 'excel')
        header = reader.next()
        has_source = len(header) > 6
        ngoodbs = 0
        ngoodlig = 0
        ndubbs = 0
        ndublig = 0
        nbadbs = 0
        nbadlig = 0
        for fields in reader:
            id, residues_to_exam_string, ligandresidues_string, binding_site_string, ligandgood, bsgood = fields[:6]
            if has_source:
                source = fields[6]
            else:
                source = 'PDB'
            bsgood = bsgood.lower()
            ligandgood = ligandgood.lower()
            residues_to_exam = residues_to_exam_string.split(';')
            ligandresidues = ligandresidues_string.split(';')
            binding_site = binding_site_string.split(';')
            if bsgood == 'true':
                ngoodbs += 1
            elif bsgood == 'false':
                nbadbs += 1
            elif bsgood == 'dubious':
                ndubbs +=1
            else:
                print "Unrecognized type: %s" % bsgood

            if ligandgood == 'true':
                ngoodlig += 1
            elif ligandgood == 'false':
                nbadlig += 1
            elif ligandgood == 'dubious':
                ndublig +=1
            else:
                print "Unrecognized type: %s" % ligandgood
            self.resultdict[id + '|' +ligandresidues[0]] = [ligandresidues, residues_to_exam, binding_site, ligandgood, bsgood, source]
        else:
            print 'Data loaded'
        csvfile.close()

        #Display results
        print "             \tGood\tDubious\tBad"
        print "Ligands     \t%s\t%s\t%s" % (ngoodlig, ndublig, nbadlig)
        print "Binding sites\t%s\t%s\t%s" % (ngoodbs, ndubbs, nbadbs)

        #Ask about which structures to look at.
        struc_d = StructureSelectDialog(values, (ngoodlig, ndublig, nbadlig, ngoodbs, ndubbs, nbadbs), csvfilename)
        wannasee = struc_d.show()

        check_good_bs = wannasee['bs']['Good']
        check_good_ligand = wannasee['ligand']['Good']
        check_bad_bs = wannasee['bs']['Bad']
        check_bad_ligand = wannasee['ligand']['Bad']
        check_dubious_bs = wannasee['bs']['Dubious']
        check_dubious_ligand = wannasee['ligand']['Dubious']

        if check_good_bs == check_good_ligand == check_bad_bs == check_bad_ligand == check_dubious_bs == check_dubious_ligand == False:
            restart = True
            return restart
        for k, v in self.resultdict.items():
            ligandresidues, residues_to_exam, binding_site, ligandgood, bsgood, source = v
            if not(\
                (check_good_bs and bsgood == 'true') or\
                (check_bad_bs and bsgood == 'false') or\
                (check_dubious_bs and bsgood == 'dubious') or\
                (check_good_ligand and ligandgood == 'true') or\
                (check_bad_ligand and ligandgood == 'false') or\
                (check_dubious_ligand and ligandgood == 'dubious')
            ):
                self.resultdict.pop(k)

        if not self.resultdict:
            print 'File without data! %s' % values.outputfile
            showWarningDialog('No binding sites to be viewed.')
            restart = True
        return restart

    def restart(self):
        if self.frame:
            self.frame.setVisible(False)
            self.frame.dispose()
            self.frame = None
        self.__init__(SettingsDialog().getValues())

class DialogShower(SwingWorker):
    def __init__(self, diag, viewer):
        SwingWorker.__init__(self)
        self.diag=diag
        self.viewer=viewer
    def doInBackground(self):
        self.diag.show(True)
        time.sleep(1)
        while self.viewer.isScriptExecuting():
            time.sleep(1)
    def done(self):
        try:
            self.diag.show(False)
            self.get()  #raise exception if abnormal completion
        except ExecutionException, e:
            raise SystemExit, e.getCause()

class StructureSelectDialog(object):
    def __init__(self, values, stats, csvfilename):
        ngoodlig, ndublig, nbadlig, ngoodbs, ndubbs, nbadbs = stats
        total = sum((ngoodlig, ndublig, nbadlig))
        self.diag = JDialog(JFrame(iconImage=vhelibsicon),size = (500, 200), title = 'Select which structures to view', modal=True)
        self.panel = JPanel(GridBagLayout())
        constraints = GridBagConstraints()
        constraints.weightx = 0.5
        constraints.weighty = 0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 0
        constraints.gridx = 0
        self.panel.add(JLabel('Data loaded from %s' % csvfilename))
        constraints.gridy += 1
        constraints.gridwidth = 2
        self.panel.add(JLabel('Select which structures to view:', ), constraints)
        constraints.gridx += 2
        self.panel.add(JLabel('(%s found)' % total), constraints)
        constraints.gridx -= 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Good Binding Sites (%s/%s)' % (ngoodbs, total), toolTipText='Check Binding Sites which met all the required conditions (%s binding sites)' % ngoodbs), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Good Ligands (%s/%s)' % (ngoodlig, total), toolTipText='Check Ligands which met all the required conditions (%s ligands)' % ngoodlig), constraints)
        constraints.gridx -= 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Bad Binding Sites (%s/%s)' % (nbadbs, total), toolTipText='Check Binding Sites which failed to meet more than %s of the required conditions (%s binding sites)' % (values.tolerance, nbadbs)), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Bad Ligands (%s/%s)' % (nbadlig, total), toolTipText='Check Ligands which failed to meet more than %s of the required conditions (%s ligands)' % (values.tolerance, nbadlig)), constraints)
        constraints.gridx -= 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Dubious Binding Sites (%s/%s)' % (ndubbs, total), toolTipText='Check Binding Sites with up to %s unmet contiditions (%s binding sites)' % (values.tolerance, ndubbs), selected=True), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Dubious Ligands (%s/%s)' % (ndublig, total), toolTipText='Check Ligands with up to %s unmet contiditions (%s ligands)' % (values.tolerance, ndublig), selected=True), constraints)
        constraints.gridy += 1
        constraints.gridx = 0
        constraints.gridwidth = 2
        self.panel.add(JButton('OK', actionPerformed=self.choose), constraints)
        constraints.gridx += 3
        self.panel.add(JButton('Cancel', actionPerformed=self.cancel), constraints)
        self.diag.add(self.panel)
        self.diag.setLocationRelativeTo(None)
        self.diag.pack()
    def show(self, boolean=True):
        self.diag.visible = boolean
        return self.getChoice()

    def getChoice(self):
        checkboxdict = {'bs':{'Good':None, 'Bad':None, 'Dubious':None}, 'ligand':{'Good':None, 'Bad':None, 'Dubious':None}}
        for c in self.panel.components:
            if 'Binding Site' in c.text:
                d = 'bs'
            elif 'Ligand' in c.text:
                d = 'ligand'
            else:
                continue
            quality = c.text.split()[0]
            if quality in checkboxdict[d]:
                checkboxdict[d][quality] = c.selected
        return checkboxdict

    def choose(self, event):
        self.diag.visible = False
    def cancel(self, event):
        for c in self.panel.components:
            if 'Binding' in c.text or 'Ligand' in c.text:
                c.selected = False
        self.choose(event)


class SettingsDialog(object):
    profiles = {'Default (PDB)':{
                            'distance':math.sqrt(rsr_analysis.inner_distance)
                            , 'rsr_lower':rsr_analysis.RSR_lower
                            , 'rsr_upper':rsr_analysis.RSR_upper
                            , 'max_owab':rsr_analysis.OWAB_max
                            , 'min_rscc':rsr_analysis.RSCC_min
                            , 'max_resolution':rsr_analysis.RESOLUTION_max
                            , 'tolerance':rsr_analysis.TOLERANCE
                            , 'min_occupancy':rsr_analysis.OCCUPANCY_min
                            , 'min_rfree':rsr_analysis.RFREE_min
                            , 'use_owab': rsr_analysis.CHECK_OWAB
                            , 'use_res': rsr_analysis.CHECK_RESOLUTION
                            , 'use_pdb_redo': False
                            , 'outputfile':  'vhelibs_analysis_default.csv'
                            , 'editable': False
                            }
                    , 'Default (PDB_REDO)':{
                            'distance':math.sqrt(rsr_analysis.inner_distance)
                            , 'rsr_lower': 0.165
                            , 'rsr_upper':rsr_analysis.RSR_upper
                            , 'max_owab':rsr_analysis.OWAB_max
                            , 'min_rscc':rsr_analysis.RSCC_min
                            , 'max_resolution':rsr_analysis.RESOLUTION_max
                            , 'tolerance':rsr_analysis.TOLERANCE
                            , 'min_occupancy':rsr_analysis.OCCUPANCY_min
                            , 'min_rfree':rsr_analysis.RFREE_min
                            , 'use_owab': False
                            , 'use_res': rsr_analysis.CHECK_RESOLUTION
                            , 'use_pdb_redo': True
                            , 'outputfile':  'vhelibs_analysis_PDB_REDO_default.csv'
                            , 'editable': False
                            }
                    ,'Iridium':{
                            'distance':5
                            , 'rsr_lower':0.1
                            , 'rsr_upper':0.3
                            , 'max_owab':50
                            , 'min_rscc':0.9
                            , 'max_resolution':5
                            , 'tolerance':rsr_analysis.TOLERANCE
                            , 'min_occupancy':1.0
                            , 'min_rfree':0 #FIXME:
                            , 'use_owab': False
                            , 'use_res': False
                            , 'outputfile':  'vhelibs_analysis_iridium.csv'
                            , 'use_pdb_redo': False
                            , 'editable': False
                            }
                    ,'Custom':{
                        'outputfile':  'vhelibs_analysis_custom.csv'
                        , 'editable': True
                        }
                    }
    def __init__(self, args=['--no-args']):
        self.values = argparser.parse_args(args)
        self.panel = JPanel(GridBagLayout())
        for k in self.profiles['Default (PDB)']:
            if k in ('use_owab', 'use_res', 'use_pdb_redo', 'editable'):
                continue
            else:
                self.__dict__[k] = JTextField(str(self.values.__dict__[k]))

        constraints = GridBagConstraints()
        constraints.weightx = 0.5
        constraints.weighty = 0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 0
        constraints.gridx = 0

        constraints.gridwidth = 2
        infolabel = JLabel('Set options and select the structures to be checked')
        self.panel.add(infolabel, constraints)
        constraints.gridwidth = 1

        constraints.gridy += 1
        tooltip = 'Each profile provides a different set of values to use. The "Custom" profile lets you edit them'
        self.panel.add(JLabel(u'Select a predefined profile', toolTipText=tooltip), constraints)
        self.distance.toolTipText=tooltip
        constraints.gridx += 1
        self.profilecbb = JComboBox(self.profiles.keys(), toolTipText=tooltip)
        self.profilecbb.selectedItem = 'Custom'
        self.profilecbb.addActionListener(make_listen(self.load_profile))
        self.panel.add(self.profilecbb, constraints)
        constraints.gridx -= 1
        ##############
        constraints.gridy += 1
        tooltip="Save the current profile to a file"
        self.panel.add(JButton('Export profile', toolTipText=tooltip, actionPerformed=self.export_current_profile), constraints)
        constraints.gridx += 1
        tooltip="Load a profile from a file"
        self.panel.add(JButton('Import profile', toolTipText=tooltip, actionPerformed=self.import_profile), constraints)
        constraints.gridx -= 1
        ##############

        ##############
        constraints.gridy += 1
        constraints.gridwidth = 2
        tooltip="Models will be loaded from PDB_REDO instead of the PDB"
        self.use_pdb_redo = JCheckBox("Use models from PDB_REDO", toolTipText=tooltip, selected=self.profiles['Default (PDB)']['use_pdb_redo'], actionPerformed=self._check_pdbredo_owab)
        self.use_pdb_redo.toolTipText=tooltip
        self.panel.add(self.use_pdb_redo, constraints)
        constraints.gridwidth = 1
        ##############

        constraints.gridy += 1
        distancetooltip = 'Residues with at least one atom within this distance from any atom of the ligand will be considered as part of the binding site'
        self.panel.add(JLabel(u'Radius (in Å)', toolTipText=distancetooltip), constraints)
        self.distance.toolTipText=distancetooltip
        constraints.gridx += 1
        self.panel.add(self.distance, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        hrsrtooltip="Ligands and binding sites with at least one 'residue' with an RSR above this value will have their score increased by 1"
        self.panel.add(JLabel('Upper cap for RSR', toolTipText=hrsrtooltip), constraints)
        constraints.gridx += 1
        self.rsr_upper.toolTipText=hrsrtooltip
        self.panel.add(self.rsr_upper, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        lrsrtooltip="Ligands and binding sites with a 'residue' with an RSR above this value will have their score increased by 1"
        self.panel.add(JLabel('Good RSR cap', toolTipText=lrsrtooltip), constraints)
        constraints.gridx += 1
        self.rsr_lower.toolTipText=lrsrtooltip
        self.panel.add(self.rsr_lower, constraints)
        constraints.gridx -= 1

        ##############
        constraints.gridy += 1
        tooltip="Ligands and residues with an RSCC below this value will have their score increased by 1"
        self.panel.add(JLabel('Lowest good RSCC', toolTipText=tooltip), constraints)
        constraints.gridx += 1
        self.min_rscc.toolTipText=tooltip
        self.panel.add(self.min_rscc, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        tooltip="Residues and ligands with an average occupancy below this value will have their score increased by 1"
        self.panel.add(JLabel('Average Occupancy', toolTipText=tooltip), constraints)
        constraints.gridx += 1
        self.min_occupancy.toolTipText=tooltip
        self.panel.add(self.min_occupancy, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        tooltip="Residues and ligands with an occupancy-weighted B-factor (OWAB) above this value will have their score increased by 1"
        self.owab_cb = JCheckBox('OWAB', toolTipText=tooltip, selected=(self.profiles['Default (PDB)']['max_owab'] != self.values.max_owab) or self.profiles['Default (PDB)']['use_owab'], actionPerformed=self._check_pdbredo_owab)
        self.panel.add(self.owab_cb, constraints)
        constraints.gridx += 1
        self.max_owab.toolTipText=tooltip
        self.panel.add(self.max_owab, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        tooltip="All residues and ligands from structures with a resolution above this value will have their score increased by 1"
        self.res_cb = JCheckBox("Resolution limit", toolTipText=tooltip, selected=(self.profiles['Default (PDB)']['max_resolution'] != self.values.max_resolution) or self.profiles['Default (PDB)']['use_res'])
        self.panel.add(self.res_cb, constraints)
        constraints.gridx += 1
        self.max_resolution.toolTipText=tooltip
        self.panel.add(self.max_resolution, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        tooltip="All residues and ligands from structures with an R-free value below this value will have their score increased by 1"
        self.panel.add(JLabel('Minimum R-free value', toolTipText=tooltip), constraints)
        constraints.gridx += 1
        self.min_rfree.toolTipText=tooltip
        self.panel.add(self.min_rfree, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        tooltip="<html><body>Each ligand and binding site will be classified as follows: <br>* Good if the score is 0.<br>* Dubious if the score is lower than the tolerance <br>* Bad if the score is higher than the tolerance</html>"
        self.panel.add(JLabel('Tolerance', toolTipText=tooltip), constraints)
        constraints.gridx += 1
        self.tolerance.toolTipText=tooltip
        self.panel.add(self.tolerance, constraints)
        constraints.gridx -= 1
        ##############

        constraints.gridy += 1
        outtt="Select or enter the name of the output file and where to save it."
        self.panel.add(JButton('Output file name', toolTipText=outtt, actionPerformed=self.selectOutFileName), constraints)
        constraints.gridx += 1
        self.outputfile.toolTipText=outtt
        self.outputfile.text = os.path.abspath(self.outputfile.text)
        self.panel.add(self.outputfile, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        constraints.gridwidth = 2
        tooltip="Data will be saved to and loaded from %s" % os.path.join(os.path.expanduser('~'), '.vhelibs_cache')
        self.use_cache = JCheckBox("Do not download already downloaded files", toolTipText=tooltip, selected=(self.profiles['Default (PDB)']['max_resolution'] != self.values.max_resolution) or self.profiles['Default (PDB)']['use_res'])
        self.use_cache.toolTipText=tooltip
        self.panel.add(self.use_cache, constraints)
        constraints.gridwidth = 1

        constraints.gridy += 1
        pdbtt="Parse structures from their PDB codes, either by entering their codes or by providing a file containing them"
        self.panel.add(JButton('Use PDB codes', toolTipText=pdbtt, actionPerformed=self.loadStructsFrom), constraints)
        constraints.gridx += 1
        uniprottt="Parse structures from their UniProtKB codes, either by entering their codes or by providing a file containing them"
        self.panel.add(JButton('Use UniProtKB codes', toolTipText=uniprottt, actionPerformed=self.loadStructsFrom), constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        slett="Writes to a file a list of the ligands that would be excluded from the analysis, i.e. solvents, ions"
        self.panel.add(JButton('Save ligands to exclude', toolTipText=slett, actionPerformed=self.writeLe), constraints)
        constraints.gridx += 1
        llett="Loads from a file a list of ligands to be excluded from the analysis"
        self.panel.add(JButton('Load ligands to exclude', toolTipText=llett, actionPerformed=self.readLe), constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        constraints.gridwidth = 2
        csvfilett = "Load a previously generated file to review its models"
        self.panel.add(JButton('Load previous results file',  toolTipText=csvfilett, actionPerformed=self.csvFileDialog), constraints)
        self.frame = JFrame(iconImage=vhelibsicon,size = (0, 0), title = TITLE)
        self.diag = JDialog(self.frame,size = (500, 200), title = TITLE, modal=True)
        self.diag.add(self.panel)
        self.frame.setLocationRelativeTo(None)
        self.diag.setLocationRelativeTo(self.frame)
        self.frame.setUndecorated(True)
        self.diag.pack()
        self.profilecbb.selectedItem = 'Default (PDB)'
        #self.load_profile(e=None, profilename='Default (PDB)')

    def _check_pdbredo_owab(self, event):
        """
        No OWAB data for PDB_REDO
        """
        if self.owab_cb.selected and self.use_pdb_redo.selected:
            if event.source is self.owab_cb:
                self.use_pdb_redo.selected = False
            else:
                self.owab_cb.selected = False

    def export_current_profile(self, event=None):
        filename = str(showOpenDialog(SimpleFileFilter('.tsv', None, 'VHELIBS profile files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
        if not filename or filename == 'None':
            return
        if not filename.endswith('.tsv'):
            filename += '.tsv'
        profile = self.profiles[self.profilecbb.selectedItem]
        for key in self.profiles['Default'].keys():
            if key not in profile:
                if key == 'use_owab':
                    profile[key] = self.owab_cb.selected
                elif key == 'use_res':
                    profile[key] = self.res_cb.selected
                elif key == 'use_pdb_redo':
                    profile[key] = self.use_pdb_redo.selected
                else:
                    profile[key] = self.__dict__[key].text
        file = open(filename, 'wb')
        csvwriter = csv.writer(file, delimiter="\t")
        for k, v in profile.items():
            csvwriter.writerow([k, v])
        file.close()

    def import_profile(self, event = None):
        filename = str(showOpenDialog(SimpleFileFilter('.tsv', None, 'VHELIBS profile files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
        if filename and filename != 'None':
            file = open(filename, 'rb')
            profilename = os.path.splitext(os.path.basename(filename))[0]
            self.profiles[profilename] = {}
            for k, v in csv.reader(file, delimiter="\t"):
                if v == 'False':
                    v = False
                elif v == 'True':
                    v = True
                self.profiles[profilename][k] = v
            file.close()
            self.profilecbb.addItem(profilename)
            self.profilecbb.selectedItem = profilename

    def load_profile(self,e=None, profilename='Default'):
        if e and e.actionCommand != u'comboBoxChanged': return
        if e:
            profilename = self.profilecbb.selectedItem
        profile = self.profiles[profilename]
        for k,  v in profile.items():
            if k == 'editable':
                editable = v
                for w in self.profiles['Default (PDB)'].keys():
                    if w == 'editable': continue
                    elif w == 'use_owab':
                        self.owab_cb.enabled = editable
                    elif w == 'use_res':
                        self.res_cb.enabled = editable
                    elif w == 'outputfile':
                        continue
                    else:
                        self.__dict__[w].enabled = editable
            elif k == 'use_owab':
                self.owab_cb.selected = v
            elif k == 'use_res':
                self.res_cb.selected = v
            elif k == 'use_pdb_redo':
                self.use_pdb_redo.selected = v
            elif k == 'outputfile':
                self.outputfile.text = os.path.abspath(v)
            else:
                self.__dict__[k].text = str(v)

    def writeLe(self, e):
        lefn = str(showOpenDialog(SimpleFileFilter('.csv', None, 'CSV files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
        if lefn and lefn != 'None':
            try:
                print 'Saved to', lefn
                rsr_analysis.cofactors.dump_lists(lefn)
            except Exception, e:
                showErrorDialog(e)

    def readLe(self, e):
        lefn = str(showOpenDialog(SimpleFileFilter('.csv', None, 'CSV files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
        if lefn and lefn != 'None' and os.path.isfile(lefn):
            try:
                print "Loading from", lefn
                rsr_analysis.cofactors.load_lists(lefn)
            except Exception, e:
                showErrorDialog(e)

    def selectOutFileName(self, event):
        outfn = str(showOpenDialog(SimpleFileFilter('.csv', None, 'CSV files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
        if outfn and outfn != 'None':
            if not outfn.endswith('.csv'):
                outfn += '.csv'
            self.outputfile.text = outfn

    def loadStructsFrom(self, event):
        idstring = None
        if 'PDB' in event.source.text:
            source = 'PDB'
            def updateval(v): self.values.pdbids = v
        else:
            source = 'UniProtKB'
            def updateval(v): self.values.swissprot = v
        label = 'Enter %s codes: (must be separated by commas, white space or tabs)' % source
        options = ('Load %s codes from file' % source, 'Enter %s codes manually' %  source, 'Cancel')
        choice = JOptionPane.showOptionDialog(None, 'Select %s codes source' % source,'', JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE, None, options, options[1])
        if choice == 1:
            idstring = JOptionPane.showInputDialog(label)
        elif choice == 0:
            file = str(showOpenDialog(None,multiselect=False))
            if not file or not os.path.isfile(file):
                showErrorDialog('%s is not a readable file' % file)
                return
            elif file == 'None':
                return
            idstring = ' '.join([line.strip() for line in open(file, 'rb')])
        else:
            return
        if idstring:
            ids = idstring.replace(',', ' ').split()
            updateval(ids)
            self.go()

    def csvFileDialog(self, event):
        csvfilter = SimpleFileFilter('.csv', None, 'CSV files')
        file = showOpenDialog(csvfilter, prefs=prefs, prefkey='loadedFiles', multiselect=False)
        if file:
            self.values.csvfile = str(file)
        self.go()

    def load_settings(self):
        self.values.distance = float(self.distance.text)
        self.values.rsr_lower = float(self.rsr_lower.text)
        self.values.rsr_upper = float(self.rsr_upper.text)
        self.values.outputfile = self.outputfile.text
        self.values.min_rscc = float(self.min_rscc.text)
        self.values.min_rfree = float(self.min_rfree.text)
        self.values.min_occupancy = float(self.min_occupancy.text)
        self.values.use_cache = self.use_cache.selected
        self.values.use_pdb_redo = self.use_pdb_redo.selected
        if self.owab_cb.selected:
            self.values.max_owab = float(self.max_owab.text)
        else:
            self.values.max_owab = None
        if self.res_cb.selected:
            self.values.max_resolution = float(self.max_resolution.text)
        else:
            self.values.max_resolution = None
        self.values.tolerance = int(self.tolerance.text)

    def isViable(self):
        return bool(self.values.csvfile or self.values.pdbidfile or self.values.pdbids or self.values.swissprot)

    def go(self):
        self.load_settings()
        if self.isViable():
            self.frame.visible = self.diag.visible = 0

    def getValues(self):
        if not self.isViable():
            self.frame.visible = self.diag.visible = 1
            if not self.isViable():
                exit(1)
        return self.values

class WaitDialog(Runnable):
    def __init__(self, parent=None, info=None, modal=False):
        self.info = info if info else '<html>Calculating binding sites and retrieving relevant data<br /> Please be patient</html>'
        i = UIManager.getIcon("OptionPane.informationIcon")
        self.icon = JLabel(i) if i is not None else JLabel(_infoicon)
        self.frame =  JFrame(iconImage=vhelibsicon)
        self.panel = JPanel()
        self.panel.add(self.icon)
        self.panel.add(JLabel(self.info))
        self.dialog = JDialog(self.frame,'Please wait', modal)
        self.dialog.add(self.panel)
        self.dialog.setLocationRelativeTo(parent)
        self.dialog.pack()
    def show(self, boolean=True):
         self.dialog.visible = boolean
    def run(self):
        self.show(True)

class AboutDialog(object):
    def __init__(self, parent=None):
        self.parentframe = parent.frame if parent else None
        urvfn = 'URV.png'
        urv = java.lang.ClassLoader.getSystemClassLoader().getResource(urvfn)
        if urv:
            u = URL(str(urv))
            self.urvicon = ImageIcon(u)
        else:
           self.urvicon=ImageIcon(urvfn)
        ctnsfn = 'logo_ctns.png'
        ctns = java.lang.ClassLoader.getSystemClassLoader().getResource(ctnsfn)
        if ctns:
            u = URL(str(ctns))
            self.ctnsicon = ImageIcon(u)
        else:
            self.ctnsicon=ImageIcon(ctnsfn)
        self.frame =  JFrame(iconImage=vhelibsicon)
        self.panel = JPanel(GridBagLayout())
        constraints = GridBagConstraints()
        constraints.weightx = 0.5
        constraints.weighty = 0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)

        constraints.gridy = 0
        constraints.gridx = 0
        constraints.gridwidth = 1
        self.panel.add(JLabel(ImageIcon(vhelibsicon)), constraints)
        constraints.gridx += 1
        constraints.gridwidth = 2

        label = dedent(u"""<html><body>
        <div style="text-align: center;"><big><big><big><span style="font-weight: bold;">VHELIBS %s</span></big></big></big><br>
        <br>
        <big><span style="font-weight: bold;">V</span>alidation <span style="font-weight: bold;">HE</span>lper for <span style="font-weight: bold;">LI</span>gands and <span style="font-weight: bold;">B</span>inding <span style="font-weight: bold;">S</span>ites</big>
        <br>
        <p style="margin-top: 0.12cm; margin-bottom: 0cm; line-height: 0.53cm;" lang="en-US"><small><small>
        <font face="Helvetica-Light, Arial Unicode MS, sans-serif"><font style="font-size: 13pt;" size="3"><small><small>Adrià
        Cereto-Massagué<sup> 1</sup>, María José Ojeda<sup>1</sup>, Robbie P. Joosten<sup>2</sup>,
        Cristina Valls<sup> 1</sup>, Miquel Mulero<sup> 1</sup>, M. Josepa
        Salvado<sup>1</sup>, Anna Arola-Arnal<sup>1</sup>, Lluís Arola<sup>1,
        3</sup>, Santiago Garcia-Vallvé<sup>1, 3</sup> and Gerard Pujadas<sup>
        1, 3,</sup></small></small></font></font></small></small></p>
        <p style="margin-top: 0.07cm; margin-bottom: 0.09cm;" ><small><small><font face="Helvetica-Light, Arial Unicode MS, sans-serif"><sup>1</sup>Grup
        de Recerca en Nutrigenòmica, Departament de Bioquímica i
        Biotecnologia, Universitat Rovira i Virgili, Campus de Sescelades, C/
        Marce&#320;lí Domingo s/n, 43007 Tarragona, Catalonia, Spain</font></small></small></p>
        <p style="margin-top: 0.07cm; margin-bottom: 0.09cm;" ><small><small><font face="Helvetica-Light, Arial Unicode MS, sans-serif"><sup>2</sup>Department of Biochemistry, Netherlands Cancer Institute, Plesmanlaan 121, 1066 CX Amsterdam, The Netherlands</font></small></small></p>
        <p style="margin-top: 0.07cm; margin-bottom: 0.09cm;" ><small><small><font face="Helvetica-Light, Arial Unicode MS, sans-serif"><sup>3</sup>Centre
        Tecnològic de Nutrició i Salut (CTNS), TECNIO, CEICS, Avinguda
        Universitat 1, 43204, Reus, Catalonia, Spain</font></small></small></p>
        <br>
        <br>
        </div>
        </body></html>""") % VHELIBS_VERSION


        self.panel.add(JLabel(label), constraints)
        constraints.gridy += 1

        constraints.gridx = 0
        constraints.gridwidth = 1
        constraints.gridheight = 2
        self.panel.add(JLabel(self.urvicon), constraints)

        constraints.gridx += 2
        self.panel.add(JLabel(self.ctnsicon), constraints)

        constraints.gridx -= 1
        constraints.gridheight = 1
        self.panel.add(JLabel(u"""<html>More information, help and documentation can be found at:</html>"""), constraints)
        constraints.gridy += 1

        website = JLabel()
        website.text = "<html><a href=\"\">urvnutrigenomica-ctns.github.com/VHELIBS</a></html>"
        website.cursor = Cursor(Cursor.HAND_CURSOR)
        website.mouseClicked = lambda e: Desktop.getDesktop().browse(URI("http://urvnutrigenomica-ctns.github.com/VHELIBS"))
        self.panel.add(website, constraints)

        self.diag = JDialog(self.frame, title = 'About VHELIBS')
        self.diag.setLocationRelativeTo(self.parentframe)
        self.diag.add(self.panel)
        self.diag.preferredSize = (660,420)
        self.diag.pack()
    def show(self):
        self.diag.visible = True

class DisplaySettingsDialog(object):
    keys = ('ligwfv', 'ligsfv', 'ligcolor', 'ligedmcolor', 'bswfv', 'bssfv', 'bscolor', 'bsedmcolor', 'rewfv', 'resfv', 'recolor', 'reedmcolor', 'edmdistance', 'sigma')
    def __init__(self, parent):
        self.parent = parent
        minimum_size = (100, 20)
        self.frame =  JFrame(iconImage=vhelibsicon)
        self.panel = JPanel(GridBagLayout())
        constraints = GridBagConstraints()
        constraints.weightx = 0.5
        constraints.weighty = 0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 0
        constraints.gridx = 1
        self.panel.add(JLabel('Ligand'), constraints)

        constraints.gridy = 1
        self.ligwfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.ligwfv, constraints)

        constraints.gridy = 2
        self.ligsfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.ligsfv, constraints)

        constraints.gridy = 3
        self.ligcolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.ligcolor, constraints)

        constraints.gridy = 4
        self.ligedmcolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.ligedmcolor, constraints)

        constraints.gridy = 0
        constraints.gridx = 2
        self.panel.add(JLabel('Binding Site'), constraints)

        constraints.gridy = 1
        self.bswfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.bswfv, constraints)

        constraints.gridy = 2
        self.bssfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.bssfv, constraints)

        constraints.gridy = 3
        self.bscolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.bscolor, constraints)

        constraints.gridy = 4
        self.bsedmcolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.bsedmcolor, constraints)

        constraints.gridy = 0
        constraints.gridx = 3
        self.panel.add(JLabel('Coordinates to exam'), constraints)

        constraints.gridy = 1
        self.rewfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.rewfv, constraints)

        constraints.gridy = 2
        self.resfv = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.resfv, constraints)

        constraints.gridy = 3
        self.recolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.recolor, constraints)

        constraints.gridy = 4
        self.reedmcolor = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.reedmcolor, constraints)

        constraints.gridx = 0
        constraints.gridy = 1
        self.panel.add(JLabel('Wireframe'), constraints)
        constraints.gridy = 2
        self.panel.add(JLabel('Spacefill'), constraints)
        constraints.gridy = 3
        self.panel.add(JLabel('Color'), constraints)
        constraints.gridy = 4
        self.panel.add(JLabel('EDM Color'), constraints)
        constraints.gridy = 5
        constraints.insets = Insets(15,3,3,3)
        self.panel.add(JLabel(u'EDM Radius (in Å)', toolTipText="Distance within which the Electron Density Map will be showed"), constraints)
        constraints.gridx = 1
        self.edmdistance = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.edmdistance, constraints)
        constraints.gridx = 2
        self.panel.add(JLabel('EDM sigma', toolTipText="Contour level of the Electron Density Map"), constraints)
        constraints.gridx = 3
        self.sigma = JTextField(minimumSize=minimum_size, preferredSize=minimum_size)
        self.panel.add(self.sigma, constraints)
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 6
        constraints.gridx = 1
        constraints.gridwidth = 2
        self.panel.add(JButton('Apply', actionPerformed=self.saveprefs), constraints)
        constraints.gridwidth = 1
        constraints.gridx = 3
        self.panel.add(JButton('Defaults', actionPerformed=self.loaddefaults), constraints)
        self.diag = JDialog(self.frame, title = 'Display Settings')
        self.diag.setLocationRelativeTo(self.parent.frame)
        self.diag.add(self.panel)
        self.diag.pack()

    def loadprefs(self):
        self.ligwfv.text = prefs.get('ligwfv', '0.2')
        self.ligsfv.text = prefs.get('ligsfv', '0.4')
        self.ligcolor.text = prefs.get('ligcolor', 'magenta')
        self.ligedmcolor.text = prefs.get('ligedmcolor', 'red')

        self.bswfv.text = prefs.get('bswfv', '0.01')
        self.bssfv.text = prefs.get('bssfv', 'off')
        self.bscolor.text = prefs.get('bscolor', 'white')
        self.bsedmcolor.text = prefs.get('bsedmcolor', 'cyan')

        self.rewfv.text = prefs.get('rewfv', '0.2')
        self.resfv.text = prefs.get('resfv', '0.4')
        self.recolor.text = prefs.get('recolor', 'temperature')
        self.reedmcolor.text = prefs.get('reedmcolor', 'yellow')

        self.edmdistance.text = prefs.get('edmdistance', '2.1')
        self.sigma.text = prefs.get('sigma', '1.0')

        self.saveprefs()

    def saveprefs(self, event=None):
        displaybs = prefbool(prefs['bindingsite'])
        displaycoords = prefbool(prefs['coordstoexam'])
        displaylig = prefbool(prefs['ligand'])

        if (prefs['edmdistance'], prefs['sigma']) != (self.edmdistance.text, self.sigma.text):
            prefs['edmdistance'] = self.edmdistance.text
            prefs['sigma'] = self.sigma.text
            self.parent.ligandresidues_IS = self.parent.residues_to_exam_IS = self.parent.binding_site_IS = None
            if displaybs:
                self.parent.actionsDict[u'toggle binding site'].selected = False
            if displaycoords:
                self.parent.actionsDict[u'toggle coordinates to exam'].selected = False
            if displaylig:
                self.parent.actionsDict[u'toggle ligand'].selected = False

        if (prefs['bswfv'], prefs['bssfv'], prefs['bscolor']) != (self.bswfv.text, self.bssfv.text, self.bscolor.text):
            prefs['bswfv'] = self.bswfv.text
            prefs['bssfv'] = self.bssfv.text
            prefs['bscolor'] = self.bscolor.text
            self.parent.actionsDict[u'toggle binding site'].selected = False
        if prefs['bsedmcolor'] != self.bsedmcolor.text:
            prefs['bsedmcolor'] = self.bsedmcolor.text
            self.parent.binding_site_IS = None
            if displaybs:
                self.parent.actionsDict[u'toggle binding site'].selected = False

        if (prefs['rewfv'], prefs['resfv'], prefs['recolor']) != (self.rewfv.text, self.resfv.text, self.recolor.text):
            prefs['rewfv'] = self.rewfv.text
            prefs['resfv'] = self.resfv.text
            prefs['recolor'] = self.recolor.text
            self.parent.actionsDict[u'toggle coordinates to exam'].selected = False
        if prefs['reedmcolor'] != self.reedmcolor.text:
            prefs['reedmcolor'] = self.reedmcolor.text
            self.parent.residues_to_exam_IS = None
            if displaycoords:
                self.parent.actionsDict[u'toggle coordinates to exam'].selected = False

        if (prefs['ligwfv'], prefs['ligsfv'], prefs['ligcolor']) != (self.ligwfv.text, self.ligsfv.text, self.ligcolor.text):
            prefs['ligwfv'] = self.ligwfv.text
            prefs['ligsfv'] = self.ligsfv.text
            prefs['ligcolor'] = self.ligcolor.text
            self.parent.actionsDict[u'toggle ligand'].selected = False
        if prefs['ligedmcolor'] != self.ligedmcolor.text:
            prefs['ligedmcolor'] = self.ligedmcolor.text
            self.parent.ligandresidues_IS = None
            if displaylig:
                self.parent.actionsDict[u'toggle ligand'].selected = False

        self.parent.actionsDict[u'toggle binding site'].selected = displaybs
        self.parent.actionsDict[u'toggle coordinates to exam'].selected = displaycoords
        self.parent.actionsDict[u'toggle ligand'].selected = displaylig


    def show(self, doit=True):
        self.loadprefs()
        self.diag.visible = doit

    def loaddefaults(self, event=None):
        for key in self.keys:
            if key in prefs:
                prefs.remove(key)
        self.loadprefs()

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
                sellist.append("[%s]%s:%s" % (resname, resnum, chain))
            except IndexError:
                print "Malformed residue string:"
                print res
    if sellist:
        return ' or '.join(sellist)
    else:
        return 'none'

def main(args=sys.argv):
    """
    """
    ttm = ToolTipManager.sharedInstance()
    ttm.reshowDelay = 0
    ttm.dismissDelay = 100000
    values = SettingsDialog(args).getValues()
    struva = StruVa(values)
    return struva

if __name__ == '__main__':
    try:
        sv = main()
    except Exception, e:
        print e
        showErrorDialog(e)
