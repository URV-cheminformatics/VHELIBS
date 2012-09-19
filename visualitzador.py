# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#

#Python stuff
import sys
class OutLog(object):
    def __init__(self, filename):
        self.file = open(filename, 'w')
    def write(self, text):
        sys.__stdout__.write(text)
        self.file.write(text)
        self.file.flush()
    def close(self):
        self.file.close()
sys.stdout = OutLog('VHELIBS_log.txt')
######One-jar magic#######
sys.path.append('__pyclasspath__/pylib')
if not sys.prefix:
    sys.prefix='.'
######One-jar magic#######
import os
import csv
import math
import time
from textwrap import dedent
from sys import exit
#Java stuff
import java
from java.net import URL
from java.lang import Runnable
from java.awt import BorderLayout, Dimension, GridLayout, GridBagLayout, GridBagConstraints, Insets
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

VHELIBS_VERSION = "0.99b"

#Own stuff
sys.argv = [arg for arg in sys.argv if __file__ not in arg]
if not len(sys.argv):
    sys.argv.append('--no-args')
import rsr_analysis
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
        self.actionsDict = {}
        self.wd = WaitDialog()
        if self.values:
            try:
                self.loadCSV()
                if not self.viewer:
                    SwingUtilities.invokeAndWait(self)
                    self.start()
            except Exception, e:
                self.e = e
                print e
                showErrorDialog('Unable to load RSR analysis results file:\n %s' % str(e))

    def run(self):
        if self.values:
            self.setupUi()
            self.setVisible(True)

    def setupUi(self):
        self.frame = JFrame("VHELIBS", iconImage=vhelibsicon, defaultCloseOperation = JFrame.EXIT_ON_CLOSE, size = (700, 410))
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
        ,'set antialiasTranslucent ON']))
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

        self.structs_cbox = JComboBox(self.resultdict.keys())
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
        self.setVisible = self.frame.setVisible
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
                self.resultdict[self.key][-1] = True
            elif bs_valid == 'bad':
                self.resultdict[self.key][-1] = False
            else:
                self.resultdict[self.key][-1] = bs_valid

            if ligand_valid == 'good':
                self.resultdict[self.key][-2] = True
            elif ligand_valid == 'bad':
                self.resultdict[self.key][-2] = False
            else:
                self.resultdict[self.key][-2] = ligand_valid
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

    def displayAbout(self, event):
        self.aboutdiag.show()
        pass

    def clean(self):
        #Clean everything
        self.execute('delete;')
        game_over = JOptionPane.showConfirmDialog(self.frame,u'Continue working with other structures?',u'No more structures to view!',JOptionPane.YES_NO_OPTION)
        if game_over == JOptionPane.OK_OPTION:
            self.restart()
        else:
            exit(0)

    def reloadStruct(self, event=None):
        if event and event.actionCommand != u'comboBoxChanged': return
        self.key = self.structs_cbox.selectedItem
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
        self.ligandresidues, self.residues_to_exam, self.binding_site,  self.ligandgood, self.bsgood = self.resultdict[self.key]
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
            self.execute(';'.join(['load "=%s"' % self.pdbid
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
                self.execute('isosurface BINDINGSITE color %s sigma %s within %s {binding_site} "=%s" mesh dots fill translucent 0.3' %\
                            (prefs.get('bsedmcolor', 'cyan'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.1'),  self.pdbid))
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
        ,'color %s' % prefs.get('recolor', 'cpk')
        ,'select none']))
        if prefs.get('coordstoexam_edm', True):
            if self.residues_to_exam_IS == 0:
                self.execute('isosurface COORDS_TO_EXAM on')
            elif not self.residues_to_exam_IS:
                self.execute('isosurface COORDS_TO_EXAM color %s sigma %s within %s {coords_to_exam} "=%s" mesh dots fill translucent 0.3' %\
                        (prefs.get('reedmcolor', 'yellow'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.1'),  self.pdbid))
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
                self.execute('isosurface LIGAND color %s sigma %s within %s {svligand} "=%s" mesh dots fill translucent 0.3' %\
                        (prefs.get('ligedmcolor', 'red'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.1'),  self.pdbid))
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
                ligandresidues, residues_to_exam, binding_site, ligandgood, bsgood = d[key]
                csvfile.writerow([key.split('|')[0], ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site), ligandgood, bsgood])
        outfile.close()

    def start(self):
        self.savedkeys = {}
        self.key = self.structs_cbox.selectedItem
        self.structs_cbox.selectedItem = self.key

    def loadCSV(self):
        values = self.values
        self.resultdict = {}
        csvfilename = values.csvfile
        if csvfilename:
            if not os.path.isfile(csvfilename):
                print 'File %s does not exist' % csvfilename
                showWarningDialog('File %s does not exist' % csvfilename)
                self.restart()
            outdir = os.path.dirname(csvfilename)
            basename = os.path.splitext(os.path.basename(csvfilename))[0]
        else:
            self.wd.show()
            datawritten = rsr_analysis.main(
                                            values.pdbidfile
                                            , pdbidslist = values.pdbids
                                            , swissprotlist =values.swissprot
                                            , rsr_upper=values.rsr_upper
                                            , rsr_lower = values.rsr_lower
                                            , distance=values.distance
                                            , outputfile = values.outputfile
                                            )
            self.wd.show(False)
            if values.no_view:
                exit(0)
            if datawritten:
                showMessageDialog('Analysis data saved to %s' % datawritten, 'Analysis completed')
                csvfilename = values.outputfile
                outdir = os.path.dirname(csvfilename)
                basename = os.path.splitext(os.path.basename(csvfilename))[0]
            else:
                showWarningDialog('No structures to be viewed.')
                self.restart()
                return
        #Ask about which structures to look at.
        struc_d = StructureSelectDialog(values)
        wannasee = struc_d.show()

        check_good_bs = wannasee['bs']['Good']
        check_good_ligand = wannasee['ligand']['Good']
        check_bad_bs = wannasee['bs']['Bad']
        check_bad_ligand = wannasee['ligand']['Bad']
        check_dubious_bs = wannasee['bs']['Dubious']
        check_dubious_ligand = wannasee['ligand']['Dubious']

        if check_good_bs == check_good_ligand == check_bad_bs == check_bad_ligand == check_dubious_bs == check_dubious_ligand == False:
            self.restart()
            return

        print('Loading data from %s...' % csvfilename)
        csvfile = open(csvfilename, 'rb')
        reader = csv.reader(csvfile, 'excel')
        for fields in reader:
            id, residues_to_exam_string, ligandresidues_string, binding_site_string, ligandgood, bsgood = fields
            if id != 'PDB ID':
                bsgood = bsgood.lower()
                ligandgood = ligandgood.lower()
                cont = True
                if not(\
                    (check_good_bs and bsgood == 'true') or\
                    (check_bad_bs and bsgood == 'false') or\
                    (check_dubious_bs and bsgood == 'dubious') or\
                    (check_good_ligand and ligandgood == 'true') or\
                    (check_bad_ligand and ligandgood == 'false') or\
                    (check_dubious_ligand and ligandgood == 'dubious')
                ):
                    continue
                residues_to_exam = residues_to_exam_string.split(';')
                ligandresidues = ligandresidues_string.split(';')
                binding_site = binding_site_string.split(';')
                self.resultdict[id + '|' +ligandresidues[0]] = [ligandresidues, residues_to_exam, binding_site, ligandgood, bsgood]
        else:
            print 'Data loaded'
        csvfile.close()
        self.checkedfilename = os.path.join(outdir, basename + '_checked.csv')
        if not self.resultdict:
            print 'File without data! %s' % self.checkedfilename
            showWarningDialog('No structures to be viewed.')
            self.restart()

    def restart(self):
        self.__init__(SettingsDialog().getvalues())

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
    def __init__(self, values):
        self.diag = JDialog(JFrame(iconImage=vhelibsicon),size = (500, 200), title = 'Select which structures to view', modal=True)
        self.panel = JPanel(GridBagLayout())
        constraints = GridBagConstraints()
        constraints.weightx = 0.5
        constraints.weighty = 0.5
        constraints.fill = GridBagConstraints.HORIZONTAL
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 0
        constraints.gridx = 0
        constraints.gridwidth = 5
        self.panel.add(JLabel('Select which structures to view'), constraints)
        constraints.gridwidth = 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Good Binding Site', toolTipText='Check Binding Sites whith all residues with an RSR < %s' % values.rsr_lower), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Good Ligand', toolTipText='Check Ligands with RSR < %s' % values.rsr_lower), constraints)
        constraints.gridx -= 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Bad Binding Site', toolTipText='Check Binding Sites with residues with RSR > %s' % values.rsr_upper), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Bad Ligand', toolTipText='Check Ligands with RSR > %s' % values.rsr_upper), constraints)
        constraints.gridx -= 2
        constraints.gridy += 1
        self.panel.add(JCheckBox('Dubious Binding Site', toolTipText='Check Binding Sites with residues with RSR between %s and %s' % (values.rsr_lower, values.rsr_upper), selected=True), constraints)
        constraints.gridx += 2
        self.panel.add(JCheckBox('Dubious Ligand', toolTipText='Check Ligands with RSR between %s and %s' % (values.rsr_lower, values.rsr_upper), selected=True), constraints)
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
    def __init__(self, args=['--no-args']):
        self.values = argparser.parse_args(args)
        self.panel = JPanel(GridBagLayout())

        if not self.values.distance:
            self.values.distance = math.sqrt(rsr_analysis.inner_distance)
        self.distance  = JTextField(str(self.values.distance))

        if not self.values.outputfile:
            self.values.outputfile = 'rsr_analysis.csv'
        self.outputfile = JTextField(str(self.values.outputfile))

        if not self.values.rsr_lower:
            self.values.rsr_lower = rsr_analysis.RSR_lower
        self.rsr_lower = JTextField(str(self.values.rsr_lower))

        if not self.values.rsr_upper:
            self.values.rsr_upper = rsr_analysis.RSR_upper
        self.rsr_upper = JTextField(str(self.values.rsr_upper))

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
        distancetooltip = 'Residues with at least one atom within this distance from any atom of the ligand will be considered as part of the binding site'
        self.panel.add(JLabel(u'Radius (in Å)', toolTipText=distancetooltip), constraints)
        self.distance.toolTipText=distancetooltip
        constraints.gridx += 1
        self.panel.add(self.distance, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        hrsrtooltip="Ligands and binding sites with at least one 'residue' with an RSR above this value will be tagged as Bad."
        self.panel.add(JLabel('Upper cap for RSR', toolTipText=hrsrtooltip), constraints)
        constraints.gridx += 1
        self.rsr_upper.toolTipText=hrsrtooltip
        self.panel.add(self.rsr_upper, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        lrsrtooltip="Ligands and binding sites with all 'residues' with a RSR below this value will be tagged as Good"
        self.panel.add(JLabel('Good RSR cap', toolTipText=lrsrtooltip), constraints)
        constraints.gridx += 1
        self.rsr_lower.toolTipText=lrsrtooltip
        self.panel.add(self.rsr_lower, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        outtt="Select or enter the name of the output file and where to save it."
        self.panel.add(JButton('Output file name', toolTipText=outtt, actionPerformed=self.selectOutFileName), constraints)
        constraints.gridx += 1
        self.outputfile.toolTipText=outtt
        self.panel.add(self.outputfile, constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        pdbtt="Parse structures from their PDB codes, either by entering their codes or by providing a file containing them"
        self.panel.add(JButton('Use PDB codes', toolTipText=pdbtt, actionPerformed=self.loadStructsFrom), constraints)
        constraints.gridx += 1
        uniprottt="Parse structures from their UniProtKB codes, either by entering their codes or by providing a file containing them"
        self.panel.add(JButton('Use UniProtKB codes', toolTipText=uniprottt, actionPerformed=self.loadStructsFrom), constraints)
        constraints.gridx -= 1

        constraints.gridy += 1
        constraints.gridwidth = 2
        csvfilett = "Load a previously generated file to review its structures"
        self.panel.add(JButton('Load previous results file',  toolTipText=csvfilett, actionPerformed=self.csvFileDialog), constraints)

        self.diag = JDialog(JFrame(iconImage=vhelibsicon),size = (500, 200), title = 'VHELIBS', modal=True)
        self.diag.add(self.panel)
        self.diag.setLocationRelativeTo(None)
        self.diag.pack()

    def selectOutFileName(self, event):
        outfn = str(showOpenDialog(SimpleFileFilter('.csv', None, 'CSV files'), prefkey='loadedFiles', prefs=prefs,multiselect=False))
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

    def isViable(self):
        return bool(self.values.csvfile or self.values.pdbidfile or self.values.pdbids or self.values.swissprot)

    def go(self):
        self.load_settings()
        if self.isViable():
            self.diag.visible = 0

    def getvalues(self):
        if not self.isViable():
            self.diag.visible = 1
            if not self.isViable():
                exit(1)
        return self.values

class WaitDialog(Runnable):
    def __init__(self, parent=None, info=None, modal=False):
        self.info = info if info else '<html>Calculating binding sites and retrieving RSR information<br /> Please be patient</html>'
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
        urv = java.lang.ClassLoader.getSystemClassLoader().getResource(ctnsfn)
        if urv:
            u = URL(str(urv))
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
        constraints.gridwidth = 3
        self.panel.add(JLabel(ImageIcon(vhelibsicon)), constraints)

        label = dedent(u"""<html><body>
        <div style="text-align: center;"><big><big><big><span style="font-weight: bold;">VHELIBS %s</span></big></big></big><br>
        <br>
        <big><span style="font-weight: bold;">V</span>alidation <span style="font-weight: bold;">HE</span>lper for <span style="font-weight: bold;">LI</span>gands and <span style="font-weight: bold;">B</span>inding <span style="font-weight: bold;">S</span>ites</big>
        <br>
        <p style="margin-top: 0.12cm; margin-bottom: 0cm; line-height: 0.53cm;" lang="en-US"><small><small>
        <font face="Helvetica-Light, Arial Unicode MS, sans-serif"><font style="font-size: 13pt;" size="3"><small><small>Adrià
        Cereto-Massagué<sup> 1</sup>, María José Ojeda<sup>1</sup>,
        Cristina Valls<sup> 1</sup>, Miquel Mulero<sup> 1</sup>, M. Josepa
        Salvado<sup>1</sup>, Anna Arola-Arnal<sup>1</sup>, Lluís Arola<sup>1,
        2</sup>, Santiago Garcia-Vallvé<sup>1, 2</sup> and Gerard Pujadas<sup>
        1, 2,</sup></small></small></font></font></small></small></p>
        <p style="margin-top: 0.07cm; margin-bottom: 0.09cm;" lang="pt-BR"><small><small><font face="Helvetica-Light, Arial Unicode MS, sans-serif"><sup>1</sup>Grup
        de Recerca en Nutrigenòmica, Departament de Bioquímica i
        Biotecnologia, Universitat Rovira i Virgili, Campus de Sescelades, C/
        Marce&#320;lí Domingo s/n, 43007 Tarragona, Catalonia, Spain</font></small></small></p>
        <p style="margin-top: 0.07cm; margin-bottom: 0.09cm;" lang="pt-BR"><small><small><font face="Helvetica-Light, Arial Unicode MS, sans-serif"><sup>2</sup>Centre
        Tecnològic de Nutrició i Salut (CTNS), TECNIO, CEICS, Avinguda
        Universitat 1, 43204, Reus, Catalonia, Spain</font></small></small></p>
        <br>
        <br>
        </div>
        </body></html>""") % VHELIBS_VERSION

        constraints.gridy = 1
        self.panel.add(JLabel(label), constraints)

        constraints.gridy = 2
        constraints.gridx = 0
        constraints.gridwidth = 1
        self.panel.add(JLabel(self.urvicon), constraints)

        constraints.gridx = 2
        self.panel.add(JLabel(self.ctnsicon), constraints)

        constraints.gridx = 1
        self.panel.add(JLabel(u"""<html><body>More information, help and documentation can be found at <br><a href="http://urvnutrigenomica-ctns.github.com/VHELIBS/">http://urvnutrigenomica-ctns.github.com/VHELIBS/</a></body></html>"""), constraints)

        self.diag = JDialog(self.frame, title = 'About VHELIBS')
        self.diag.setLocationRelativeTo(self.parentframe)
        self.diag.add(self.panel)
        self.diag.pack()
        self.diag.size = (734,524)
    def show(self):
        self.diag.visible = True

class DisplaySettingsDialog(object):
    keys = ('ligwfv', 'ligsfv', 'ligcolor', 'ligedmcolor', 'bswfv', 'bssfv', 'bscolor', 'bsedmcolor', 'rewfv', 'resfv', 'recolor', 'reedmcolor', 'edmdistance', 'sigma')
    def __init__(self, parent):
        self.parent = parent
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
        self.ligwfv = JTextField()
        self.panel.add(self.ligwfv, constraints)

        constraints.gridy = 2
        self.ligsfv = JTextField()
        self.panel.add(self.ligsfv, constraints)

        constraints.gridy = 3
        self.ligcolor = JTextField()
        self.panel.add(self.ligcolor, constraints)

        constraints.gridy = 4
        self.ligedmcolor = JTextField()
        self.panel.add(self.ligedmcolor, constraints)

        constraints.gridy = 0
        constraints.gridx = 2
        self.panel.add(JLabel('Binding Site'), constraints)

        constraints.gridy = 1
        self.bswfv = JTextField()
        self.panel.add(self.bswfv, constraints)

        constraints.gridy = 2
        self.bssfv = JTextField()
        self.panel.add(self.bssfv, constraints)

        constraints.gridy = 3
        self.bscolor = JTextField()
        self.panel.add(self.bscolor, constraints)

        constraints.gridy = 4
        self.bsedmcolor = JTextField()
        self.panel.add(self.bsedmcolor, constraints)

        constraints.gridy = 0
        constraints.gridx = 3
        self.panel.add(JLabel('Coordinates to exam'), constraints)

        constraints.gridy = 1
        self.rewfv = JTextField()
        self.panel.add(self.rewfv, constraints)

        constraints.gridy = 2
        self.resfv = JTextField()
        self.panel.add(self.resfv, constraints)

        constraints.gridy = 3
        self.recolor = JTextField()
        self.panel.add(self.recolor, constraints)

        constraints.gridy = 4
        self.reedmcolor = JTextField()
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
        self.edmdistance = JTextField()
        self.panel.add(self.edmdistance, constraints)
        constraints.gridx = 2
        self.panel.add(JLabel('EDM sigma', toolTipText="Contour level of the Electron Density Map"), constraints)
        constraints.gridx = 3
        self.sigma = JTextField()
        self.panel.add(self.sigma, constraints)
        constraints.insets = Insets(3,3,3,3)
        constraints.gridy = 6
        constraints.gridx = 1
        constraints.gridwidth = 2
        self.panel.add(JButton('Save', actionPerformed=self.saveprefs), constraints)
        constraints.gridwidth = 1
        constraints.gridx = 3
        self.panel.add(JButton('Defaults', actionPerformed=self.loaddefaults), constraints)
        self.diag = JDialog(self.frame, title = 'Display Settings')
        self.diag.setLocationRelativeTo(self.parent.frame)
        self.diag.add(self.panel)
        self.diag.pack()

    def loadprefs(self):
        self.ligwfv.text = prefs.get('ligwfv', '0.1')
        self.ligsfv.text = prefs.get('ligsfv', '0.2')
        self.ligcolor.text = prefs.get('ligcolor', 'magenta')
        self.ligedmcolor.text = prefs.get('ligedmcolor', 'red')

        self.bswfv.text = prefs.get('bswfv', '0.01')
        self.bssfv.text = prefs.get('bssfv', 'off')
        self.bscolor.text = prefs.get('bscolor', 'white')
        self.bsedmcolor.text = prefs.get('bsedmcolor', 'cyan')

        self.rewfv.text = prefs.get('rewfv', '0.1')
        self.resfv.text = prefs.get('resfv', '0.2')
        self.recolor.text = prefs.get('recolor', 'cpk')
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
    values = SettingsDialog(args).getvalues()
    struva = StruVa(values)
    return struva

if __name__ == '__main__':
    try:
        sv = main()
    except Exception, e:
        print e
        showErrorDialog(e)
