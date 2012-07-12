# -*- coding: utf-8 -*-
#
#   Copyright 2012 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#

#Python stuff
import sys
######One-jar magic#######
sys.path.append('__pyclasspath__/pylib')
print sys.path
print sys.prefix
if not sys.prefix:
    sys.prefix='.'
######One-jar magic#######
import os
import csv
import math
from sys import exit
#Java stuff
from java.awt import BorderLayout, Dimension, GridLayout, GridBagLayout, GridBagConstraints, Insets
from java.awt.event import ItemEvent
from javax.swing import JFrame, JPanel, JButton, JOptionPane, JTextField, JCheckBox, JLabel, UIManager, JDialog

systemla = UIManager.getSystemLookAndFeelClassName()
_infoicon = UIManager.getIcon("OptionPane.informationIcon")
if systemla == u'javax.swing.plaf.metal.MetalLookAndFeel':
    try:
        UIManager.setLookAndFeel(u'com.sun.java.swing.plaf.gtk.GTKLookAndFeel')
    except Exception,  e:
        print e
else:
    UIManager.setLookAndFeel(systemla)

#Jython-specific stuff
from swingutils.preferences import getUserPrefs, PreferencesAdapter
from swingutils.dialogs.filechooser import showOpenDialog, SimpleFileFilter
from swingutils.dialogs.basic import showErrorDialog, showWarningDialog, showMessageDialog
prefs = getUserPrefs('struva')

#Jmol stuff
from org.jmol.adapter.smarter import SmarterJmolAdapter
from org.jmol.api import JmolViewer
from org.openscience.jmol.app.jmolpanel import AppConsole

#Own stuff
from WrapJOptionPane import JOptionPane2
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
    if string.lower() == u'true':
        return True
    elif string.lower() == u'false':
        return False
    else:
        raise TypeError(string + ' cannot be booleaned!')

###Define useful classes###
class Console(AppConsole):
    def __init__(self, viewer, panel, buttons, parseCommand):
        self.parseCommand = parseCommand
        AppConsole.__init__(self, viewer, panel, buttons)

    def notifyEnabled(self, callbacktype):
        if str(callbacktype) == 'SYNC':
            return True
        return AppConsole.notifyEnabled(self, callbacktype)

    def notifyCallback(self, callbacktype, data):
        if str(callbacktype) == 'SYNC':
            self.parseCommand(data)
        else:
            return AppConsole.notifyCallback(self, callbacktype, data)

class JmolPanel(JPanel):
    def __init__(self, preferredSize):
        self.currentSize = Dimension()
        self.preferredSize = preferredSize
        self.viewer = JmolViewer.allocateViewer(self, SmarterJmolAdapter(), None, None, None, None, None)

    def paint(self, g):
        self.getSize(self.currentSize)
        self.viewer.renderScreenImage(g, self.currentSize.width, self.currentSize.height)

class StruVa(object):
    actions = (u'good', u'bad', u'dubious', u'list', u'help', u'options', u'toggle ligand', u'toggle binding site', u'toggle residues to exam',)
    def __init__(self, values):
        self.actionsDict = {}
        if values:
            self.loadCSV(values)
            self.setupUi()
            self.setVisible(True)
            self.console.sendConsoleMessage(self.helpmsg)
            self.start()

    def setupUi(self):
        frame = JFrame("Structure Validation Helper", defaultCloseOperation = JFrame.EXIT_ON_CLOSE, size = (700, 410))
        self.optionsdiag = OptionsDialog()
        contentPane = frame.contentPane
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
        self.console = Console(jmolPanel.viewer, panel2,"History Variables State Clear Help", self.parseCommand)
        jmolPanel.viewer.jmolCallbackListener = self.console
        panelc.gridwidth = 1
        panelc.gridheight = 2
        panelc.gridx = 2
        panelc.gridy = 0
        panelc.weightx = 0
        panelc.weighty = 0.5
        panelc.fill = GridBagConstraints.VERTICAL
        panel.add(panel2, panelc)
        contentPane.add(panel)
        self.execute('wireframe only')
        self.execute('wireframe off')
        self.execute('set bondMode OR')
        self.execute('set syncScript ON')
        self.execute('set antialiasDisplay ON')
        self.execute('set antialiasTranslucent ON')
        #buttonPanelLayout = GridLayout(3, 2, 3, 3)
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
        buttonPanel = JPanel(GridBagLayout())
        later = []
        for action in self.actions:
            caction = action[0].upper() + action[1:]
            if 'toggle' in action:
                cbpanel.add(JCheckBox(caction, itemStateChanged=self.nextStruct))
                if 'ligand' in action:
                    checked = prefs.get('ligand_edm', False)
                elif 'binding' in action:
                    checked = prefs.get('bindingsite_edm', False)
                elif 'exam' in action:
                    checked = prefs.get('restoexam_edm', True)
                later.append(JCheckBox(caction.replace('Toggle', 'EDM for'), prefbool(checked), itemStateChanged=self.nextStruct))
            self.actionsDict[action] = JButton(caction, actionPerformed=self.nextStruct)
            #buttonPanel.add(self.actionsDict[action], constraints)
            self.execute('function %s () {}' % action.replace(' ', '_'))
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
        #
        constraints.gridx = 0
        constraints.gridy = 0
        buttonPanel.add(self.actionsDict[u'good'], constraints)
        #
        constraints.gridx = 1
        constraints.gridy = 0
        buttonPanel.add(self.actionsDict[u'bad'], constraints)
        #
        constraints.gridx = 2
        constraints.gridy = 0
        buttonPanel.add(self.actionsDict[u'dubious'], constraints)
        #
        constraints.gridx = 0
        constraints.gridy = 1
        buttonPanel.add(self.actionsDict[u'list'], constraints)
        #
        constraints.gridx = 1
        constraints.gridy = 1
        buttonPanel.add(self.actionsDict[u'help'], constraints)
        #
        constraints.gridx = 2
        constraints.gridy = 1
        buttonPanel.add(self.actionsDict[u'options'], constraints)
        #
        buttonPanel.setVisible(True)
        panel2.add(buttonPanel, BorderLayout.NORTH)
        self.panel = panel
        self.setVisible = frame.setVisible

    def parseCommand(self, data):
        command = unicode(data[1].split('##')[0])[:-1].lower()
        if command in self.actions:
            self.nextStruct(text=command)

    def nextStruct(self, event = None, text = None):
        if event:
            text = event.source.text.lower()
        if not text:
            return
        #self.console.sendConsoleEcho("%s has been executed" % text)
        needreload = False
        if self.resultdict:
            ltext = text.lower()
            if ltext in self.resultdict:
                self.key = text.upper()
                self.pdbid = self.key.split('|')[0]
                self.reloadStruct()
            elif 'toggle' in ltext:
                checked = event.getStateChange() == ItemEvent.SELECTED
                if 'ligand' in ltext:
                    prefs['ligand'] = checked
                    self.displayLigand(checked)
                elif 'binding' in ltext or 'site' in ltext:
                    prefs['bindingsite'] = checked
                    self.displayBindingSite(checked)
                elif 'exam' in ltext:
                    prefs['restoexam'] = checked
                    self.displayResToExam(checked)
            elif 'edm for' in ltext:
                checked = event.getStateChange() == ItemEvent.SELECTED
                if 'ligand' in ltext:
                    prefs['ligand_edm'] = checked
                    self.displayLigand(self.actionsDict[u'toggle ligand'].selected)
                elif 'binding' in ltext or 'site' in ltext:
                    prefs['bindingsite_edm'] = checked
                    self.displayBindingSite(self.actionsDict[u'toggle binding site'].selected)
                elif 'exam' in ltext:
                    prefs['restoexam_edm'] = checked
                    self.displayResToExam(self.actionsDict[u'toggle residues to exam'].selected)
            elif ltext == 'list':
                self.console.sendConsoleEcho( '\n'.join(self.resultdict.keys()))
            elif ltext == 'help':
                self.console.sendConsoleMessage(self.helpmsg)
            elif ltext in ('good', 'bad', 'dubious'):
                self.updateOutFile(text)
                self.saveWIP()
                if self.resultdict:
                    self.key = self.resultdict.iterkeys().next()
                    self.pdbid = self.key.split('|')[0]
                    self.reloadStruct()
                else:
                    self.key = None
                    self.pdbid = None
                    self.clean()
            elif ltext.strip() == 'options':
                self.optionsdiag.show()
        else:
            self.clean()

    def clean(self):
        #Neteja-ho tot
        try:
            os.remove(self.wipfilename)
        except:
            pass
        self.execute('delete;')
        self.console.sendConsoleMessage( "No more structures to view!")
        self.console.sendConsoleMessage( 'Do you want to exit now?')
        self.console.sendConsoleMessage( 'you can still manually load more structures and do other operations using the Jmol console')

    def reloadStruct(self):
        self.wd.show()
        #Neteja
        showbs = prefs.get('bindingsite', True)
        showre = prefs.get('restoexam', True)
        showlig = prefs.get('ligand', True)
        self.actionsDict[u'toggle binding site'].selected = False
        self.actionsDict[u'toggle residues to exam'].selected = False
        self.actionsDict[u'toggle ligand'].selected = False
        prefs['bindingsite'] = showbs
        prefs['restoexam'] = showre
        prefs['ligand'] = showlig
        self.execute('delete')
        self.ligandresidues, self.residues_to_exam, self.binding_site = self.resultdict[self.key]
        self.ligandresidues_IS = self.residues_to_exam_IS = self.binding_site_IS = None
        if self.ligandresidues == ['']:
            self.console.sendConsoleMessage( 'Structure without ligands!')
            self.console.sendConsoleMessage( 'Skipping it')
            self.nextStruct(text='Skip')
        else:
            #load in Jmol
            try:
                self.execute('load "=%s"' % self.pdbid)
                self.execute('select all')
                self.execute('wireframe only')
                self.execute('wireframe off')
                self.execute('select none')
                self.actionsDict[u'toggle binding site'].selected = prefbool(prefs['bindingsite'])
                self.actionsDict[u'toggle residues to exam'].selected = prefbool(prefs['restoexam'])
                self.actionsDict[u'toggle ligand'].selected = prefbool(prefs['ligand'])
                self.execute('zoom 0')
            except Exception,  e:
                self.console.sendConsoleMessage("ERROR: " + str(e))
                showErrorDialog(e)
        self.console.sendConsoleEcho( "\n####################################")
        self.console.sendConsoleEcho( "Viewing structure %s" % self.key)
        self.console.sendConsoleEcho( "####################################\n")
        self.wd.show(False)

    def displayBindingSite(self, visible=True):
        if not self.binding_site:
            return
        if not visible:
            self.execute('select(binding_site)')
            self.execute('wireframe only')
            self.execute('wireframe off')
            self.execute('select none')
            if self.binding_site_IS:
                self.execute('isosurface BINDINGSITE off')
                self.binding_site_IS = 0
            return
        binding_site_selection = reslist_to_sel(self.binding_site)
        self.execute('define binding_site (%s)' % binding_site_selection)
        self.execute('select(binding_site)')
        self.execute('wireframe %s' % prefs.get('bswfv','0.01'))
        self.execute('spacefill %s' % prefs.get('bssfv', 'off'))
        self.execute('color %s' % prefs.get('bscolor', 'white'))
        self.execute('select none')
        if prefs.get('bindingsite_edm', False):
            if self.binding_site_IS == 0:
                self.execute('isosurface BINDINGSITE on')
            elif not self.binding_site_IS:
                self.execute('isosurface BINDINGSITE color %s sigma %s within %s {binding_site} "=%s" mesh nofill' %\
                            (prefs.get('bsedmcolor', 'cyan'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.0'),  self.pdbid))
            self.binding_site_IS = 1
        elif self.binding_site_IS:
            self.execute('isosurface BINDINGSITE off')
            self.binding_site_IS = 0

    def displayResToExam(self, visible=True):
        if not self.residues_to_exam:
            return
        if not visible:
            self.execute('select(res_to_exam)')
            self.execute('wireframe only')
            self.execute('wireframe off')
            self.execute('select none')
            if self.residues_to_exam_IS:
                self.execute('isosurface RES_TO_EXAM off')
                self.residues_to_exam_IS = 0
            return
        exam_residues_selection = reslist_to_sel(self.residues_to_exam)
        self.execute('define res_to_exam (%s)' % exam_residues_selection)
        self.execute('select(res_to_exam)' )
        self.execute('wireframe %s' % prefs.get('rewfv', '0.1'))
        self.execute('spacefill %s' % prefs.get('resfv', '0.2'))
        self.execute('color %s' % prefs.get('recolor', 'cpk'))
        self.execute('select none')
        if prefs.get('restoexam_edm', True):
            if self.residues_to_exam_IS == 0:
                self.execute('isosurface RES_TO_EXAM on')
            elif not self.residues_to_exam_IS:
                self.execute('isosurface RES_TO_EXAM color %s sigma %s within %s {res_to_exam} "=%s" mesh nofill' %\
                        (prefs.get('reedmcolor', 'yellow'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.0'),  self.pdbid))
            self.residues_to_exam_IS = 1
        elif self.residues_to_exam_IS:
            self.execute('isosurface RES_TO_EXAM off')
            self.residues_to_exam_IS = 0

    def displayLigand(self, visible=True):
        if not self.ligandresidues:
            return
        if not visible:
            self.execute('select(svligand)')
            self.execute('wireframe only')
            self.execute('wireframe off')
            self.execute('select none')
            if self.ligandresidues_IS:
                self.execute('isosurface LIGAND off')
                self.ligandresidues_IS = 0
            return
        ligands_selection = reslist_to_sel(self.ligandresidues)
        self.execute('define svligand (%s)' % ligands_selection)
        self.execute('select(svligand)')
        self.execute('wireframe %s' % prefs.get('ligwfv', '0.1'))
        self.execute('spacefill %s' % prefs.get('ligsfv', '0.2'))
        self.execute('color %s' % prefs.get('ligcolor', 'magenta'))
        self.execute('select none')
        self.execute('center svligand')
        if prefs.get('ligand_edm', False):
            if self.ligandresidues_IS == 0:
                self.execute('isosurface LIGAND on')
            elif not self.ligandresidues_IS:
                self.execute('isosurface LIGAND color %s sigma %s within %s {svligand} "=%s" mesh nofill' %\
                        (prefs.get('ligedmcolor', 'red'), prefs.get('sigma', '1.0'), prefs.get('edmdistance', '2.0'),  self.pdbid))
            self.ligandresidues_IS = 1
        elif self.ligandresidues_IS:
            self.execute('isosurface LIGAND off')
            self.ligandresidues_IS = 0

    def saveWIP(self):
        outfile = open(self.wipfilename, 'wb')
        csvfile = csv.writer(outfile)
        csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
        for key in self.resultdict:
            ligandresidues, residues_to_exam, binding_site = self.resultdict[key]
            csvfile.writerow([key.split('|')[0], ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
        outfile.close()

    def updateOutFile(self, filetype):
        filetype = filetype.lower()
        if filetype not in self.filesdict:
            raise TypeError('Unknown destination file')
        ligandresidues, residues_to_exam, binding_site = self.resultdict.pop(self.key)
        row = [self.key.split('|')[0], ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)]
        file = open(self.filesdict[filetype], 'ab')
        writer = csv.writer(file, self.dialect)
        writer.writerow(row)
        file.flush()
        file.close()

    def start(self):
        self.key = self.resultdict.iterkeys().next()
        self.pdbid= self.key.split('|')[0]
        self.reloadStruct()

    def loadCSV(self, values):
        self.wd = WaitDialog()
        self.resultdict = {}
        csvfilename = values.csvfile
        if csvfilename:
            if not os.path.isfile(csvfilename):
                print 'File %s does not exist' % csvfilename
                exit(1)
            outdir = os.path.dirname(csvfilename)
            basename = os.path.splitext(os.path.basename(csvfilename))[0]
            _wipfile = os.path.join(outdir, basename + '_wip.csv~')
            if csvfilename.endswith('_wip.csv~'):
                basename = basename.replace('_wip', '')
            elif os.path.isfile(_wipfile):
                ans=JOptionPane.showConfirmDialog(None, "Would you like to continue with the validation \nof the structures from the file %s?" % _wipfile)
                if ans == 2:
                    exit(0)
                elif ans ==0:
                    csvfilename = _wipfile
        else:
            self.wd.show()
            datawritten, goodfilename = rsr_analysis.main(
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
            if goodfilename:
                showMessageDialog('Structures below the specified RSR values were\nconsidered good enough and saved to %s' % goodfilename, 'Good structures saved')
            if datawritten:
                csvfilename = values.outputfile
                outdir = os.path.dirname(csvfilename)
                basename = os.path.splitext(os.path.basename(csvfilename))[0]
                _wipfile = os.path.join(outdir, basename + '_wip.csv~')
            else:
                showWarningDialog('No structures to be viewed.')
                main(['--no-args'])
        print('Loading data from %s...' % csvfilename)
        csvfile = open(csvfilename, 'rb')
        try:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
        except csv.Error,  e:
            print e
            dialect = 'excel'
            print 'using default dialect: %s' % dialect
        csvfile.seek(0)
        self.dialect = dialect
        reader = csv.reader(csvfile, dialect)
        for id, residues_to_exam_string, ligandresidues_string, binding_site_string in reader:
            if id != 'PDB ID':
                residues_to_exam = residues_to_exam_string.split(';')
                ligandresidues = ligandresidues_string.split(';')
                binding_site = binding_site_string.split(';')
                self.resultdict[id + '|' +ligandresidues[0]] = ligandresidues, residues_to_exam, binding_site
        else:
            print 'Data loaded'
        if not self.resultdict:
            print 'File without data!'
            exit(1)
        csvfile.close()
        goodfilename = os.path.join(outdir, basename + '_good.csv')
        badfilename = os.path.join(outdir, basename + '_bad.csv')
        dubiousfilename = os.path.join(outdir, basename + '_dubious.csv')
        self.filesdict = {'good':goodfilename, 'bad':badfilename, 'dubious':dubiousfilename}
        self.wipfilename = os.path.join(outdir, basename + '_wip.csv~')
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


class SettingsDialog(object):
    def __init__(self, values, parent=None):
        self.values = values

        distance = JTextField()
        if not self.values.distance:
            self.values.distance = math.sqrt(rsr_analysis.inner_distance)
        distance.setText(str(self.values.distance))

        outputfile = JTextField()
        if not self.values.outputfile:
            self.values.outputfile = 'rsr_analysis.csv'
        outputfile.setText(str(self.values.outputfile))

        rsr_lower = JTextField()
        if not self.values.rsr_lower:
            self.values.rsr_lower = rsr_analysis.RSR_lower
        rsr_lower.setText(str(self.values.rsr_lower))

        rsr_upper = JTextField()
        if not self.values.rsr_upper:
            self.values.rsr_upper = rsr_analysis.RSR_upper
        rsr_upper.setText(str(self.values.rsr_upper))

        message = ['Distance', distance,  'Output file name', outputfile,  'Lowest RSR', rsr_lower,  'Highest RSR', rsr_upper]

        self.pane = JOptionPane2(message, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION)
        self.dialog = self.pane.createDialog(parent, "Options")
        self.dialog.visible = True
        self.values.distance = float(distance.getText())
        self.values.rsr_lower = float(rsr_lower.getText())
        self.values.rsr_upper = float(rsr_upper.getText())
        self.values.outputfile = outputfile.getText()
        #print name.getText()

class WaitDialog(object):
    def __init__(self):
        i = UIManager.getIcon("OptionPane.informationIcon")
        self.icon = JLabel(i) if i is not None else JLabel(_infoicon)
        self.frame =  JFrame()
        self.panel = JPanel()
        self.panel.add(self.icon)
        self.panel.add(JLabel('<html>Calculating binding sites and retrieving RSR information<br /> Please be patient</html>'))
        self.dialog = JDialog(self.frame,'Please wait', False)
        self.dialog.add(self.panel)
        self.dialog.size = (375, 135)
    def show(self, boolean=True):
         self.dialog.visible = boolean

class OptionsDialog(object):
    def __init__(self):
        self.frame =  JFrame(size = (500, 200), title = 'Options')
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
        self.panel.add(JLabel('Residues to exam'), constraints)

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
        self.panel.add(JLabel('EDM Distance'), constraints)
        constraints.gridx = 1
        self.edmdistance = JTextField()
        self.panel.add(self.edmdistance, constraints)
        constraints.gridx = 2
        self.panel.add(JLabel('EDM sigma'), constraints)
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

        self.frame.add(self.panel)

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

        self.edmdistance.text = prefs.get('edmdistance', '2.0')
        self.sigma.text = prefs.get('sigma', '1.0')

    def saveprefs(self, event=None):
        prefs['ligwfv'] = self.ligwfv.text
        prefs['ligsfv'] = self.ligsfv.text
        prefs['ligcolor'] = self.ligcolor.text
        prefs['ligedmcolor'] = self.ligedmcolor.text

        prefs['bswfv'] = self.bswfv.text
        prefs['bssfv'] = self.bssfv.text
        prefs['bscolor'] = self.bscolor.text
        prefs['bsedmcolor'] = self.bsedmcolor.text

        prefs['rewfv'] = self.rewfv.text
        prefs['resfv'] = self.resfv.text
        prefs['recolor'] = self.recolor.text
        prefs['reedmcolor'] = self.reedmcolor.text

        prefs['edmdistance'] = self.edmdistance.text
        prefs['sigma'] = self.sigma.text

    def show(self, doit=True):
        self.loadprefs()
        self.frame.visible = doit

    def loaddefaults(self, event=None):
        for key in prefs.keys():
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
    values = argparser.parse_args(args)
    while not (values.csvfile or values.pdbidfile or values.pdbids or values.swissprot) :
        options = ['Load CSV file', 'Enter PDB IDs', 'Enter Swissport IDs', 'Tweak options', 'Cancel']
        choice = JOptionPane.showOptionDialog(None, 'Select what to do','Select what to do', JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE, None, options, options[0])
        option = options[choice]
        if option == options[0]:
            csvfilter = SimpleFileFilter('.csv', None, 'CSV files')
            file = showOpenDialog(csvfilter, prefs=prefs, prefkey='loadedFiles', multiselect=False)
            if file:
                values.csvfile = str(file)
        elif option in options[1:3]:
            idstring = JOptionPane.showInputDialog(option)
            if idstring:
                ids = idstring.replace(',', ' ').split()
            else:
                ids = []
            if option == options[1]:
                values.pdbids = ids
            elif option == options[2]:
                values.swissprot = ids
        elif option == options[3]:
            s = SettingsDialog(values)
            values = s.values
        elif option == options[4]:
            exit(0)
    struva = StruVa(values)
    return struva

if __name__ == '__main__':
    sv = main()
