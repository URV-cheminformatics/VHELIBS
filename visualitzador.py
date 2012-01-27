#-*- coding: utf-8 -*-

import sys
import os
import csv
import urllib2
import tarfile
from java.awt import Frame,  Panel, BorderLayout, FlowLayout, Label, TextField
from java.lang import System
import astex.MoleculeViewer as MoleculeViewer
#print 'astex importat'

PDBbase = "http://www.rcsb.org/pdb/files/%s.pdb.gz"
edsurl = "http://eds.bmc.uu.se/eds/dfs/PDB2/PDB1/PDB1.tar.gz"

datadir = 'data'
if not os.path.isdir(datadir):
    os.mkdir(datadir)

def reslist_to_sel(reslist):
    sellist = []
    for res in reslist:
        resname = res[:3]
        chain = res[4]
        resnum = res[5:].strip()
        sellist.append('(chain %s and residue %s)' % (chain, resnum))
    return ' or '.join(sellist)

def main():
    if not len(sys.argv) or os.path.abspath(__file__) == os.path.abspath(sys.argv[-1]):
        print "Falta un argument: el fitxer de resultats"
        return
    print 'Carregant dades...',
    resultdict = {}
    csvfilename = sys.argv[-1]
    if not os.path.isfile(csvfilename):
        print 'El fitxer %s no existeix' % csvfilename
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
            if ans.lower().strip() in ('sí', 'si', 'yes', 'ja', 'da', 'bai', 'oui', 'oc', 'òc','jes', 'yeah', 'sim', 'ok', 'Oook'):
                csvfilename = basename + '_wip.csv'
                answered = True
            elif ans.lower().strip() in ('no', 'non', 'nein', 'nope', 'ez', 'ne', 'pas', 'não', 'nao', 'Eeek'):
                answered = True
            else:
                ans = raw_input()
    csvfile = open(csvfilename, 'rb')
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
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
        exit()
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
    moleculeViewer.moleculeRenderer.execute('set selectcount off; set symmetry off;')

    #moleculeViewer.moleculeRenderer.execute('molecule load pdb_files/P00491/2A0W.pdb.gz pdb_files/P00491/2A0W.pdb.gz;')


    ###
    if not os.path.isdir('PDB'):
        os.mkdir('PDB')
    pdbid = None
    needreload = False
    while resultdict:
        if not pdbid:
            pdbid = resultdict.keys()[0]
            needreload = True
        if needreload:
            #Netejar
            moleculeViewer.moleculeRenderer.execute('molecule remove *;')
            moleculeViewer.moleculeRenderer.execute('map remove *;')
            #Descarregar pdb
            localpdb = os.path.join(datadir, pdbid.lower(), pdbid + '.pdb.gz')
            if not os.path.isfile(localpdb):
                if not os.path.isdir(os.path.join(datadir, pdbid.lower())):
                    os.mkdir(os.path.join(datadir, pdbid.lower()))
                url = PDBbase % pdbid
                print 'Descarregant %s' % url
                remotepdb = urllib2.urlopen(url)
                localfile = open(localpdb,'wb')
                localfile.write(remotepdb.read())
                localfile.close()
                remotepdb.close()
            relmap = os.path.join(pdbid.lower(), pdbid.lower() + '.omap')
            localmap = os.path.join(datadir, relmap)
            #Visualitzar a l'astex
            moleculeViewer.moleculeRenderer.execute('molecule load %s %s;' % (pdbid, localpdb))
            moleculeViewer.moleculeRenderer.repaint()
            exam_residues_selection = []
            binding_site_selection = []
            ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
            ligands_selection = reslist_to_sel(ligandresidues)
            binding_site_selection = reslist_to_sel(binding_site)
            exam_residues_selection = reslist_to_sel(residues_to_exam)
            moleculeViewer.moleculeRenderer.execute('display lines off all;')
            moleculeViewer.moleculeRenderer.execute('display sticks on %s;' % exam_residues_selection)
            moleculeViewer.moleculeRenderer.execute('color cyan %s;' % exam_residues_selection)
            moleculeViewer.moleculeRenderer.execute('display lines on %s;' % binding_site_selection)
            moleculeViewer.moleculeRenderer.execute('display spheres on %s;' % ligands_selection)
            moleculeViewer.moleculeRenderer.execute('color magenta %s;' % ligands_selection)
            moleculeViewer.moleculeRenderer.repaint()
            moleculeViewer.moleculeRenderer.execute('select %s;' % exam_residues_selection)
            selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
            moleculeRenderer.setCenter(selectedAtoms)
            moleculeViewer.moleculeRenderer.repaint()

            try:
                if not os.path.isfile(localmap):
                    genurl = 'http://eds.bmc.uu.se/cgi-bin/eds/gen_zip.pl?pdbCode=' + pdbid.lower()
                    generate = urllib2.urlopen(genurl)
                    generate.read()
                    generate.close()
                    url = edsurl.replace('PDB1', pdbid.lower()).replace('PDB2', pdbid[1:3].lower())
                    print 'Descarregant %s' % url
                    remotearchive = urllib2.urlopen(url)
                    localarchive = os.path.join(datadir, pdbid.lower(), os.path.basename(url))
                    archivefile = open(localarchive, 'wb')
                    archivefile.write(remotearchive.read())
                    remotearchive.close()
                    archivefile.close()
                    archivefile = tarfile.open(localarchive)
                    archivefile.extract(relmap, datadir)
                moleculeViewer.moleculeRenderer.execute('map load %s %s;' % (pdbid, localmap))
                moleculeViewer.moleculeRenderer.execute('select %s or %s or %s;' % exam_residues_selection, ligands_selection, binding_site_selection)
                selectedAtoms = moleculeRenderer.getSelectedOrLabelledAtoms()
                moleculeRenderer.clipMaps(None, selectedAtoms, True)
                moleculeViewer.moleculeRenderer.repaint()
            except Exception,  e:
                print 'Impossible carregar el mapa de densitat electronica per %s' % pdbid
            needreload = False
        print "\n####################################"
        print "Visualitzant l'estructura %s" % pdbid
        print "####################################\n"
        print "Com fer servir l'OpenAstexViewer:"
        print "http://openastexviewer.net/web/interface.html"
        #Descarregar estructura
        #Carregar estructura
        #Fer visible només els lligands i els residus d'interès
        print '> Possibles ordres: '
        print "good : Considera l'estructura com a bona, la desa a %s" % goodfilename
        print "bad : Considera l'estructura com a incorrecta, la desa a %s" % badfilename
        print "dubious : Considera l'estructura com a dubtosa, la desa a %s" % dubiousfilename
        print "<pdbid> : carrega l'estructura <pdbid>"
        print "list : mostra la llista d'estructures"
        print ">>>",
        inp = raw_input()
        #inp = 'bad'
        print
        if inp.upper() in resultdict:
            pdbid = inp.upper()
            needreload = True
        elif inp.lower() == 'list':
            print resultdict.keys()
        elif inp.lower() in ('good', 'bad', 'dubious'):
            ligandresidues, residues_to_exam, binding_site = resultdict.pop(pdbid)
            writerdict[inp.lower()].writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
            ###
            outfile = open(basename + '_wip.csv', 'wb')
            csvfile = csv.writer(outfile)
            csvfile.writerow(['PDB ID', "Residues to exam", "Ligand Residues", "Binding Site Residues"])
            for pdbid in resultdict:
                ligandresidues, residues_to_exam, binding_site = resultdict[pdbid]
                csvfile.writerow([pdbid, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site)])
            outfile.close()
            pdbid = None
        else:
            print "### Ordre o codi PDB no reconeguts: '%s'. Torna-ho a provar ###" % inp

       #Netejar-ho tot

    print "Enhorabona! S'han acabat les estructures!"
