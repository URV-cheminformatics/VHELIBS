#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2011 - 2015 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, sys, urllib2, csv, itertools, math, contextlib
try:
    if sys.platform.startswith('java'):
        import java
        #Do appropiate things for jython
        import multithreading as multiprocessing
        #Load optimized java version
        import PdbAtomJava as PdbAtom
    else:
        java = None
        #CPython
        import multiprocessing
        #Use cython optimized version
        from cPdbAtom import PdbAtom
except:
    #Fallback to pure python
    from PdbAtom import PdbAtom

import PDBfiles, EDS_parser, pdb_redo
import cofactors

PDB_REDO = False
CHECK_OWAB = False
OWAB_max = 50
CHECK_RESOLUTION = False
USE_DPI = False
USE_RDIFF = False
RESOLUTION_max = 3.5
RSR_upper = 0.4
RSR_lower = 0.24
RSCC_min = 0.9
RFREE_max = 1
RDIFF_max = 0.05
DPI_max = 0.42
OCCUPANCY_min = 1.0
TOLERANCE = 2
inner_distance = 4.5**2
STATS = False
titles = ['PDB ID', "Coordinates to exam", "Ligand Residues", "Binding Site Residues", "Good Ligand", "Good Binding Site", "Model source", 'Ligand Score', 'Binding Site Score']
stat_titles = ['rFree', 'rWork']
###Create the argument parser###
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--pdbids', nargs='+', default=[], type=str, metavar='PDBID', help='list of PDB ids')
parser.add_argument('-s','--swissprot', nargs='+', default=[], type=str, metavar='SP_AN [or SP_NAME]', help='list of Swiss-Prot protein names or accession numbers')
parser.add_argument('-u','--rsr-upper', type=float, default=RSR_upper, metavar='FLOAT', help='set maximum RSR value for each residue (residues with a higher RSR will be discarded)')
parser.add_argument('-l','--rsr-lower', type=float, default=RSR_lower, metavar='FLOAT', help='set minimum RSR value for each residue (residues with a lower RSR value will be directly considered right)')
parser.add_argument('-b','--max-owab', type=float, default=None, metavar='FLOAT', help='set maximum OWAB (Occupancy-weighted B-factor) per residue')
parser.add_argument('-R','--min-rscc', type=float, default=RSCC_min, metavar='FLOAT', help='set minimum RSCC per residue')
parser.add_argument('-O','--min-occupancy', type=float, default=OCCUPANCY_min, metavar='FLOAT', help='set minimum average occupancy per residue')
parser.add_argument('-F','--max-rfree', type=float, default=RFREE_max, metavar='FLOAT', help='set maximum R-free for the structure')
parser.add_argument('-M','--max-rdiff', type=float, default=None, metavar='FLOAT', help='maximum difference between R and R-free, to avoid overfit models')
parser.add_argument('-D','--max-DPI', type=float, default=None, metavar='FLOAT', help='maximum model DPI')
parser.add_argument('-r','--max-resolution', type=float, default=None, metavar='FLOAT', help='set maximum resolution (in Å) below which to consider Good models')
parser.add_argument('-T','--tolerance', type=int, default=TOLERANCE, metavar='INT', help='set maximum number of non-met criteria of Dubious structures')
parser.add_argument('-d','--distance', type=float, default=math.sqrt(inner_distance), metavar='Å', help='consider part of the binding sites all the residues nearer than this to the ligand (in Å)')
parser.add_argument('-f','--pdbidfile', metavar='PATH', type=unicode, default=None, required=False, help='text file containing a list of PDB ids, one per line')
parser.add_argument('-o','--outputfile', metavar='PATH', type=unicode, default='vhelibs_analysis.csv', required=False, help='output file name')
parser.add_argument('-w','--writeexcludes', metavar='PATH', type=unicode, default=None, required=False, help='Write current excluded HET ids to a file')
parser.add_argument('-e','--excludesfile', metavar='PATH', type=unicode, default=None, required=False, help='Override excluded HET ids with the ones provided in this file')
parser.add_argument('-C','--use-cache', required=False, action='store_true', help="Use cached EDS and PDB data if available for the analysis, otherwise cache it.")
parser.add_argument('-P','--use-pdb-redo', required=False, action='store_true', help="Use models from PDB_REDO instead of PDB.")
parser.add_argument('-V','--verbose', required=False, action='store_true', help="Enable verbose output.")
parser.add_argument('-S','--include-stats', required=False, action='store_true', help="Include model data in output file")
#######################

def dbg(string):
    print(string)

def dummy(*args): pass

SERVICELOCATION="http://www.rcsb.org/pdb/rest/customReport"
columns = [
"experimentalTechnique"
,"rFree"
,"rWork"
,"refinementResolution"
,"unitCellAngleAlpha"
,"unitCellAngleBeta"
,"unitCellAngleGamma"
,"lengthOfUnitCellLatticeA"
,"lengthOfUnitCellLatticeB"
,"lengthOfUnitCellLatticeC"
]
QUERY_TPL = "?pdbids=%s&customReportColumns={}&service=wsfile&format=csv".format(",".join(columns))

def average_occ(residue_atoms):
    return sum([atom.occupancy for atom in residue_atoms])/len(residue_atoms)

def dpi(a, b, c, alpha, beta, gamma, natoms, reflections, rfree):
    dbg("a={}".format(a))
    dbg("b={}".format(b))
    dbg("c={}".format(c))
    dbg("alpha={}".format(alpha))
    dbg("beta={}".format(beta))
    dbg("gamma={}".format(gamma))
    cosa = math.cos(math.radians(alpha))
    dbg("cosa={}".format(cosa))
    cosb = math.cos(math.radians(beta))
    dbg("cosb={}".format(cosb))
    cosg = math.cos(math.radians(gamma))
    dbg("cosg={}".format(cosg))
    V = a*b*c*math.sqrt(1-cosa**2-cosb**2-cosg**2 + 2*cosa*cosb*cosg)
    dbg("V={}".format(V))
    dbg("reflections={}".format(reflections))
    dbg("natoms={}".format(natoms))
    dbg("rfree={}".format(rfree))
    return 1.28*(natoms**(1.0/2))*(V**(1.0/3))*(reflections**(-5.0/6))*rfree

def get_custom_report(pdbids_list):
    urlstring = SERVICELOCATION + QUERY_TPL % ','.join(pdbids_list)
    with contextlib.closing(urllib2.urlopen(urlstring)) as urlhandler:
        if urlhandler.code != 200:
            print urlhandler.msg
            raise IOError(urlhandler.msg)
        reader = csv.reader(urlhandler)
        result = {}
        reader.next()
        for row in reader:
            rowdict = {}
            if row[1] == "X-RAY DIFFRACTION":
                for n, column in enumerate(columns[1:]):
                    n2 = n+2
                    if row[n2] or row[n2] == 0:
                        rowdict[column] = float(row[n2])
            result[row[0]] = rowdict
    return result

def get_sptopdb_dict():
    """
    Returns a dictionary containing the pdb entries for each Swiss-Prot entry
    """
    url = "http://www.uniprot.org/docs/pdbtosp.txt"
    sptopdb_dict = {}
    temppdbdict = {}
    reader = urllib2.urlopen(url)
    pdbid = None
    for line in reader:
        if not len(line) > 28:
            continue
        if line[:28].strip():
            line_stripped = line.strip()
            if line_stripped:
                pdbid_candidate = line_stripped.split()[0]
                if len(pdbid_candidate) == 4 and pdbid_candidate.isupper():
                    pdbid = pdbid_candidate
                    if pdbid not in temppdbdict:
                        temppdbdict[pdbid] = set()
        if pdbid:
            spinfo = line[28:].strip()
            temppdbdict[pdbid].update([item.strip() for item in spinfo.split(',') if item])
    reader.close()
    for pdbid in temppdbdict:
        for sp_id in temppdbdict[pdbid]:
            if sp_id not in sptopdb_dict:
                sptopdb_dict[sp_id] = set()
            sptopdb_dict[sp_id].add(pdbid)
    return sptopdb_dict

def parse_binding_site(argtuple):
    """
    argtuple = (pdbid, )
    """
    pdbid = argtuple[0]
    rsr_upper, rsr_lower = RSR_upper, RSR_lower
    ligand_residues = set()
    good_rsr = set()
    dubious_rsr = set()
    bad_rsr = set()
    res_atom_dict = {}
    ligand_res_atom_dict = {}
    notligands = {}
    links = []
    struc_dict = {}
    reflections = 0
    if PDB_REDO:
        if not argtuple[1]:
            return  (pdbid, "Not in PDB_REDO")
        edd_dict = pdb_redo.get_ED_data(pdbid)
        if not edd_dict:
            return  (pdbid, "Not in PDB_REDO")
        struc_dict['rFree'] = argtuple[1].get('RFFIN',-100)
        struc_dict['rWork'] = argtuple[1].get('RFIN', 1000)
        if USE_DPI:
            reflections = float(argtuple[1].get('NREFCNT', 0))
            a = argtuple[1].get('AAXIS', 0)
            b = argtuple[1].get('BAXIS', 0)
            c = argtuple[1].get('CAXIS', 0)
            alpha = argtuple[1].get('ALPHA', 0)
            beta = argtuple[1].get('BETA', 0)
            gamma = argtuple[1].get('GAMMA', 0)
        resolution = argtuple[1].get('RESOLUTION', 0)
    else:
        if not argtuple[1]:
            dbg("Model not obtained by X-ray crystallography")
            return  (pdbid, "Model not obtained by X-ray crystallography")
        pdbdict, edd_dict = EDS_parser.get_EDS(pdbid)
        if not edd_dict:
            dbg("No EDS data available for %s, it will be discarded" % pdbid)
            return  (pdbid, "No EDS data available")
        struc_dict['rFree'] = argtuple[1].get('rFree', -100)
        struc_dict['rWork'] = argtuple[1].get('rWork', 1000)
        if USE_DPI:
            a = argtuple[1].get('lengthOfUnitCellLatticeA', 0)
            b = argtuple[1].get('lengthOfUnitCellLatticeB', 0)
            c = argtuple[1].get('lengthOfUnitCellLatticeC', 0)
            alpha = argtuple[1].get('unitCellAngleAlpha', 0)
            beta = argtuple[1].get('unitCellAngleBeta', 0)
            gamma = argtuple[1].get('unitCellAngleGamma', 0)
        resolution = argtuple[1].get('refinementResolution', 0) if CHECK_RESOLUTION else 1714
    if CHECK_RESOLUTION:
        struc_dict['Resolution'] = resolution
    if USE_RDIFF:
        Rdiff = struc_dict['rWork'] - struc_dict['rFree']
        struc_dict['Rdiff'] = Rdiff
    pdbfilepath = PDBfiles.get_pdb_file(pdbid.upper(), PDB_REDO)
    natoms = 0
    if pdbfilepath.endswith('.gz'):
        pdbfile = gzip.GzipFile(pdbfilepath)
    else:
        pdbfile = open(pdbfilepath, 'r')
    try:
        error = False
        for line in pdbfile:
            line = line.strip()
            label = line[:6].strip()
            if label in ("ATOM", "HETATM"):
                atom = PdbAtom(line)
                if atom.occupancy == 1.0:
                    natoms += 1
                if label == "ATOM" and inner_distance: #Don't care about protein when distance = 0
                    try:
                        res_atom_dict[atom.residue].add(atom)
                    except KeyError:
                        res_atom_dict[atom.residue] = set([atom, ])
                elif label == 'HETATM':
                    if atom.hetid == 'HOH':
                        continue #Skip waters
                    if (atom.hetid in cofactors.ligand_blacklist) or (atom.hetid in cofactors.metals):
                        try:
                            res_atom_dict[atom.residue].add(atom)
                        except KeyError:
                            res_atom_dict[atom.residue] = set([atom, ])
                        notligands[atom.residue] = "Blacklisted ligand"
                        continue
                    ligand_residues.add(atom.residue)
                    if not atom.residue in ligand_res_atom_dict:
                        ligand_res_atom_dict[atom.residue] = set()
                    ligand_res_atom_dict[atom.residue].add(atom)
            elif label == 'LINK':
                try:
                    dist = float(line[73:78].strip())
                except ValueError:
                    dbg("bogus LINK distance")
                    dist = 1714 #Distance will be calculated when looking for the binding site
                links.append((line[17:27],  line[47:57], float(dist))) #distance
            elif USE_DPI and label == 'REMARK':
                if line[9] == '3' and "NUMBER OF REFLECTIONS" in line:
                    try:
                        reflections = int(line.split(":")[1].strip())
                    except:
                        return  (pdbid, "number of reflections")
    except IOError, error:
        dbg(pdbfilepath)
        dbg(error)
    finally:
        pdbfile.close()
        if error:
            return  (pdbid, str(error))
    if USE_DPI:
        if reflections:
            struc_dict["DPI"] = dpi(a, b, c, alpha, beta, gamma, natoms, reflections, struc_dict["rFree"])
        else:
            struc_dict["DPI"] = 9999
        dbg("DPI is {}".format(struc_dict["DPI"]))
    #Now let's prune covalently bound ligands
    alllinksparsed = False
    while not alllinksparsed:
        for res1,  res2,  blen in list(links):
            checklink = 0
            ligres = sres = None
            if res1 in res_atom_dict:
                checklink +=1
                sres,  ligres = res1, res2
                dbg('Binding to the sequence: %s -> %s' % (ligres, sres))
            if res2 in res_atom_dict:
                checklink +=1
                sres,  ligres = res2, res1
                dbg('Binding to the sequence: %s -> %s' % (ligres, sres))
            if checklink == 2:
                links.remove((res1,  res2,  blen))
                break
            if res1[:3].strip() in cofactors.metals:
                ligres = res2
            elif res2[:3].strip() in cofactors.metals:
                ligres = res1
            if checklink == 1:
                if not blen or blen >= 2.1: #Disulfide bonds are about 2.05;
                    dbg('Bond distance big enough (%s) between %s and %s' % (blen, res1,  res2))
                    continue
                if (res1[:3].strip() in cofactors.metals) or (res2[:3].strip() in cofactors.metals):
                    dbg('Ignoring metal bonds: %s - %s' % (res1, res2))
                    continue
                if (res1[:3].strip() in cofactors.ligand_blacklist) or (res2[:3].strip() in cofactors.ligand_blacklist):
                    notligands[ligres] = "Covalently bound to a blacklisted ligand"
                else:
                    notligands[ligres] = "Covalently bound to the sequence"
                links.remove((res1,  res2,  blen))
                dbg('%s is not a ligand!' % ligres)
                if ligres in ligand_residues:
                    ligand_residues.remove(ligres)
                    dbg('%s removed from ligand residues' % ligres)
                if ligres in ligand_res_atom_dict:
                    res_atom_dict[ligres] = ligand_res_atom_dict.pop(ligres)
                    dbg('%s atoms added to protein' % ligres)
                break
        else:
            alllinksparsed = True

    if not ligand_residues:
        dbg('%s has no ligands!' % pdbid)
        return (pdbid, "no ligands found")
    ligands = group_ligands(ligand_residues, links)
    ligand_scores = []
    for ligand in ligands:
        ligand_score = 0
        for res in list(ligand):
            residue_dict = edd_dict.get(res, None)
            if residue_dict:
                residue_dict['occupancy'] = average_occ(ligand_res_atom_dict[res])
            score,  reason = classificate_residue(res, edd_dict.get(res, None), struc_dict, good_rsr, dubious_rsr, bad_rsr)
            if reason and score >= 1000:
                notligands[res] = reason
                ligand.remove(res)
            ligand_score = max(ligand_score, score)
        ligand_scores.append(ligand_score)
    binding_sites_found = False
    while not binding_sites_found:
        ligand_bs_list = []
        for ligand, ligand_score in zip(ligands, ligand_scores):
            if not ligand:
                continue
            bs = get_binding_site(ligand, ligand_score, good_rsr, bad_rsr, dubious_rsr, pdbid, res_atom_dict, ligands, ligand_res_atom_dict, rsr_upper, rsr_lower, edd_dict, struc_dict, notligands)
            if len(bs) == 1:
                for ligandres in ligand:
                    if ligandres in ligand_residues:
                        ligand_residues.remove(ligandres)
                        dbg('{} removed from ligand residues because {}'.format(ligandres, bs[0]))
                    if ligandres in ligand_res_atom_dict:
                        res_atom_dict[ligandres] = ligand_res_atom_dict.pop(ligandres)
                        dbg('{} atoms added to protein because {}'.format(ligandres, bs[0]))
                    if ligandres not in notligands:
                        notligands[ligandres] = bs[0]
                continue
            ligand_bs_list.append(bs)
        else:
            break
    dbg("done")
    return (pdbid, ligand_bs_list, notligands, struc_dict)

def classificate_residue(residue, residue_dict, struc_dict, good_rsr, dubious_rsr, bad_rsr):
    score = 0
    reason = None
    if not residue_dict:
        score +=1000
        reason = "No data for %s" % residue
        dbg(reason)
    else:
        rscc = residue_dict['RSCC']
        if RSCC_min > rscc:
            score += 1
        if CHECK_OWAB:
            owab = residue_dict['OWAB']
            if not 1 < owab < OWAB_max:
                score +=1
        occ = residue_dict['occupancy']
        if occ > 1:
            score +=1000
            reason = "Occupancy above 1"
        elif occ < OCCUPANCY_min:
            score +=1
        rsr = residue_dict['RSR']
        if rsr > RSR_lower:
            score += 1
            if rsr > RSR_upper:
                score += 1
    if struc_dict:
        rFree = struc_dict['rFree']
        if rFree > RFREE_max:
            score += 1
        if 0 > rFree:
            score +=1000
            reason = "No rFree data for %s" % residue
        if CHECK_RESOLUTION:
            resolution = struc_dict.get('Resolution', 10)
            if resolution > RESOLUTION_max:
                score += 1
            if resolution == 10:
                score +=1000
                reason = "No resolution data for %s" % residue
        if USE_RDIFF:
            Rdiff = struc_dict.get('Rdiff', 100)
            if Rdiff > RDIFF_max:
                score += 1
            elif Rdiff >= 90:
                score +=1000
                reason = "No reliable rFree and rWork data for %s" % residue
        if USE_DPI:
            DPI = struc_dict.get('DPI', -1)
            if not (DPI < DPI_max):
                score += 1
            if 0 > DPI:
                score +=1000
                reason = "No reliable structural data for %s" % residue
    else:
        if USE_DPI or USE_RDIFF or CHECK_RESOLUTION:
            score +=1000
            reason = "No structural data for %s" % residue
            dbg(reason)
    if score == 0:
        good_rsr.add(residue)
    elif score > TOLERANCE:
        bad_rsr.add(residue)
    else:
        dubious_rsr.add(residue)
    dbg("%s: %s" %(residue, score))
    return score, reason

def group_ligands(ligand_residues, links):
    """
    Group all the ligand residues into molecules
    """
    ligands = []

    ligand_links = []
    linked_ligand_res = set()
    for res1,  res2,  blen in links:
        if res1 in ligand_residues and res2 in ligand_residues:
            ligand_links.append((res1,  res2,  blen))
            linked_ligand_res.add(res1)
            linked_ligand_res.add(res2)

    while ligand_links:
        for  res1,  res2,  blen in ligand_links:
            for ligand in ligands:
                n = ''
                if res1 in ligand:
                    n = ligands.index(ligand)
                    ligand.add(res2)
                elif res2 in ligand:
                    n = ligands.index(ligand)
                    ligand.add(res1)
                if n != '':
                    ligands[n] = ligand
                    ligand_links.remove((res1,  res2,  blen))
                    break
            else:
                ligand = set()
                ligand.add(res1)
                ligand.add(res2)
                ligands.append(ligand)
                ligand_links.remove((res1,  res2,  blen))
                break

    for lres in ligand_residues:
        for ligand in ligands:
            #present = False
            if lres in ligand:
                break
        else:
            ligands.append(set([lres, ]))

    all_ligands_parsed = False
    while not all_ligands_parsed:
        all_ligands_parsed = True
        needbreak = False
        for ligand in ligands:
            if needbreak:
                needbreak = False
                break
            for lres in list(ligand):
                if lres in linked_ligand_res:
                    for ligand2 in ligands:
                        if ligand != ligand2:
                            dbg('Checking if %s is in %s'  % (lres, ligand2))
                            if lres in ligand2:
                                ligands.remove(ligand2)
                                n = ligands.index(ligand)
                                ligand.update(ligand2)
                                ligands[n] = ligand
                                all_ligands_parsed = False
                                needbreak = True
                                break
                    else:
                        all_ligands_parsed = True
                else:
                    dbg( '%s is a single ligand' % lres)
                    all_ligands_parsed = True
    return ligands


def get_binding_site(ligand, ligand_score, good_rsr, bad_rsr, dubious_rsr, pdbid, res_atom_dict, ligands, ligand_res_atom_dict, rsr_upper, rsr_lower, edd_dict, struc_dict, notligands):
    """
    Get the binding site residues for the provided ligand and return them in a tuple
    """
    inner_binding_site = set()
    for ligandres in ligand:
        if ligandres in notligands:
            reason = notligands[ligandres]
            return [reason]
        for res in res_atom_dict:
            for atom in res_atom_dict[res]:
                for ligandatom in ligand_res_atom_dict[ligandres]:
                    distance = atom | ligandatom
                    if distance <= inner_distance:
                        if distance < 2.1:
                            hetid = ligandres[:3].strip()
                            if hetid in cofactors.ligand_blacklist:
                                reason = ("Covalently bound to a blacklisted ligand", )
                                dbg('{} is not a ligand! ({})'.format(ligand, reason))
                                notligands[ligandres] = reason
                                return [reason]
                            elif hetid in cofactors.metals:
                                reason = ("Covalently bound to the sequence", )
                                dbg('{} is not a ligand! ({})'.format(ligand, reason))
                                notligands[ligandres] = reason
                                return [reason]
                        inner_binding_site.add(atom.residue)
                        break
        for l in ligands:
            if l == ligand:
                continue
            for lres in l:
                for latom in ligand_res_atom_dict[lres]:
                    for ligandatom in ligand_res_atom_dict[ligandres]:
                        distance = latom | ligandatom
                        if distance <= inner_distance:
                            inner_binding_site.add(lres)
                            break
    bad_occupancy = [ligres for ligres in ligand if edd_dict.get(ligres, {'occupancy':0})['occupancy'] < 1]
    bs_score = 0
    for res in inner_binding_site:
        residue_dict = edd_dict.get(res, None)
        resatoms = res_atom_dict.get(res, None)
        if residue_dict and resatoms:
            if resatoms:
                occ = average_occ(resatoms)
                residue_dict['occupancy'] = occ
                if occ < 1:
                    bad_occupancy.append(res)
            else:
                dbg("No data available for %s!" % res)
                dbg(res) in ligand_res_atom_dict
                residue_dict['occupancy'] = 0
                bad_occupancy.append(res)
        else:
            bad_occupancy.append(res)
        score,  reason = classificate_residue(res, edd_dict.get(res, None), struc_dict, good_rsr, dubious_rsr, bad_rsr)
        bs_score = max(bs_score, score)
    rte = inner_binding_site.union(ligand).difference(good_rsr)
    ligandgood = validate(ligand, good_rsr, bad_rsr, dubious_rsr, pdbid)
    bsgood = validate(inner_binding_site, good_rsr, bad_rsr, dubious_rsr, pdbid)
    return ligand, inner_binding_site, rte, ligandgood, bsgood, bad_occupancy, ligand_score, bs_score

def validate(residues, good_rsr, bad_rsr, dubious_rsr, pdbid):
    if residues <= good_rsr:
        return True
    if residues.intersection(bad_rsr):
        return False
    else:
        for residue in residues:
            if residue in dubious_rsr:
                return 'Dubious'
        else:
            dbg("Unclassified residues for %s:\n" % pdbid)
            dbg(residues)
            dbg(residues.intersection(dubious_rsr))
            dbg(residues.intersection(bad_rsr))
            dbg(residues.intersection(good_rsr))
            dbg('\n')
            return('Dubious')
    return '???'

def results_to_csv(results, outputfile):
    """
    Writes the output of parse_binding_site's to a csv file
    """
    datawritten = False
    csvtitles = titles
    if STATS:
        csvtitles += stat_titles
    with open(outputfile, 'wb') as outfile:
        csvfile = csv.writer(outfile)
        csvfile.writerow(csvtitles)
        with open(os.path.splitext(outputfile)[0] + '_rejected.txt', 'w') as rejectedfile:
            badoc_fn = os.path.splitext(outputfile)[0] + '_low_occupancy.txt'
            with open(badoc_fn, 'w') as badoc_file:
                badoc_w = csv.writer(badoc_file, delimiter="\t")
                badoc_w.writerow(['Binding Site', 'Low Occupancy Residues'])
                erase_badoc = True
                dbg('Calculating...')
                for restuple in results:
                    restuplelen = len(restuple)
                    if  restuplelen == 2:
                        pdbid, reason = restuple
                        rejectedfile.write('%s:\t%s\n' % (pdbid, reason))
                        continue
                    elif restuplelen == 4:
                        pdbid, ligand_bs_list, notligands, struc_dict = restuple
                        for nonligand, reason in notligands.iteritems():
            #                resname = nonligand[:3].strip()
                            line = '{}:\t{} {}\n'.format(pdbid, nonligand, reason)
                            rejectedfile.write( line)
                        for data in ligand_bs_list:
                            ligandresidues, binding_site, residues_to_exam, ligandgood, bsgood, bad_occupancy, ligand_score, bs_score = data
                            id = pdbid
                            if not ligandresidues:
                                dbg('%s has no actual ligands, it will be discarded' % pdbid)
                            else:
                                if PDB_REDO:
                                    source = 'PDB_REDO'
                                else:
                                    source = 'PDB'
                                row = [id, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site), ligandgood, bsgood, source,  ligand_score, bs_score]
                                if STATS:
                                    for k in stat_titles:
                                        row.append(struc_dict[k])
                                csvfile.writerow(row)
                                outfile.flush()
                                datawritten = bool(outputfile)
                                if bad_occupancy:
                                    erase_badoc = False
                                    badoc_w.writerow([id + '|' + list(ligandresidues)[0], ';'.join(bad_occupancy)])
            if erase_badoc and os.path.isfile(badoc_fn):
                os.remove(badoc_fn)
    if not datawritten:
        os.remove(outputfile)
    else:
        print "Results written to {}".format(outputfile)
    return datawritten

def main(values):
    global CHECK_OWAB, OWAB_max, CHECK_RESOLUTION, RESOLUTION_max, RSCC_min, TOLERANCE
    global RSR_upper, RSR_lower, RFREE_max, OCCUPANCY_min, PDB_REDO
    global USE_DPI, USE_RDIFF, RDIFF_max, DPI_max
    global dbg, STATS
    filepath =  values.pdbidfile
    if values.include_stats:
        STATS = True
    if not values.verbose:
        dbg = dummy
    if not values.max_owab is None:
        CHECK_OWAB = True
        OWAB_max = values.max_owab
    if values.use_pdb_redo:
        PDB_REDO = True
        CHECK_OWAB = False
    if not values.rsr_upper > values.rsr_lower:
        dbg('%s is higher than %s!' % (values.rsr_lower, values.rsr_upper))
        raise ValueError
    else:
        RSR_upper = values.rsr_upper
        RSR_lower = values.rsr_lower
    distance=values.distance
    writeexcludes = values.writeexcludes
    excludesfile = values.excludesfile
    TOLERANCE = values.tolerance
    if not values.min_occupancy is None:
        OCCUPANCY_min = values.min_occupancy
    if not values.min_rscc is None:
        RSCC_min = values.min_rscc
    if not values.max_rfree is None:
        RFREE_max = values.max_rfree
    if not values.max_resolution is None:
        CHECK_RESOLUTION = True
        RESOLUTION_max = values.max_resolution
        if STATS:
            stat_titles.append("Resolution")
    if not values.max_rdiff is None:
        USE_RDIFF = True
        RDIFF_max = values.max_rdiff
        if STATS:
            stat_titles.append("Rdiff")
    if not values.max_DPI is None:
        USE_DPI = True
        DPI_max = values.max_DPI
        if STATS:
            stat_titles.append("DPI")
    if values.use_cache:
        dbg('Using cache')
        cachedir = os.path.join(os.path.expanduser('~'), '.vhelibs_cache')
        if not os.path.isdir(cachedir):
            os.makedirs(cachedir)
        PDBfiles.CACHEDIR = cachedir
    if not (values.pdbids or values.swissprot or values.pdbidfile):
        return False
    if distance != None:
        global inner_distance
        inner_distance = distance**2
    pdblist = values.pdbids
    if excludesfile:
        cofactors.load_lists(excludesfile)
        dbg("Loading hetids to exclude from %s" % excludesfile)
    if writeexcludes:
        cofactors.dump_lists(writeexcludes)
        dbg('List of excluded Hetids written to %s' % writeexcludes)
    if values.swissprot:
        sptopdb_dict = get_sptopdb_dict()
        for swissprot_id in values.swissprot:
            for key in sptopdb_dict:
                if swissprot_id in key:
                    pdblist = itertools.chain(pdblist, sptopdb_dict[key])
    if filepath:
        pdblistfile = open(filepath, 'rb')
        pdb_ids_from_file = [line.strip() for line in pdblistfile.read().replace(',', '\n').replace('\t', '\n').split() if line.strip()]
        pdblist = itertools.chain(pdblist, pdb_ids_from_file)
        pdblistfile.close()
    pdblist = [pdbid.lower() for pdbid in pdblist]
    #get_custom_report
    if not PDB_REDO:
        npdbs = len(pdblist)
        maxn = 1000
        if npdbs > maxn:
            pdbids_extra_data_dict = {}
            p = multiprocessing.Pool(multiprocessing.cpu_count())
            for d in p.imap(get_custom_report, [pdblist[i:i+maxn] for i in xrange(0, npdbs, maxn)]):
                pdbids_extra_data_dict.update(d)
        else:
            pdbids_extra_data_dict = get_custom_report(pdblist)
    else:
        pdbids_extra_data_dict = pdb_redo.get_pdbredo_data(pdblist)
    argsarray = [(pdbid, pdbids_extra_data_dict.get(pdbid.upper(), {})) for pdbid in pdblist if pdbid]
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    results = pool.imap(parse_binding_site, argsarray)
    #results = (parse_binding_site(argstuple) for argstuple in argsarray)
    datawritten = results_to_csv(results, values.outputfile)
    pool.terminate()
    pool.join()
    return datawritten
