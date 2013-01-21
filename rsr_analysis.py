#/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2011 - 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, gzip, sys, urllib2, csv, itertools, math
try:
    if sys.platform.startswith('java'):
        import java
        #Do appropiate things for jython
        import multithreading as multiprocessing
        #Load optimized java version
        import PdbAtomJava as PdbAtom
        from sys import exit
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
RESOLUTION_max = 3.5
RSR_upper = 0.4
RSR_lower = 0.24
RSCC_min = 0
RFREE_min = 0
OCCUPANCY_min = 1.0
TOLERANCE = 1
inner_distance = 4.5**2
titles = ['PDB ID', "Coordinates to exam", "Ligand Residues", "Binding Site Residues", "Good Ligand", "Good Binding Site", "Model source"]

###Create the argument parser###
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--pdbids', nargs='+', default=[], type=str, metavar='PDBID', help='list of PDB ids')
parser.add_argument('-s','--swissprot', nargs='+', default=[], type=str, metavar='SP_AN [or SP_NAME]', help='list of Swiss-Prot protein names or accession numbers')
parser.add_argument('-u','--rsr_upper', type=float, default=RSR_upper, metavar='FLOAT', help='set maximum RSR value for each residue (residues with a higher RSR will be discarded)')
parser.add_argument('-l','--rsr_lower', type=float, default=RSR_lower, metavar='FLOAT', help='set minimum RSR value for each residue (residues with a lower RSR value will be directly considered right)')
parser.add_argument('-b','--max_owab', type=float, default=None, metavar='FLOAT', help='set maximum OWAB (Occupancy-weighted B-factor) per residue')
parser.add_argument('-R','--min_rscc', type=float, default=RSCC_min, metavar='FLOAT', help='set minimum RSCC per residue')
parser.add_argument('-O','--min_occupancy', type=float, default=OCCUPANCY_min, metavar='FLOAT', help='set minimum average occupancy per residue')
parser.add_argument('-F','--min_rfree', type=float, default=RFREE_min, metavar='FLOAT', help='set minimum R-free for the structure')
parser.add_argument('-r','--max_resolution', type=float, default=None, metavar='FLOAT', help='set maximum resolution (in Å) below which to consider Good models')
parser.add_argument('-T','--tolerance', type=int, default=TOLERANCE, metavar='INT', help='set maximum number of non-met criteria of Dubious structures')
parser.add_argument('-d','--distance', type=float, default=math.sqrt(inner_distance), metavar='Å', help='consider part of the binding sites all the residues nearer than this to the ligand (in Å)')
parser.add_argument('-f','--pdbidfile', metavar='PATH', type=unicode, default=None, required=False, help='text file containing a list of PDB ids, one per line')
parser.add_argument('-o','--outputfile', metavar='PATH', type=unicode, default='vhelibs_analysis.csv', required=False, help='output file name')
parser.add_argument('-w','--writeexcludes', metavar='PATH', type=unicode, default=None, required=False, help='Write current excluded HET ids to a file')
parser.add_argument('-e','--excludesfile', metavar='PATH', type=unicode, default=None, required=False, help='Override excluded HET ids with the ones provided in this file')
parser.add_argument('-C','--use-cache', required=False, action='store_true', help="Use cached EDS and PDB data if available for the analysis, otherwise cache it.")
parser.add_argument('-P','--use-pdb-redo', required=False, action='store_true', help="Use models from PDB_REDO instead of PDB.")
#######################

def dbg(string):
    print(string)
    return 0

SERVICELOCATION="http://www.rcsb.org/pdb/rest/customReport"
QUERY_TPL = "?pdbids=%s&customReportColumns=experimentalTechnique,rFree&service=wsfile&format=csv"

def average_occ(residue_atoms):
    return sum([atom.occupancy for atom in residue_atoms])/len(residue_atoms)

def get_custom_report(pdbids_list):
    urlstring = SERVICELOCATION + QUERY_TPL % ','.join(pdbids_list)
    urlhandler = urllib2.urlopen(urlstring)
    if urlhandler.code != 200:
        print urlhandler.msg
        raise IOError(urlhandler.msg)
    reader = csv.reader(urlhandler)
    result = {}
    header = reader.next()
    rowlen = len(header)
    for row in reader:
        rowdict = {}
        if row[1] == "X-RAY DIFFRACTION":
            rowdict['rFree'] = float(row[2]) if row[2] else 0
        result[row[0]] = rowdict
    urlhandler.close()
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
    if PDB_REDO:
        if not argtuple[1]:
            return  (pdbid, "Not in PDB_REDO")
        edd_dict = pdb_redo.get_ED_data(pdbid)
        if not edd_dict:
            return  (pdbid, "Not in PDB_REDO")
        edd_dict['rFree'] = float(argtuple[1].get('RFFIN', '0'))
        resolution = float(argtuple[1].get('RESOLUTION', '0'))
        edd_dict['Resolution'] = resolution
    else:
        if not argtuple[1]:
            dbg("Model not obtained by X-ray crystallography")
            return  (pdbid, "Model not obtained by X-ray crystallography")
        pdbdict, edd_dict = EDS_parser.get_EDS(pdbid)
        if not pdbdict[pdbid.lower()]:
            dbg("No EDS data available for %s, it will be discarded" % pdbid)
            return  (pdbid, "No EDS data available")
        edd_dict['rFree'] = argtuple[1].get('rFree', 0)
        resolution = 0 if CHECK_RESOLUTION else 1714
#    for key, value in argtuple[1].items():
#        edd_dict[key] = value
    pdbfilepath = PDBfiles.get_pdb_file(pdbid.upper(), PDB_REDO)
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
            elif resolution == 0 and label == 'REMARK':
                if line[9] == '2' and len(line) > 10:
                    try:
                        resolution = float(line[23:30])
                        edd_dict['Resolution'] = resolution
                    except ValueError:
                        return  (pdbid, "resolution not found")
    except IOError, error:
        dbg(pdbfilepath)
        dbg(error)
    finally:
        pdbfile.close()
        if error:
            return  (pdbid, str(error))
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
    ligands_res = set()
#    for ligand in ligands:
#        ligands_res.update(ligand)
#    ligdiff = ligand_residues.difference(ligands_res)
#    if ligdiff:
#        dbg("!!!Ligand residues without ligand:")
#        dbg("\n".join(ligdiff))
    for ligand in ligands:
        for res in ligand:
            residue_dict = edd_dict.get(res, None)
            if residue_dict:
                residue_dict['occupancy'] = average_occ(ligand_res_atom_dict[res])
            if 1337 == classificate_residue(res, edd_dict, good_rsr, dubious_rsr, bad_rsr):
                notligands[res] = "Occupancy above 1"

    binding_sites_found = False
    while not binding_sites_found:
        ligand_bs_list = []
        for ligand in ligands:
            bs = get_binding_site(ligand, good_rsr, bad_rsr, dubious_rsr, pdbid, res_atom_dict, ligands, ligand_res_atom_dict, rsr_upper, rsr_lower, edd_dict)
            if len(bs) == 1:
                notligands[ligand] = bs[0]
                dbg('%s is not a ligand!' % ligand)
                if ligand in ligand_residues:
                    ligand_residues.remove(ligand)
                    dbg('%s removed from ligand residues' % ligand)
                if ligand in ligand_res_atom_dict:
                    res_atom_dict[ligand] = ligand_res_atom_dict.pop(ligand)
                    dbg('%s atoms added to protein' % ligand)
                break
            ligand_bs_list.append(bs)
        else:
            break
    dbg("done")
    return (pdbid, ligand_bs_list, notligands)

def classificate_residue(residue, edd_dict, good_rsr, dubious_rsr, bad_rsr):
    score = 0
    residue_dict = edd_dict.get(residue, None)
    if not residue_dict:
        bad_rsr.add(residue)
        return 0
    rscc = residue_dict['RSCC']
    if RSCC_min > rscc:
        score -= 1
    if CHECK_OWAB:
        owab = residue_dict['OWAB']
        if not 1 < owab < OWAB_max:
            score -=1
    occ = residue_dict['occupancy']
    if occ > 1:
        bad_rsr.add(residue)
        return 1337
    elif occ < OCCUPANCY_min:
        score -=1
    rsr = residue_dict['RSR']
    rFree = edd_dict['rFree']
    if rFree < RFREE_min:
        score -= 1
    if CHECK_RESOLUTION:
        resolution = edd_dict.get('Resolution', 0)
        if resolution > RESOLUTION_max:
            score -= 1
    if rsr <= RSR_upper:
        if rsr <= RSR_lower:
            score +=1
    else:
        score -= 1
    if score > 0:
        good_rsr.add(residue)
    elif score <= -TOLERANCE:
        bad_rsr.add(residue)
    else:
        dubious_rsr.add(residue)
    return 0

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
            present = False
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


def get_binding_site(ligand, good_rsr, bad_rsr, dubious_rsr, pdbid, res_atom_dict, ligands, ligand_res_atom_dict, rsr_upper, rsr_lower, edd_dict):
    """
    Get the binding site residues for the provided ligand and return them in a tuple
    """
    inner_binding_site = set()
    for ligandres in ligand:
        for res in res_atom_dict:
            for atom in res_atom_dict[res]:
                for ligandatom in ligand_res_atom_dict[ligandres]:
                    distance = atom | ligandatom
                    if distance <= inner_distance:
                        if distance < 2.1:
                            hetid = ligandres[:3].strip()
                            if hetid in cofactors.ligand_blacklist:
                                reason = ("Covalently bound to a blacklisted ligand", )
                                return reason
                            elif hetid in cofactors.metals:
                                reason = ("Covalently bound to the sequence", )
                                return reason
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
                            #classificate_residue(lres, edd_dict, good_rsr, dubious_rsr, bad_rsr)
                            break
    for res in inner_binding_site:
        residue_dict = edd_dict.get(res, None)
        resatoms = res_atom_dict.get(res, None)
        if residue_dict and resatoms:
            if resatoms:
                residue_dict['occupancy'] = average_occ(resatoms)
            else:
                dbg("No data available for %s!" % res)
                print res in ligand_res_atom_dict
                residue_dict['occupancy'] = 0
        classificate_residue(res, edd_dict, good_rsr, dubious_rsr, bad_rsr)
    rte = inner_binding_site.union(ligand).difference(good_rsr)
    ligandgood = validate(ligand, good_rsr, bad_rsr, dubious_rsr, pdbid)
    bsgood = validate(inner_binding_site, good_rsr, bad_rsr, dubious_rsr, pdbid)
    return ligand, inner_binding_site, rte, ligandgood, bsgood

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
    outfile = open(outputfile, 'wb')
    rejectedfile = open(os.path.splitext(outputfile)[0] + '_rejected.txt', 'w')
    csvfile = csv.writer(outfile)
    csvfile.writerow(titles)
    dbg('Calculating...')
    datawritten = False
    for restuple in results:
        restuplelen = len(restuple)
        if  restuplelen == 2:
            pdbid, reason = restuple
            rejectedfile.write('%s:\t%s\n' % (pdbid, reason))
            continue
        elif restuplelen == 3:
            pdbid, ligand_bs_list, notligands = restuple
            for nonligand in notligands:
                resname = nonligand[:3].strip()
                line = '%s:\t%s %s\n' %(pdbid, nonligand, notligands[nonligand])
                rejectedfile.write( line)
            for ligandresidues, binding_site, residues_to_exam, ligandgood, bsgood in ligand_bs_list:
                id = pdbid
                if not ligandresidues:
                    dbg('%s has no actual ligands, it will be discarded' % pdbid)
                else:
                    if PDB_REDO:
                        source = 'PDB_REDO'
                    else:
                        source = 'PDB'
                    csvfile.writerow([id, ';'.join(residues_to_exam), ';'.join(ligandresidues),';'.join(binding_site), ligandgood, bsgood, source])
                    outfile.flush()
                    datawritten = bool(outputfile)
    outfile.close()
    rejectedfile.close()
    if not datawritten:
        os.remove(outputfile)
    return datawritten

def main(values):
    global CHECK_OWAB, OWAB_max, CHECK_RESOLUTION, RESOLUTION_max, RSCC_min, TOLERANCE
    global RSR_upper, RSR_lower, RFREE_min, OCCUPANCY_min, PDB_REDO
    filepath =  values.pdbidfile
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
    if not values.min_rfree is None:
        RFREE_min = values.min_rfree
    if not values.max_resolution is None:
        CHECK_RESOLUTION = True
        RESOLUTION_max = values.max_resolution
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
    pdblist = list(pdblist)
    #get_custom_report
    if not PDB_REDO:
        npdbs = len(pdblist)
        maxn = 1595
        if npdbs > maxn:
            pdbids_extra_data_dict = {}
            steps = npdbs / maxn
            remnant = npdbs % maxn
            p = multiprocessing.Pool()
            for d in p.imap(get_custom_report, (pdblist[i:i+maxn] for i in xrange(0, npdbs, maxn))):
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

if __name__ == "__main__":
    r = parse_binding_site(('3dzu', 0.4, 0.24))
    print r
