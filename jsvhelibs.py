import pyodide_http
pyodide_http.patch_all()
import requests
from pyscript import when, display, window, document
from pyweb import pydom
import js
from js import icn3d, cfg, console, createObject
#from pyscript.js_modules import jquery, jquery_ui, three, icn3d
from pyodide.ffi import to_js, create_proxy
createObject(create_proxy(globals()), "python")
import pandas as pd
import math
import rsr_analysis

rsr_analysis.parser.add_argument('--no-args', required=False, action='store_true')

pdb_map_command = "set map 2fofc sigma 0.2 file dsn6 | https://edmaps.rcsb.org/maps/{0}_2fofc.dsn6"
redo_map_command = 'set map 2fofc sigma 0.2 file mtz | https://pdb-redo.eu/db/{0}/{0}_final.mtz'
redo_struc_command = "load url https://pdb-redo.eu/db/{0}/{0}_final.cif | type mmcif"
pdb_struc_command = "load url https://files.rcsb.org/view/{0}.cif | type mmcif"

cfg.showmenu = False
cfg.simplemenu = True
cfg.showcommand = True
cfg.imageonly = False
cfg.closepopup = True

cmd_template = """LOAD_CMD

select all; style proteins nothing; style chemicals nothing; style nucleotides nothing; style ions nothing; style ntbase nothing; style sidec nothing; style water nothing; set surface nothing; glycans cartoon no; set membrane 0; set glycan 0


select {0} | binding_site
saved atoms binding_site
style proteins lines
style nucleotides lines
style ntbase lines
style sidec lines
style chemicals lines
style ions lines
color white
center selection
zoom selection

select {2} | name svligand
saved atoms svligand
style proteins ball and stick
style chemicals ball and stick
style ions ball and stick
style nucleotides ball and stick
style ntbase ball and stick
style sidec ball and stick
color magenta
center selection
zoom selection

select {1} | coords_to_exam
saved atoms coords_to_exam

style proteins ball and stick
style chemicals ball and stick
style ions ball and stick
style nucleotides ball and stick
style ntbase ball and stick
style sidec ball and stick
color b factor percentile

MAP_CMD

select {2} 
center selection
zoom selection
clear all
"""

profiles = {'Default (PDB)':{
                        'distance':math.sqrt(rsr_analysis.INNER_DISTANCE)
                        , 'rsr_lower':rsr_analysis.RSR_lower
                        , 'rsr_upper':rsr_analysis.RSR_upper
                        , 'max_owab':rsr_analysis.OWAB_max
                        , 'min_rscc':rsr_analysis.RSCC_min
                        , 'max_resolution':rsr_analysis.RESOLUTION_max
                        , 'tolerance':rsr_analysis.TOLERANCE
                        , 'min_occupancy':rsr_analysis.OCCUPANCY_min
                        , 'max_rfree':rsr_analysis.RFREE_max
                        , 'use_owab': rsr_analysis.CHECK_OWAB
                        , 'use_res': rsr_analysis.CHECK_RESOLUTION
                        , 'use_pdb_redo': False
                        , 'outputfile':  'vhelibs_analysis_default_PDB.csv'
                        , 'editable': False
                        , 'use_rdiff': rsr_analysis.USE_RDIFF
                        , 'max_rdiff': rsr_analysis.RDIFF_max
                        , 'use_DPI': rsr_analysis.USE_DPI
                        , 'max_DPI': rsr_analysis.DPI_max
                        }
                , 'Default (PDB_REDO)':{
                        'distance':math.sqrt(rsr_analysis.INNER_DISTANCE)
                        , 'rsr_lower': 0.165
                        , 'rsr_upper':rsr_analysis.RSR_upper
                        , 'max_owab':rsr_analysis.OWAB_max
                        , 'min_rscc':rsr_analysis.RSCC_min
                        , 'max_resolution':rsr_analysis.RESOLUTION_max
                        , 'tolerance':rsr_analysis.TOLERANCE
                        , 'min_occupancy':rsr_analysis.OCCUPANCY_min
                        , 'max_rfree':rsr_analysis.RFREE_max
                        , 'use_owab': False
                        , 'use_res': rsr_analysis.CHECK_RESOLUTION
                        , 'use_pdb_redo': True
                        , 'outputfile':  'vhelibs_analysis_default_PDB_REDO.csv'
                        , 'editable': False
                        , 'use_rdiff': rsr_analysis.USE_RDIFF
                        , 'max_rdiff': rsr_analysis.RDIFF_max
                        , 'use_DPI': rsr_analysis.USE_DPI
                        , 'max_DPI': rsr_analysis.DPI_max
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
                        , 'max_rfree':0 #FIXME:
                        , 'use_owab': False
                        , 'use_res': False
                        , 'outputfile':  'vhelibs_analysis_iridium.csv'
                        , 'use_pdb_redo': False
                        , 'editable': False
                        , 'use_rdiff': rsr_analysis.USE_RDIFF
                        , 'max_rdiff': rsr_analysis.RDIFF_max
                        , 'use_DPI': rsr_analysis.USE_DPI
                        , 'max_DPI': rsr_analysis.DPI_max
                        }
                ,'Custom':{
                    'outputfile':  'vhelibs_analysis_custom.csv'
                    , 'editable': True
                    }
                }

def translate_single(res):
    print(res)
    resname = res[:3]
    chain = res[4]
    resnum = res[5:].strip()
    if not resnum.isdigit():
        for char in resnum:
            if not char.isdigit():
                resnum = resnum.replace(char, '')
    ret = ".{}:{}".format(chain, resnum)
    print(ret)
    return ret

def translate_multi(y):
    return " or ".join([translate_single(x) for x in y.split(";") if x.strip()])

class VHELIBS():
    def __init__(self, rows=[]):
        self.rows = rows
        self.pos = 0
        self.cfg = cfg
        #self.pdbid = pdbid
        #print(self.pdbid)
        self.load_pdb = when("click", "#load_pdb")(self.load_pdb)
        self.load_pdb_redo = when("click", "#load_pdb_redo")(self.load_pdb_redo)
        self.next_struct = when("click", "#do_something")(self.next_struct)
        self.analysis = when("click", "#start")(self.analysis)
        if rows:
            self.icn3dui = icn3d.iCn3DUI.new(cfg)
            #await self.next_struct()
        else:
            pass
            print("Calen resultats")
        
        #print(dir(self.icn3dui))
    async def analysis(self, event=None):
        args=['--no-args', '-C']
        values = rsr_analysis.parser.parse_args(args)
        values.pdbids = [p.strip() for p in document.querySelector("#pdbids").value.replace(",", " ").split() if p.strip()]
        values.swissprot = [p.strip() for p in document.querySelector("#upids").value.replace(",", " ").split() if p.strip()]
        values.use_pdb_redo = document.querySelector("#pdbredo").checked
        #values.pdbids = ["3dzu"]
        self.rows, self.badoc, self.rejected = rsr_analysis.main(values)
        self.pos = 0
        await self.next_struct(event)

    async def next_struct(self, event=None):
        print("Next!")
        if self.pos >= len(self.rows):
            display("Finished")
            return
        row = self.rows[self.pos]
        self.pos += 1
        pdbid, residues_to_exam_string, ligandresidues_string, binding_site_string, ligandgood, bsgood, source = row[:7]
        self.pdbid = pdbid
        print(pdbid)
        self.template = cmd_template.format(translate_multi(binding_site_string), translate_multi(residues_to_exam_string), translate_multi(ligandresidues_string))
        if source == "PDB":
            await self.load_pdb(event)
        elif source == "PDB_REDO":
            await self.load_pdb_redo(event)

    async def load(self):
        if hasattr(self, "icn3dui"):
            print("Not re-showing 3D structure")
            #print(self.icn3dui)
            self.icn3dui.icn3d.loadScriptCls.loadScript(self.cfg.command, False)
        else:
            self.icn3dui = icn3d.iCn3DUI.new(self.cfg)
            print("showing 3D structure")
            self.icn3dui.show3DStructure()
            js.viewer = self.icn3dui
        #self.ic = self.icn3dui.icn3d
        
        
    async def load_pdb(self, event=None):
        print("Loading PDB")
        print(self.pdbid)
        pdbid = self.pdbid
        command = self.template.replace("LOAD_CMD", pdb_struc_command).replace("MAP_CMD", pdb_map_command).format(pdbid)
        print(command)
        self.cfg.command = command
        await self.load()
        #await self.ic.loadScriptCls.loadScript(command, False)
        print("done")
    async def load_pdb_redo(self, event=None):
        print("Loading PDB_REDO")
        pdbid = self.pdbid
        command = self.template.replace("LOAD_CMD", redo_struc_command).replace("MAP_CMD", redo_map_command).format(pdbid)
        self.cfg.command = command
        print(command)
        await self.load()
        #await self.ic.loadScriptCls.loadScript(command, False)
        print("done")
    
    

vhelibs = VHELIBS()

@when("click", "#test")
async def do_something(arg):
    print("hola")
    pdbids = document.querySelector("#pdbids").value
    upids = document.querySelector("#upids").value
    pdbredo = document.querySelector("#pdbredo").checked
    print(pdbids)
    print(upids)
    print(pdbredo)
    #await vhelibs.ic.loadScriptCls.loadScript(testcmd, False)

print("Ready")
# requests.get("https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/dz/3dzu.cif.gz")
