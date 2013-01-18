from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

cflags=["-O3"]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("rsr_analysis", ["rsr_analysis.py","rsr_analysis.pxd"],extra_compile_args=cflags)
                  ,Extension("cPdbAtom", ["cPdbAtom.pyx"],extra_compile_args=cflags)
                  ,Extension("PDBfiles", ["PDBfiles.py", "PDBfiles.pxd"],extra_compile_args=cflags)
                  ,Extension("EDS_parser", ["EDS_parser.py", "EDS_parser.pxd"],extra_compile_args=cflags)
                  ,Extension("cofactors", ["cofactors.py","cofactors.pxd"],extra_compile_args=cflags)
                  ,Extension("pdb_redo", ["pdb_redo.py","pdb_redo.pxd"],extra_compile_args=cflags)
                  ]
)
