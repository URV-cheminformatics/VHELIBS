from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("rsr_analysis", ["rsr_analysis.py","rsr_analysis.pxd"],extra_compile_args=["-O3"])
                  ,Extension("cPdbAtom", ["cPdbAtom.pyx"],extra_compile_args=["-O3"])
                  ]
)
