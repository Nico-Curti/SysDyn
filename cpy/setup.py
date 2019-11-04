# -*- coding: utf-8 -*-
#!/usr/bin/env python

# python setup.py develop --user
# cython Brusselator.pyx -a

try:
  from setuptools import setup
  from setuptools import Extension
  from setuptools import find_packages

except ImportError:
  from distutils.core import setup
  from distutils.core import Extension
  from distutils.core import find_packages

from Cython.Distutils import build_ext
from distutils.sysconfig import customize_compiler

class _build_ext (build_ext):
  '''
  Custom build type
  '''

  def build_extensions (self):

    customize_compiler(self.compiler)

    try:
      self.compiler.compiler_so.remove('-Wstrict-prototypes')

    except (AttributeError, ValueError):
      pass

    build_ext.build_extensions(self)


ext_modules = [ Extension(name='fastBrusselator',
                          sources=['Brusselator.pyx'],
                          libraries=['m'],
                          extra_compile_args = ['-g0',
                                                '-Ofast',
                                                '-Wno-unused-function'],
                          language='c++'
                          )]

setup(
        name='Brusselator',
        cmdclass = {'build_ext': _build_ext},
        ext_modules = ext_modules
      )
