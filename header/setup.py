from distutils.core import setup, Extension

whole_module = Extension('_whole', sources=['whole_wrap.cxx'])

setup(name='whole', version='0.1',
      author='Silvester Sabathiel',
      description="""Example SWIG Module.""",
      ext_modules=[whole_module], py_modules=['whole'])