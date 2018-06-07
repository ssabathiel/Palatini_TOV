from distutils.core import setup, Extension

whole_module = Extension('_test_main', sources=['test_main_wrap.cxx','test_main.cpp'])

setup(name='test_main', version='0.1',
      author='Silvester Sabathiel',
      description="""Example SWIG Module.""",
      ext_modules=[whole_module], py_modules=['test_main'])
