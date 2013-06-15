# -*- encoding: utf-8 -*-


__package_name__='AGlib'
__author__="Artur Glavic"
__copyright__="Copyright 2013"
__license__="GPL v3"
__email__="artur.glavic@gmail.com"
__author_email__=__email__
__url__="http://"
__description__='''Python library with general functions often used in scattering.'''
__version__='0.1'


__scripts__=[]
__py_modules__=[]
__package_dir__={}
__packages__=['aglib']
__package_data__={}
__data_files__=[]
__requires__=['numpy', 'scipy']

from distutils.core import setup

#### Run the setup command with the selected parameters ####
setup(name=__package_name__,
      version=__version__,
      description=__description__,
      author=__author__,
      author_email=__email__,
      url=__url__,
      scripts=__scripts__,
      py_modules=__py_modules__,
      packages=__packages__,
      package_dir=__package_dir__,
      package_data=__package_data__,
      data_files=__data_files__,
      requires=__requires__, #does not do anything
     )
