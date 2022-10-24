from setuptools import setup, find_packages

setup(name='btllib',
      version='1.4.9',
      packages=find_packages(include=['btllib']),
      package_data={'btllib': ['_btllib.so']},
      url='https://github.com/bcgsc/btllib')
