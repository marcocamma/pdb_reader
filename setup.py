from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='pdb_reader',
      version='0.0.1',
      packages=["pdb_reader",],
      zip_safe=False,
      install_requires=required,
      )
