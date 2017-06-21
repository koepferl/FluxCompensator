#from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name='FluxCompensator',
    version='0.1.0',
    author='Christine Maria Koepferl and Thomas Robitaille',
    author_email='koepferl@usm.lmu.de',
    #packages=['fluxcompensator', 'fluxcompensator.utils', 'fluxcompensator.database', 'fluxcompensator.database.extinction', 'fluxcompensator.database.fieldstars', 'fluxcompensator.database.missions', 'fluxcompensator.database.missions.PSF', 'fluxcompensator.database.missions.filter'],
    packages=find_packages(),
    package_data={'fluxcompensator.database.missions.PSF': ['*.fits'],
                  'fluxcompensator.database.missions.filter': ['*.txt'],
                  'fluxcompensator.database.extinction': ['*.txt'],
                  'fluxcompensator.database.fieldstars': ['*.txt'],
                  'fluxcompensator.tests': ['*.txt', '*.fits', '*.rtout '],
                  'fluxcompensator.database.fieldstars': ['*.txt']},
    #url='http://pypi.python.org/pypi/FluxCompensator/',
    install_requires = ['numpy>=1.10.4',
                        'astropy>=1.1.1',
                        'matplotlib>=1.5.1',
                        'scipy>=0.17.0',
                        'pytest>=2.9.2',
                        'hyperion>=0.9.7',
                        'ipython'],
    license='LICENSE.txt',
    description='Create "realistic" synthetic observations from continuum radiative transfer simulations making those directly comparable to real observations.',)
