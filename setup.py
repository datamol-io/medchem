#!/usr/bin/env python

import os
import glob

try:
    from setuptools import setup
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install

from distutils.command.build import build
from subprocess import call
from multiprocessing import cpu_count

BASEPATH = os.path.dirname(os.path.abspath(__file__))
MEDCHEM_PATH = os.path.join(BASEPATH, 'medchem')

# define project information
NAME = 'medchem'
VERSION = '1.0.0'
DESCRIPTION = 'Molecule filtering code for medchem'

class LillyBuild(build):
    def run(self):
        # run original build code
        build.run(self)

        # build lilly medchem
        build_path = os.path.abspath(MEDCHEM_PATH)

        cmd = [
            'make',
        ]

        try:
            cmd.append('-j%d' % cpu_count())
        except NotImplementedError:
            print('Unable to determine number of CPUs. Using single threaded make.')

        options = [
            'DEBUG=n',
        ]
        cmd.extend(options)

        targets = ['all']
        cmd.extend(targets)

        target_files = list(glob.glob(os.path.join(build_path, 'build', '*'))) + list(glob.glob(os.path.join(build_path, 'lilly', 'lib', '*')))

        def compile():
            call(cmd, cwd=MEDCHEM_PATH)

        self.execute(compile, [], 'Compiling lilly medchem')

        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, self.build_lib)


class LillyInstall(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)

        # install executable
        self.copy_tree(self.build_lib, self.install_lib)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name=NAME,
    version=VERSION,
    description='LillyMedchem filtering rules',
    author='InVivo AI',
    author_email='emmanuel@invivoai.com',
    install_requires=["numpy"],
    maintainer='Emmanuel Noutahi',
    maintainer_email='emmanuel@invivoai.com',
    long_description=read('README.md'),
    packages=['medchem'],
    license='Not Open Source',
    scripts=glob.glob('bin/*'),
    include_package_data=True,
    python_requires=">=3.6.8",  # Python version restrictions
    cmdclass={
        'build': LillyBuild,
        'install': LillyInstall,
    }
)