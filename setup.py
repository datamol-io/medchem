#!/usr/bin/env python

import shutil
import os
from sys import platform

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
MEDCHEM_PATH = os.path.join(BASEPATH, 'medchem/')

# define project information
NAME = 'medchem'
VERSION = '1.0.0'
DESCRIPTION = 'Molecule filtering code for medchem'

class LillyBuild(build):
    def run(self):
        # run original build code
        build.run(self)

        # build lilly medchem
        build_path = os.path.abspath(self.build_temp)

        cmd = [
            'make',
            'OUT=' + build_path,
        ]

        try:
            cmd.append('-j%d' % cpu_count())
        except NotImplementedError:
            print 'Unable to determine number of CPUs. Using single threaded make.'

        options = [
            'DEBUG=n',
            'ENABLE_SDL=n',
        ]
        cmd.extend(options)

        targets = ['python']
        cmd.extend(targets)

        if platform == 'darwin':
            target_path = 'OSX64_PYTHON'
        else:
            target_path = 'UNIX_PYTHON'

        target_files = [
            os.path.join(build_path, target_path, 'bin', 'xcsoar.so')
        ]

        def compile():
            call(cmd, cwd=XCSOAR_PATH)

        self.execute(compile, [], 'Compiling xcsoar')

        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, self.build_lib)


class XCSoarInstall(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)

        # install XCSoar executables
        self.copy_tree(self.build_lib, self.install_lib)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
