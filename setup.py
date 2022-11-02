#!/usr/bin/env python

import os
import glob
import shutil
from setuptools import setup, find_packages
from setuptools import Command
from setuptools.command.install import install
from setuptools.command.develop import develop
from distutils.command.build import build
from subprocess import call
from multiprocessing import cpu_count

BASEPATH = os.path.dirname(os.path.abspath(__file__))
MEDCHEM_PATH = os.path.join(BASEPATH, "medchem/lilly")

# define project information
NAME = "medchem"
VERSION = "1.3.8"
DESCRIPTION = "Molecule filtering code for medchem"


def build_make(targets="all", options=["DEBUG=n"]):
    if isinstance(targets, str):
        targets = [targets]
    if isinstance(options, str):
        options = [options]
    cmd = [
        "make",
    ]
    if "all" not in targets:
        try:
            cmd.append("-j%d" % cpu_count())
        except NotImplementedError:
            print("Unable to determine number of CPUs. Using single threaded make.")

    cmd.extend(options)
    cmd.extend(targets)

    def compile():
        call(cmd, cwd=MEDCHEM_PATH)

    return compile


class RunMake(Command):
    description = "Run makefile for C/C++ source"
    user_options = [
        # The format is (long option, short option, description).
        ("target=", "t", "Target for the make"),
    ]

    def initialize_options(self):
        self.target = None

    def finalize_options(self):
        if self.target is None:
            self.target = "clean"
        if self.target not in ["clean", "uninstall", "all"]:
            raise Exception(
                f"Wrong values {self.target} for target. Expect one of ['clean', 'uninstall', 'all']"
            )

    def run(self):
        # run original install code
        # install executable
        self.execute(build_make(targets=self.target), [], "Cleaning C/C++ files")


class MedChemInstall(install):
    def run(self):
        # build lilly medchem
        self.execute(build_make(targets="all"), [], "Building C/C++ files")
        # run original build code
        install.run(self)

    def do_egg_install(self):
        self.execute(build_make(targets="all"), [], "Building C/C++ files")
        install.do_egg_install(self)


class MedChemBuild(build):
    def run(self):
        # run original build code
        build.run(self)
        # build lilly medchem
        build_path = os.path.abspath(MEDCHEM_PATH)

        self.execute(build_make(targets="all"), [], "Compiling lilly medchem")
        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)
        if not self.dry_run:
            for target in ["build", "lib"]:
                target_dir = os.path.join(build_path, target)
                output = os.path.join(self.build_lib, "medchem", "lilly", target)
                self.copy_tree(target_dir, output)


class MedChemDev(develop):
    def run(self):
        # build lilly medchem
        self.execute(build_make(targets="all"), [], "Building C/C++ files")
        # run original build code
        develop.run(self)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name=NAME,
    description="Medchem filtering",
    author="InVivo AI",
    author_email="emmanuel@valencediscovery.com",
    maintainer="Emmanuel Noutahi",
    maintainer_email="emmanuel@valencediscovery.com",
    long_description=read("README.md"),
    packages=find_packages(),
    version=VERSION,
    cmdclass={
        "run_make": RunMake,
        "build": MedChemBuild,
        "develop": MedChemDev,
        "install": MedChemInstall,
        "pyinstall": install,
    },
    license="Not Open Source",
    entry_points={
        "console_scripts": [
            "chemfilter=medchem.cli:run",
        ],
    },
    include_package_data=True,
    python_requires=">=3.6.8",  # Python version restrictions
)
