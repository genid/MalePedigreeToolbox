from setuptools import setup, find_packages
from pathlib import Path
import re


with open("README.rst") as readme:
    long_description = readme.read()


def get_version():
    """Get version number from __init__.py"""
    version_file = Path(__file__).resolve().parent / "MalePedigreeToolbox" / "__init__.py"
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read_text(), re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


def get_requirements():
    requirements = []
    with open('requirements.txt') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements


setup(
    name='male_pedigree_toolbox',
    version=get_version(),
    long_description=long_description,
    packages=find_packages(),
    url='https://github.com/genid/MalePedigreeToolbox.git',
    entry_points={
        'console_scripts': [
            'mpt=MalePedigreeToolbox.main:main',
            'mpt_gui=MalePedigreeToolbox.gui.main_gui:mpt_gui'
        ],
    },
    license='MIT',
    author='Bram van Wersch and Diego Montiel Gonzalez',
    author_email='b.vanwersch@erasmusmc.nl',
    description='tools for getting information from pedigress',
    install_requires=get_requirements(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
