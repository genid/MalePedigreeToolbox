from setuptools import setup, find_packages
from pathlib import Path
import re


with open("README.md") as readme:
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


setup(
    name='male_pedigree_toolbox',
    version=get_version(),
    long_description=long_description,
    packages=find_packages(),
    url='https://github.com/bramvanwersch/male_pedigree_toolbox.git',
    entry_points={
        'console_scripts': [
            'mpt=MalePedigreeToolbox.main:main',
            'mpt_gui=MalePedigreeToolbox.gui.main_gui:mpt_gui'
        ],
    },
    license='MIT',
    author='bramv',
    author_email='b.vanwersch@erasmusmc.nl',
    description='tools for getting information from pedigress'
)
