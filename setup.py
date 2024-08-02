from setuptools import setup, find_packages
from gfase_maxcut_solver import __version__, __description__, __url__

with open('README.md') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if x]

setup(
    name='gfase_maxcut_solver',
    version=__version__,
    packages=find_packages(),
    author="Jia-Yuan Zhang",
    author_email="zhangjiayuan@genomics.cn",
    description=__description__,
    long_description=long_description,
    url=__url__,
    classifiers=[
        # Trove classifiers
        # Full list at https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
    install_requires=required,
)