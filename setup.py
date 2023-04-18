"""Setup script"""
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
import pytest


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

install_requires = [
    'MOODS-python'
]

tests_require = [
    'pytest >= 7.2.0, < 8',
    'coverage',
    'pylint == 2.15.5',
    'mypy == 0.982',
]


class PyTest(TestCommand):
    """Run tests"""

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name="MiniMotif",
    version="0.0.1",
    author="Hannah Augustijn",
    author_email="hannah.augustijn@wur.nl",
    description="MiniMotif: a tool for transcription factor binding site detection in bacteria",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/HAugustijn/MiniMotif",
    license='GNU Affero General Public License v3.0',

    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3.0',
        'Operating System :: Linux or MacOS',
    ],
    python_requires='>=3.7',
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data={'minimotif': ['data/*.joblib', 'data/*.sh'],},
    install_requires=install_requires,
    cmdclass={'test': PyTest},
    tests_require=tests_require,
    extras_require={'testing': tests_require,},
)