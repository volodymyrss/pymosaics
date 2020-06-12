from setuptools import setup
import sys

setup_requires = ['setuptools >= 30.3.3', 'setuptools-git-version']

if {'pytest', 'test', 'ptr'}.intersection(sys.argv):
    setup_requires.append('pytest-runner')


setup(description="pymosaic-fits",
      long_description=open('README.md').read(),
      version='0.4.1',
      include_package_data=True,
      setup_requires=setup_requires)
