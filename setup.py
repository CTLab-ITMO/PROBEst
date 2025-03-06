from setuptools import setup, find_packages
import json
import os

def read_pipenv_dependencies(fname):
    """Pipfile.lock dependencies"""
    filepath = os.path.join(os.path.dirname(__file__), fname)
    with open(filepath) as lockfile:
        lockjson = json.load(lockfile)
        return [dependency for dependency in lockjson.get('default')]

if __name__ == '__main__':
    setup(
        name='PROBESt',
        version=os.getenv('PACKAGE_VERSION', '0.1.3'),
        package_dir={'': 'src'},
        packages=find_packages('src', include=[
            'PROBESt*'
        ]),
        description='PROBESt: package for nucleotide probes generation',
        install_requires=[
              *read_pipenv_dependencies('Pipfile.lock'),
        ]
    )