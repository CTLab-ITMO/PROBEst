from setuptools import setup, find_packages
import os

setup(
    name='PROBESt',
    version='0.1.4',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.10',
    author='CTLab',
    author_email='dvsmutin@itmo.ru',
    description='PROBESt: package for nucleotide probes generation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/CTLab-ITMO/PROBESt',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts': [
            'probebase=PROBESt.scripts.probebase:main',
            'genome_operations=PROBESt.genome_operations:main',
        ],
    },
)