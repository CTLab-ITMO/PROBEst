from setuptools import setup, find_packages
import os

def read_requirements():
    """Read requirements from requirements.txt."""
    with open('requirements.txt') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='PROBESt',
    version='0.1.4',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.8',
    install_requires=read_requirements(),
    author='Your Name',
    author_email='your.email@example.com',
    description='PROBESt: package for nucleotide probes generation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/PROBESt',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
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