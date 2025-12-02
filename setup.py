from setuptools import setup, find_packages
import os

def read_requirements():
    """Read requirements from requirements.txt file."""
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    with open(requirements_path, 'r') as f:
        requirements = []
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if line and not line.startswith('#'):
                requirements.append(line)
        return requirements

setup(
    name='PROBESt',
    version='0.2.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.12',
    install_requires=read_requirements(),
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
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts': [
            'probebase=PROBESt.scripts.probebase:main',
            'genome_operations=PROBESt.genome_operations:main',
        ],
    },
)