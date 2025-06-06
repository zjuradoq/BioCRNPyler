from setuptools import setup

# Get the long description from the README file
with open('README.md') as fp:
    long_description = fp.read()

setup(
    name='biocrnpyler',
    version='1.2.0',
    author='Build-A-Cell',
    url='https://github.com/BuildACell/biocrnpyler/',
    description='A chemical reaction network compiler for generating large biological circuit models',
    long_description=long_description,
    packages=['biocrnpyler'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
    ],
    install_requires=[
          "python-libsbml",
          ],
    extras_require = { 
        "all": [
            "numpy<2.0",
            "matplotlib",
            "networkx",
            "bokeh>=1.4.0",
            "fa2_modified",
            "jupyter",
            "pytest",
            "pytest-cov",
            "nbconvert"
            ],
        },
    tests_require=["pytest", "pytest-cov"],
    setup_requires=["pytest-runner"],
    python_requires='>=3.6',
    keywords="SBML synthetic systems biology modeling Chemical Reaction Network CRN model",
    project_urls={
    'Documentation': 'https://readthedocs.org/projects/biocrnpyler/',
    'Funding': 'http://www.cds.caltech.edu/~murray/wiki/index.php?title=Developing_Standardized_Cell-Free_Platforms_for_Rapid_Prototyping_of_Synthetic_Biology_Circuits_and_Pathways',
    'Source': 'https://github.com/BuildACell/biocrnpyler',
    'Tracker': 'https://github.com/BuildACell/BioCRNPyler/issues',
    }, # give an installation message about bioscrape / roadrunner.
)
