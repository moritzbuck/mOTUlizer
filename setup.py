import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mOTUlizer", # Replace with your own username
    version="0.0.3",
    author="Moritz Buck",
    author_email="moritz.buck@slu.se",
    description="making OTUs from genomes, and stats on them. and maybe even core-genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/moritzbuck/0039_mOTUlizer/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords='bioinformatics clustering metagenomics microbial-genomics genomics',
    python_requires='>=3.6',
    scripts = ['mOTUlizer/bin/mOTUlize.py','mOTUlizer/bin/mOTUpan.py']

)
