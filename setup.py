from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    "pymummer>=0.10.3",
    'biopython>=1.70',
]

desc = """Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). Compose synthetic MITES from complete DNA transposons."""

setup(name='TE-splitter',
      version='0.1.0',
      description=desc,
      url='https://github.com/Adamtaranto/TE-splitter',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['tesplitter'],
      classifiers=pypi_classifiers,
      keywords=["Transposon","LTR","TIR","MITE","TE","retrotransposon"],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'tsplit-LTR=tesplitter.cmd_LTR:main',
            'tsplit-TIR=tesplitter.cmd_TIR:main',
        ],
    },
    )
