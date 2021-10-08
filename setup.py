#!/usr/bin/env python
"""
MultiReport for Quartet Proteomics QC
"""

from setuptools import setup, find_packages

version = '0.3.0'

setup(
  name = 'quartet_proteome_report',
  version = version,
  author = 'Yaqing Liu',
  author_email = 'liuyaqing@outlook.com',
  description = 'MultiReport for Quartet Proteomics QC.',
  long_description = __doc__,
  keywords = 'bioinformatics',
  url = 'https://github.com/clinico-omics/quartet-proteome-report',
  download_url = 'https://github.com/clinico-omics/quartet-proteome-report/releases',
  license = 'MIT',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'multiqc==1.9',
    'plotly==4.9.0',
    'pandas==1.1.0'
  ],
  entry_points = {
    'multiqc.modules.v1': [
      'data_generation_information = quartet_proteome_report.modules.data_generation_information:MultiqcModule',
      'conclusion = quartet_proteome_report.modules.conclusion:MultiqcModule',
      'snr = quartet_proteome_report.modules.snr:MultiqcModule',
      'correlation = quartet_proteome_report.modules.correlation:MultiqcModule',
      'supplementary = quartet_proteome_report.modules.supplementary:MultiqcModule'
    ],
    'multiqc.hooks.v1': [
      'execution_start = quartet_proteome_report.custom_code:quartet_proteome_report_execution_start'
    ],
    'multiqc.cli_options.v1': [
      'disable_plugin = quartet_proteome_report.cli:disable_plugin'
    ],
    'multiqc.templates.v1': [
      'report_templates = quartet_proteome_report.templates.default'
    ]
  },
  classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Environment :: Web Environment',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Programming Language :: JavaScript',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Visualization',
  ],
)