#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
  def __init__(self):
        
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
      return None
    
    # Initialise the parent module Class object
    super(MultiqcModule, self).__init__(
      name='Performance Conclusion',
      target='Performance Conclusion',
      #anchor='conclusion',
      #href='https://github.com/clinico-omics/quartet-proteome-report',
      info=' is an report module to show the overall data quality.'
    )

    # Find and load any input files for conclusion
    table_summary = []
    for f in self.find_log_files('conclusion/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
    
    content = pd.read_csv(f_p, sep = "\t")
    keys = content.columns.to_list()
    for index,row in content.iterrows():
      table_summary.append(dict(zip(keys, row)))
    
    table_summary_dic = {}
    for i in table_summary:
      key = i['Quality Metrics']
      pop_i = i.pop('Quality Metrics')
      table_summary_dic[key] = i

    if len(table_summary_dic) != 0:
      self.plot_summary_table('conclusion_summary', table_summary_dic)
    else:
      log.debug('No file matched: conclusion - conclusion_table.tsv')
  
  def plot_summary_table(self, id, data, title='', section_name='', description=None):
    """ Create the HTML for Performance Conclusion """
    
    headers = OrderedDict()
    headers['Value'] = {
      'title': 'Value',
      'description': 'Value',
      'scale': False,
      'format': '{0:.2f}'
    }

    headers['Historical value(mean ± SD)'] = {
      'title': 'Historical Value',
      'description': 'Historical value(mean ± SD)',
      'scale': False,
      'format': '{0:.2f}'
    }

    headers['Rank'] = {
      'title': 'Rank',
      'description': 'Rank',
      'scale': False,
      'format': '{:.0f}'
    }
    
    table_config = {
      'namespace': 'conclusion_summary',
      'id': id,
      'table_title': '',
      'col1_header': 'Quality Metrics',
      'no_beeswarm': True,
      'sortRows': False
    }

    # Add a report section with the table
    self.add_section(
      name = section_name if section_name else '',
      anchor = id + '_anchor',
      description = description if description else '',
      plot = table.plot(data, headers, table_config)
    )