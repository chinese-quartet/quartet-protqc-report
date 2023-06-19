#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
import logging
from multiqc import config
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
      name='Data Generation Information',
      target='The basic information',
      info=' about the proteomics data.'
    )
    
    # Find and load any input files for general_information
    information = []
    for f in self.find_log_files('general_information/information'):
      information = eval(f['f'])
    
    if len(information) != 0:
      self.plot_information('general_information', information)
    else:
      log.debug('No file matched: general_information - general_information.txt')

  def plot_information(self, id, data, title='', section_name='', description=None, helptext=None):
    html_data = ["<dl class='dl-horizontal'>"]
    for k,v in data.items():
      line = " <dt style='text-align:left; width: 250px'>{}</dt>\n <dd>{}</dd>".format(k,v)
      html_data.append(line)
    html_data.append("</dl>")

    html = '\n'.join(html_data)
    
    self.add_section(
      name = '',
      anchor = '',
      description = '',
      plot = html
    )