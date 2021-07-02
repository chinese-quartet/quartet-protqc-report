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
            #anchor='data_generation_information',
            #href='https://github.com/clinico-omics/quartet-proteome-report',
            info=' about the proteomics data.'
        )

        # Find and load any input files for data_generation_information
        for f in self.find_log_files('data_generation_information/information'):
            information = eval(f['f'])
        
        if len(information) != 0:
            self.plot_information('data_generation_information', information)
        else:
            log.debug('No file matched: data_generation_information - data_generation_information.txt')

    def plot_information(self, id, data, title='', section_name='', description=None, helptext=None):
        html_data = ["<dl class='dl-horizontal'>"]
        for k,v in data.items():
            line = "        <dt style='text-align:left;margin-top:1ex'>{}</dt><dd>{}</dd>".format(k,v)
            html_data.append(line)
        html_data.append("    </dl>")

        html = '\n'.join(html_data)
        
        self.add_section(
            name = '',
            anchor = '',
            description = '',
            plot = html
        )