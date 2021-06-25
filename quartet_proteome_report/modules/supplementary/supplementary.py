#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
import base64
import logging
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

def read_image(image):
  with open(image, "rb") as image_file:
    encoded_string = base64.b64encode(image_file.read())
    return encoded_string.decode('utf-8')

class MultiqcModule(BaseMultiqcModule):
  def __init__(self):
        
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
      return None
    
    # Initialise the parent module Class object
    super(MultiqcModule, self).__init__(
      name='Supplementary',
      target='',
      #anchor='supplementary',
      #href='https://github.com/clinico-omics/quartet-proteome-report',
      #info=''
    )
    
    html = '''
      <!-- Methods -->
      <div class='methods'>
        <div class='small-12 columns'>
        <h3 class='section-header black'>Methods</h3>
        <p>SNR was established to characterize the ability of a platform or lab or batch, which was able to distinguish intrinsic differences among distinct biological sample groups (“signal”) from variations in technical replicates of the same sample group ("noise").</p>

        <p>Relative expression data (fold changes) were calculated for a total of six pairs of sample-to-sample comparisons (D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, and M8/F7) among the Quartet samples. In order to improve the reliability of the reference values, genes that were satisfied with thresholds of p < 0.05 and detectable across the two samples in each sample pair were used.</p>
        </div>
      </div>

      <!-- Disclaimer -->
      <div class='disclaimer'>
        <div class='small-12 columns'>
        <h3 class='section-header black'>Disclaimer</h3>
        <p>This quality control report is only for this specific test data set and doesn’t represent an evaluation of the business level of the sequencing company. This report is only used for scientific research, not for clinical or commercial use. We don’t bear any economic and legal liabilities for any benefits or losses (direct or indirect) from using the results of this report.</p>
        </div>
      </div>
      '''

    self.add_section(
      name = '',
      anchor = '',
      description = '',
      plot = html
    )