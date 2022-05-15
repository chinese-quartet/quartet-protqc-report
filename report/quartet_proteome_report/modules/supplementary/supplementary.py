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
    )
    
    html = '''
      <!-- Methods -->
      <div class='methods'>
        <div class='small-12 columns'>
        <h3 class='section-header black'>Methods</h3>
        <p>The QC pipeline starts from the expression profiles at peptide/protein levels, and enables to calculate 6 metrics. A Total score is the geometric mean of the linearly normalized values of these metrics.</p>
        <p><b>1. Number of features</b>: We expect as many proteins (mapped to gene symbols) as possible for downstreaming analyses.</p>
        <p><b>2. Missing percentage (%)</b>: Too many missing values interfere with comparability. This metric is calculated globally.</p>
        <p><b>3. Coefficient of variantion (CV, %)</b>: A CV value is calculated to indicate the dispersion within replicates feature by feature.</p>
        <p><b>4. Absolute Correlation</b>: Pearson correlation reflects overall reproducibility within replicates. We calculate correlation coefficients between each two replicates within each biological sample (D5, D6, F7, M8), and take the median as the final value for absolute correlation.</p>
        <p><b>5. Signal-to-Noise Ratio (SNR)</b>: SNR is established to characterize the ability of a platform or lab or batch, which is able to distinguish intrinsic differences among distinct biological sample groups ("signal") from variations in technical replicates of the same sample group ("noise").</p>
        <p><b>6. Relative Correlation with Reference Datasets (RC)</b>: RC is used for assessment of quantitative consistency with the reference dataset at relative levels. For shotgun proteomics, quantitation at peptide levels is theoretically more reliable. Therefore, the reference dataset is established by benchmarking the relative expression values (log2FCs), for each peptide sequence of each sample pair (D5/D6, F7/D6, M8/D6), in historical datasets at peptide levels. We calculate relatively qualified (satisfied with thresholds of p < 0.05) log2FCs of the queried data, for overlapped peptides with the reference dataset, as the input for the assessment of quantitative consistency. Then RC value is Pearson correlation coefficient between the test dataset and the reference dataset.</p>
        </div>
      </div>
      
      <!-- Contact us -->
      <div class='contact'>
        <div class='small-12 columns'>
        <h3 class='section-header black'>Contact us</h3>
          <b>Fudan University Pharmacogenomics Research Center</b>
          <li>Project manager: Quartet Team</li>
          <li>Email: quartet@fudan.edu.cn</li>
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