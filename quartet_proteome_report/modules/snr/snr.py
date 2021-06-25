#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
from typing import List
import pandas as pd

from multiqc import config
from multiqc.plots import scatter
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
      name='Signal-to-Noise Ratio (SNR)',
      target='SNR',
      #anchor='snr',
      #href='https://github.com/clinico-omics/quartet-proteome-report',
      info=' is established to characterize the power in discriminating multiple groups. The PCA plot is used to visualise the metric.'
    )

    # Find and load any input files for snr
    self.snr_pca_data = dict()
    self.quartet_cats = list()
    self.quartet_colors = {'D5':'#4CC3D9', 'D6':'#7BC8A4', 'F7':'#FFC65D', 'M8':'#F16745'}

    snr_pca_table = []
    for f in self.find_log_files('snr/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
    
    content = pd.read_csv(f_p, sep = "\t")
    self.quartet_cats = list(set(content["group"].to_list()))
    keys = content.columns.to_list()
    for index,row in content.iterrows():
      snr_pca_table.append(dict(zip(keys, row)))
    
    for i in snr_pca_table:
      key = i['sample_id']
      pop_i = i.pop('sample_id')
      self.snr_pca_data[key] = i
    
    # Now add a PCA plot
    if len(self.snr_pca_data) != 0:
      self.pca_plot()
    else:
      log.debug('No file matched: snr - pca_table.txt')

  def pca_plot(self):
    data = OrderedDict()

    # cycle over samples and add PC coordinates to data dict
    for s_name, d in self.snr_pca_data.items():
      if "PC1" in d and "PC2" in d:
        data[s_name] = {
          "x": d["PC1"],
          "y": d["PC2"],
          "color": self.quartet_colors[d["group"]],
        }
    
    # generate section and plot
    if len(data) > 0:
      pconfig = {
        "id": "snr_pca_plot",
        "title": "Principal components of Quartet samples (D5, D6, F7, M8)",
        "xlab": "PC1",
        "ylab": "PC2",
        "marker_size": 8,
        "marker_line_width": 0,
      }

      self.add_section(
        name="",
        description="",
        anchor="snr-pca",
        plot=scatter.plot(data, pconfig)
      )

