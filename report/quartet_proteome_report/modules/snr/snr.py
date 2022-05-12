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
import plotly.express as px
import plotly.figure_factory as ff
from quartet_proteome_report.utils.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
  def __init__(self):
        
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
      return None
    
    # Initialise the parent module Class object
    super(MultiqcModule, self).__init__(
      name='Signal-to-Noise Ratio'
    )
    
    snr_pca_df = pd.DataFrame()
    for f in self.find_log_files('snr/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      snr_pca_df = pd.read_csv(f_p, sep = "\t")
    
    # Now add a PCA plot
    if len(snr_pca_df) != 0:
      self.plot_pca("snr-pca", snr_pca_df)
    else:
      log.debug('No file matched: snr - pca_table.txt')
  
  ### Function: Plot the scatter plot
  def plot_pca(self, id, fig_data, title=None, section_name=None, description=None, helptext=None):
    
    fig_data = fig_data[["Sample.ID", "Sample", "PC1", "PC2"]]
    fig_data['PC1'] = fig_data['PC1'].map(lambda x: ('%.3f') % x)
    fig_data['PC2'] = fig_data['PC2'].map(lambda x: ('%.3f') % x)

    fig = px.scatter(fig_data, 
          x = 'PC1', y = 'PC2',
          title = title, 
          color = 'Sample',
          color_discrete_map={"D5": "#00ACC6", "D6": "#5BAF89", "F7": "#FFB132", "M8": "#E8633B"},
          hover_data={'PC1': ':.3f', 'PC2': ':.3f', 'Sample.ID': True},
          render_mode = 'svg')
    
    fig.update_traces(marker=dict(size=15, opacity=1, line_color='white', line_width=0.5))
    fig.update_layout(yaxis_title='PC1',
                      xaxis_title='PC2',
                      font=dict(family="Arial, sans-serif", size=12.5, color="black"),
                      template="plotly_white", 
                      # xaxis_range = [-tick, tick], 
                      # yaxis_range = [-tick, tick],
                      margin=dict(l=150, r=150, t=10, b=10)
                      )
    
    html = plotly_plot(fig, {
          'id': id + '_plot',
          'data_id': id + '_data',
          'title': title,
          'auto_margin': False
          })
    
    # Add a report section with the scatter plot
    self.add_section(
        name="",
        description = """
        SNR is established to characterize the power in discriminating multiple groups. The PCA plot is used to visualise the metric.<br>
        Points are coloured as follows: 
        <span style="color: #00ACC6;"><b>D5</b></span>, 
        <span style="color: #5BAF89;"><b>D6</b></span>, 
        <span style="color: #FFB132;"><b>F7</b></span>, 
        <span style="color: #E8633B;"><b>M8</b></span>.""",
        anchor="correlation-scatter",
        plot = html
    )