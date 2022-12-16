#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
from typing import List
import pandas as pd
import math

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
      name='Correlation with Reference Datasets',
    )
    
    # Find and load any input files for correlation
    corr_df = pd.DataFrame()
    for f in self.find_log_files('correlation/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])

      corr_df = pd.read_csv(f_p, sep = "\t")
    
    # Now add a Scatter plot
    if len(corr_df) != 0:
      self.plot_rc("correlation-scatter", corr_df)
    else:
      log.debug('No file matched: correlation - corr_table.tsv')
  
  ### Function: Plot the scatter plot
  def plot_rc(self, id, fig_data, title=None, section_name=None, description=None, helptext=None):
    fig_data = fig_data[['logFC.Test', 'Sample.Pair', 'logFC.Reference', "Sequence"]]
    fig_data.sort_values('Sample.Pair', inplace=True, ascending=True)
    fig_data['logFC.Test'] = fig_data['logFC.Test'].map(lambda x: ('%.3f') % x)
    fig_data['logFC.Reference'] = fig_data['logFC.Reference'].map(lambda x: ('%.3f') % x)
    
    fig_data[['logFC.Test', 'logFC.Reference']] = fig_data[['logFC.Test', 'logFC.Reference']].astype('float')
    min_value = min([fig_data['logFC.Test'].min(), fig_data['logFC.Reference'].min()])
    max_value = max([fig_data['logFC.Test'].max(), fig_data['logFC.Reference'].max()])
    
    tick = max(abs(min_value), abs(max_value))

    fig = px.scatter(fig_data, 
          x = 'logFC.Test', y = 'logFC.Reference',
          title = title, 
          color = 'Sample.Pair',
          color_discrete_map={"D5/D6": "#00ACC6", "F7/D6": "#FFB132", "M8/D6": "#E8633B"},
          hover_data={'logFC.Test': ':.3f', 'logFC.Reference': ':.3f', 'Sequence': True},
          render_mode = 'svg')
    
    fig.update_traces(marker=dict(size=10, opacity=0.5))
    fig.update_layout(yaxis_title='logFC.Test',
                      xaxis_title='logFC.Reference',
                      font=dict(family="Arial, sans-serif", size=12.5, color="black"),
                      template="plotly_white", 
                      xaxis_range = [-tick, tick], 
                      yaxis_range = [-tick, tick],
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
        description="""
        Relative correlation with reference datasets metric which was representing the numerical consistency of the relative expression profiles.
        """,
        anchor="correlation-scatter",
        plot = html
    )
