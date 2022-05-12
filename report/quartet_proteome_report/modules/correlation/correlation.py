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
from quartet_dnaseq_report.utils.plotly import plot as plotly_plot


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
    self.corr_pca_data = dict()

    corr_df = pd.DataFrame()
    for f in self.find_log_files('correlation/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])

      corr_df = pd.read_csv(f_p, sep = "\t")
      # keys = content.columns.to_list()
      # for index,row in content.iterrows():
      #   corr_pca_table.append(dict(zip(k  eys, row)))
      
      # for i in corr_pca_table:
      #   key = i['Name']
      #   pop_i = i.pop('Name')
      #   self.corr_pca_data[key] = i
    
    # Now add a Scatter plot
    if len(self.corr_df) != 0:
      self.plot_rc("correlation-scatter", corr_df)
    else:
      log.debug('No file matched: correlation - pca_table.tsv')
  
  # def pca_plot(self):
  #   data = OrderedDict()

  #   # cycle over samples and add PC coordinates to data dict
  #   for s_name, d in self.corr_pca_data.items():
  #     data[s_name] = {
  #       "x": d["logFC.Reference"],
  #       "y": d["logFC.Test"],
  #       "color": "rgb(54,100,139,0.5)" # stellblue4
  #     }
    
  #   # generate section and plot
  #   if len(data) > 0:
  #     pconfig = {
  #       "id": "correlation_plot",
  #       "title": "Relative Correlation with Reference Datasets (RC)",
  #       "xlab": "Reference Datasets",
  #       "ylab": "Test Dataset",
  #       "marker_size": 5,
  #       "marker_line_width": 0,
  #       "xmax": 6,
  #       "xmin": -6, 
  #       "ymax": 6,
  #       "ymin": -6, 
  #       "square": True,
  #     }

  #     self.add_section(
  #       name="",
  #       description="",
  #       anchor="correlation-scatter",
  #       plot=scatter.plot(data, pconfig)
  #     )

  ### Function: Plot the scatter plot
  def plot_rc(self, id, fig_data, title=None, section_name=None, description=None, helptext=None):
    
    fig_data['logFC.Test'] = fig_data['logFC.Test'].map(lambda x: ('%.3f') % x)
    fig_data['logFC.Reference'] = fig_data['logFC.Reference'].map(lambda x: ('%.3f') % x)
    
    fig = px.scatter(fig_data, 
          x = 'logFC.Test', y = 'logFC.Reference',
          # symbol = 'Group',
          # symbol_map = {"PCR": 0, "PCR-free": 0, "Queried": 18},
          title = title, 
          color = 'Sample Pair', 
          color_discrete_map={"D5/D6": "#2f5c85", "F7/D6": "#7ba1c7", "M8/D6": "#bb1616"},
          marginal_y='box', marginal_x='box', 
          hover_data={'logFC.Test': ':.3f', 'logFC.Reference': ':.3f', 'Batch': True})
    
    fig.update_traces(marker=dict(size=10, line_color='white', line_width=0.5))
    fig.update_layout(xaxis_title='logFC.Test',
                      yaxis_title='logFC.Reference',
                      font=dict(family="Arial, sans-serif",
                                size=12.5,
                                color="black"),
                      template="simple_white")
    
    html = plotly_plot(fig, {
          'id': id + '_plot',
          'data_id': id + '_data',
          'title': title,
          'auto_margin': True
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
