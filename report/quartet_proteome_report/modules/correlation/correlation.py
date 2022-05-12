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
    # corr_data = dict(); corr_pca_table=[]
    corr_df = pd.DataFrame()
    for f in self.find_log_files('correlation/table'):
      f_p = '%s/%s' % (f['root'], f['fn'])

      corr_df = pd.read_csv(f_p, sep = "\t")
      # keys = corr_df.columns.to_list()
      # for index,row in corr_df.iterrows():
      #   corr_pca_table.append(dict(zip(keys, row)))
      
      # for i in corr_pca_table:
      #   key = i['Name']
      #   pop_i = i.pop('Name')
      #   corr_data[key] = i
    
    # Now add a Scatter plot
    if len(corr_df) != 0:
      self.plot_rc("correlation-scatter", corr_df)
      # self.plot_rc1("correlation-scatter", corr_data)
    else:
      log.debug('No file matched: correlation - pca_table.tsv')
  
  # def plot_rc1(self, id, corr_data):
  #   data = OrderedDict()

  #   # cycle over samples and add PC coordinates to data dict
  #   for s_name, d in corr_data.items():
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
  #       "xlab": "logFC.Reference",
  #       "ylab": "logFC.Test",
  #       "marker_size": 5,
  #       "marker_line_width": 0,
  #       "xmax": 9,
  #       "xmin": -9,
  #       "ymax": 9,
  #       "ymin": -9
  #     }
      
  #     self.add_section(
  #       name="",
  #       description="",
  #       anchor="correlation-scatter",
  #       plot=scatter.plot(data, pconfig)
  #     )

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
    print(-tick, tick)

    fig = px.scatter(fig_data, 
          x = 'logFC.Test', y = 'logFC.Reference',
          title = title, 
          color = 'Sample.Pair',
          color_discrete_map={"D5/D6": "#00ACC6", "F7/D6": "#FFB132", "M8/D6": "#E8633B"},
          # color_discrete_map={"D5/D6": "#541F1B", "F7/D6": "#EA6E6C", "M8/D6": "#E9BD84"},
          # color_discrete_map={"D5/D6": "#7470AF", "F7/D6": "#C96728", "M8/D6": "#4B9B7A"},
          # color_discrete_map={"D5/D6": "#0D0881", "F7/D6": "#F1E258", "M8/D6": "#BC5078"},
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
