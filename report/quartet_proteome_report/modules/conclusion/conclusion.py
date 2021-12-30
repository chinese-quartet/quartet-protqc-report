#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd

from multiqc import config
from multiqc.plots import table, heatmap
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_proteome_report.modules.plotly import plot as plotly_plot


# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
  def __init__(self):
        
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
      return None
    
    # Initialise the parent module Class object
    super(MultiqcModule, self).__init__(
      name='Assessment Summary'
    )

    ### Conclusion Table
    table_summary = []
    for f in self.find_log_files('conclusion/conclusion_table'):
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


    ### Performance Score
    quality_score_df = pd.DataFrame()
    for f in self.find_log_files('conclusion/rank_table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      quality_score_df = pd.read_csv(f_p, sep = "\t")
      # Sort the dataframe by total score
      quality_score_df.sort_values('Total', inplace=True, ascending=False)
      quality_score_df[['COR']] = quality_score_df[['COR']]*100
      #query_batch_name = 'Lot2_test'
      #indicator = 'Total'
      
      if quality_score_df.shape[0] != 0:
        self.plot_quality_score('plot_quality_score', quality_score_df)
      else:
        log.debug('No file matched: conclusion - rank_table.tsv')


  ### Conclusion Table
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
    
    headers['Performance'] = {
      'title': 'Performance',
      'description': 'Performance',
      'scale': False,
      'format': '{:.0f}',
      "cond_formatting_rules": {
        "green": [{"s_eq": "high"}],
        "lightgreen": [{"s_eq": "mid-high"}],
        "orange": [{"s_eq": "mid-low"}],
        "red": [{"s_eq": "low"}]
        },
      "cond_formatting_colours": [
        {"green": "#0f9115"},
        {"lightgreen": "#70c402"},
        {"orange": "#d97c11"},
        {"red": "#b80d0d"}
        ]
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
      name = 'Evaluation Metrics',
      anchor = id + '_anchor',
      description = """The submitted data to be tested can be divided into 4 levels based on the Quartile Index 
      of the metrics scores by comparing with historical batches: 
        <span style="color: #b80d0d;">low</span>, 
        <span style="color: #d97c11;">mid-low</span>, 
        <span style="color: #70c402;">mid-high</span>, 
        <span style="color: #0f9115;">high</span>.""",
      plot = table.plot(data, headers, table_config)
    )

  ### Plot: quality score
  def plot_quality_score(self, id, quality_score_df, title=None, section_name=None, description=None, helptext=None):
    # After transposing, there are 3 rows and n columns
    final_data = quality_score_df[["SNR_normalized", "COR_normalized", "Total"]].T.values.tolist()
    final_xcats = quality_score_df['batch'].to_list()
    final_ycats = ['Normalized SNR', 'Normalized COR', 'Total Score']

    pconfig = {
      "id": id,
      "xTitle": "Batch",
      "yTitle": "Evaluation Metrics",
      "decimalPlaces": 0,
      "square": False,
      "xcats_samples": False,
    }

    self.add_section(
      name="Performance Score",
      anchor= id + '_anchor',
      description='''
      <p>Scores of evaluation metrics for the current batch and all historical batches assessed.</p>
      <p>Please note that the results shown here are <span style="color: #b80d0d;">Normalized Scores</span> for all batches in each metric.
      The name of your data is <span style="color: #b80d0d;">QUERIED DATA</span>.</p>
      ''',
      plot=heatmap.plot(final_data, final_xcats, final_ycats, pconfig),
    )

