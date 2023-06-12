#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
from string import printable
import pandas as pd
import numpy as np
import os

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
    # Add to self.css and self.js to be included in template
    self.css = {
        "assets/css/rank.css": os.path.join(
            os.path.dirname(__file__), "assets", "css", "rank.css"
        )
    }

    ### Cutoff Table
    cutoff_table = pd.DataFrame()
    for f in self.find_log_files('conclusion/cutoff_table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      cutoff_table = pd.read_csv(f_p, sep = "\t")
    if cutoff_table.empty:
      log.debug('No file matched: conclusion - cutoff_table.tsv')


    ### Conclusion Table
    table_summary_dic = {}
    for f in self.find_log_files('conclusion/conclusion_table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      content = pd.read_csv(f_p, sep = "\t").replace(r'\\u00b1', '±', regex=True)
      content.columns = ['Quality Metrics', 'Value', 'Historical Value (mean ± SD)', 'Rank', 'Performance']
      table_summary_dic = content.set_index('Quality Metrics').T.to_dict()
    if len(table_summary_dic) != 0:
      self.plot_summary_table('conclusion_summary', table_summary_dic, cutoff_table)
    else:
      log.debug('No file matched: conclusion - conclusion_table.tsv')


    ### Performance Score
    quality_score_df = pd.DataFrame()
    for f in self.find_log_files('conclusion/rank_table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      quality_score_df = pd.read_csv(f_p, sep = "\t")
      # Sort the dataframe by total score
      quality_score_df.sort_values('Total', inplace=True, ascending=True)
    if quality_score_df.shape[0] != 0:
      self.plot_quality_score('plot_quality_score', quality_score_df)
    else:
      log.debug('No file matched: conclusion - rank_table.tsv')


  ### Function 1: Evaluation metrics
  def plot_summary_table(self, id, table_data, score_bar, title='', section_name='', description=None, helptext=None):
    """ Create the HTML for Performance Conclusion """
    Q1 = format(score_bar.loc[1, 'Percentile'], '.2f')
    Q2 = format(score_bar.loc[2, 'Percentile'], '.2f')
    Q3 = format(score_bar.loc[3, 'Percentile'], '.2f')
    # Calculate percentage
    total = table_data['Total Score']['Rank'].split('/')[1]
    bad_len = 20; bad = "%.2f%s" % (bad_len, '%')
    fair_len = 30; fair = "%.2f%s" % (fair_len, '%')
    good_len = 30; good = "%.2f%s" % (good_len, '%')
    great_len = 20; great = "%.2f%s" % (great_len, '%')
    # Queried data arrow and score
    # score_bar = score_bar.loc[5, 'Cut-off']
    score = format(score_bar.loc[5, 'Percentile'], '.3f')
    queried = "%.2f%s" % (float(score_bar.loc[5, 'Cut-off'].strip('%'))*2, '%') #(((total-float(score_bar.strip('%')))*2/total + 1/total) * 100, '%')
    if float(score) <= 1:
      queried = "0%"
      score = 1
    elif float(score) >= 10:
      queried = "200%"
      score = 10
    # Position of ticks
    tick_Q1 = "%.2f%s" % (bad_len-0.9, '%')
    tick_Q2 = "%.2f%s" % (bad_len+fair_len-1, '%')
    tick_Q3 = "%.2f%s" % (bad_len+fair_len+good_len-1, '%')
    # print(tick_Q1, tick_Q2, tick_Q3)
    overview_html = """
    <!-- Arrow -->
    <div class="arrow" style="width: {queried}; margin-top:10px; height: 35px;">
      <svg class="lower-tangle" transform="translate(0 18)"></svg>
      <span class="lower-label" style="margin-bottom: 25px;"><b> {score} </b></span>
    </div>
    
    <!-- Progress bar -->
    <div class="progress">
      <div class="progress-bar progress-bar-bad" style="width: {bad}" data-toggle="tooltip" title="" data-original-title=""><b>Bad</b></div>
      <div class="progress-bar progress-bar-fair" style="width: {fair}" data-toggle="tooltip" title="" data-original-title=""><b>Fair</b></div>
      <div class="progress-bar progress-bar-good" style="width: {good}" data-toggle="tooltip" title="" data-original-title=""><b>Good</b></div>
      <div class="progress-bar progress-bar-great" style="width: {great}" data-toggle="tooltip" title="" data-original-title=""><b>Great</b></div>
    </div>
    
    <!-- Scale interval -->
    <span style="float:left; left:0%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">1</span>
    <span style="float:left; left:{tick_Q1}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q1}</span>
    <span style="float:left; left:{tick_Q2}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q2}</span>
    <span style="float:left; left:{tick_Q3}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q3}</span>
    <span style="float:left; left:99%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">10</span>
    <br>
    """.format(queried=queried, score=score, bad=bad, fair=fair, good=good, great=great, Q1=Q1, Q2=Q2, Q3=Q3, tick_Q1=tick_Q1, tick_Q2=tick_Q2, tick_Q3=tick_Q3)
    
    headers = OrderedDict()
    headers['Value'] = {
      'title': 'Value',
      'description': 'Value',
      'scale': False,
      'format': '{0:.3f}'
    }

    headers['Historical Value (mean ± SD)'] = {
      'title': 'Historical Value',
      'description': 'Historical value(mean ± SD)',
      'scale': False,
      'format': '{0:.3f}'
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
        "green": [{"s_eq": "Great"}],
        "lightgreen": [{"s_eq": "Good"}],
        "orange": [{"s_eq": "Fair"}],
        "red": [{"s_eq": "Bad"}]
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
    
    table_html = table.plot(table_data, headers, table_config)

    # Add a report section with the table
    self.add_section(
      name = 'Evaluation Metrics',
      anchor = id + '_anchor',
      description = """
      The performance of the submitted data will be graded as <span style="color: #b80d0d;font-weight:bold">Bad</span>, <span style="color: #d97c11;font-weight:bold">Fair</span>, <span style="color: #70c402;font-weight:bold">Good</span>, or <span style="color: #0f9115;font-weight:bold">Great</span> based on the ranking by comparing the total score with the historical datasets.<br>
      The total score is the geometric mean of the scaled values of the number of features, missing percentage, absolute correlation, coefficient of variantion, Signal-to-Noise Ratio (SNR), and relative correlation with reference datasets (RC).
      """,
      plot = overview_html + '\n' + table_html,
      helptext = helptext if helptext else '''
      **Evaluation metrics:**
      
      * The total score is the geometric mean of the scaled values of the number of features, missing percentage, absolute correlation, coefficient of variantion, Signal-to-Noise Ratio (SNR), and relative correlation with reference datasets (RC).
      * For better comparison and presentation, the total score was scaled to the interval [1, 10], with the worst dataset being 1 and the best dataset scoring 10.

      **Four levels of performance:**
      
      Based on the scaled total score, the submitted data will be ranked together with all Quartet historical datasets. The higher the score, the higher the ranking. After this, the performance levels will be assigned based on their ranking ranges.

      * _Bad_ - the bottom 20%.
      * _Fair_ - between bottom 20% and median 50%.
      * _Good_ - between median 50% and top 20%.
      * _Great_ - the top 20%.
      '''
    )
  
  ### Function 2: Historical scores
  def plot_quality_score(self, id, quality_score_df, title=None, section_name=None, description=None, helptext=None):
    quality_score_df = quality_score_df.replace('QUERIED DATA', 'Queried_Data')#.replace(np.nan, 0)
    quality_score_df = quality_score_df.drop(columns=['Total'])

    metrics = ["Batch", "Absolute Correlation", "Coefficient of variantion (CV, %)", "Missing percentage (%)", "Number of features", "Relative Correlation with Reference Datasets (RC)", "Signal-to-Noise Ratio (SNR)", "Total_norm"]
    metrics[-1] = 'Total Score'
    quality_score_df.columns = metrics
    final_data = quality_score_df[metrics[1:]].T.values.tolist()
    final_xcats = quality_score_df['Batch'].to_list()
    final_ycats = metrics[1:]
    
    pconfig = {
      "id": id,
      "xTitle": "Performance of batches gradually increases from left to right",
      "yTitle": "Evaluation metrics",
      "decimalPlaces": 3,
      "square": False,
      "xcats_samples": False,
      "reverseColors": True,
      "height": 251,
      "borderWidth": 0
    }
    
    self.add_section(
      name="Performance Score",
      anchor= id + '_anchor',
      description='''
      Scores of evaluation metrics for the current batch and all historical batches assessed.<br>Please note that the results shown here are <span style="background-color: transparent;font-weight:bold;">scaled values</span> for all batches in each metric. The name of your data is <span style="background-color: transparent;font-weight:bold;">Queried_Data</span>. The white colour block means a <span style="background-color: transparent;font-weight:bold;">NULL</span> value, although 0.00 is shown.
      ''',
      plot=heatmap.plot(final_data, final_xcats, final_ycats, pconfig),
    )

