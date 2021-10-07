#!/usr/bin/env python

""" Quartet Proteomics Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import plotly.express as px
import plotly.offline as py
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import table
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
      name='Performance Conclusion',
      target='Performance Conclusion',
      #anchor='conclusion',
      #href='https://github.com/clinico-omics/quartet-proteome-report',
      info=' is an report module to show the overall data quality.'
    )
      
    ### Performance Score
    quality_score_df = pd.DataFrame()
    for f in self.find_log_files('conclusion/rank_table'):
      f_p = '%s/%s' % (f['root'], f['fn'])
      quality_score_df = pd.read_csv(f_p, sep = "\t")
      # Sort the dataframe by total score
      #quality_score_df.sort_values('Total', inplace=True, ascending=False)
      query_batch_name = 'Lot2_test'
      indicator = 'Total'

      #test_score = quality_score_df[quality_score_df.batch == 'Lot2_test'].Total.to_list()[0]
      #quality_score_list = quality_score_df['Total'].values.tolist()
      
      if quality_score_df.shape[0] != 0:
        self.plot_quality_score('plot_quality_score', quality_score_df, query_batch_name, indicator)
      else:
        log.debug('No file matched: conclusion - rank_table.tsv')

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


  ### Plot: quality score
  def plot_quality_score(self, id, quality_score_df, query_batch_name, indicator, title=None, section_name=None, description=None, helptext=None):
    quality_score_df.sort_values(indicator, inplace=True, ascending=True)
    performance_score = list(quality_score_df[indicator])
    batch = list(quality_score_df['batch'])
    query_score = list(quality_score_df.loc[quality_score_df['batch'] == query_batch_name][indicator])[0]
    
    fig = px.imshow([performance_score] , x=[performance_score][0], y=['Score'], template="simple_white")

    fig.update_traces(dict(showscale=False,
                      coloraxis=None,
                      colorscale='RdYlGn'),
                      selector={'type': 'heatmap'})

    fig.update_layout(showlegend=False, 
                      annotations=[
                        dict(x=query_score,
                             y=-3,
                             xref="x",
                             yref="y",
                             text=str(query_score)[0:5],
                             textangle=0,
                             showarrow=True,
                             font=dict(family="Arial, sans-serif", size=45, color="black"),
                             align="center",
                             arrowhead=1,
                             arrowsize=4,
                             arrowwidth=3,
                             arrowside="end",
                             arrowcolor="grey",
                             ax=0,
                             ay=-60,
                             yshift=-145)])

    fig.update_xaxes(ticks="outside",
                        tickwidth=2,
                        tickcolor='black',
                        ticklen=10,
                        showline=True,
                        linewidth=2,
                        linecolor='black',
                        tickfont=dict(family='Arial, sans-serif', color='black', size=20))

    fig.update_yaxes(linecolor = 'white', zeroline = False, showline=False, showticklabels=False, showgrid=False, ticks='')
    fig.update_layout(height=500,margin=dict(t=0,b=70))

    html = plotly_plot(fig, {
              'id': id + '_plot',
              'data_id': id + '_data',
              'title': '',
              'auto_margin': True})
    

    # Add a report section with the scatter plot
    self.add_section(
      name='Overall Performance',
      anchor=id + '_anchor',
      description=description if description else
      'Performance metrics and thresholds using Quartet reference metabolite',
      helptext=helptext if helptext else '''
      This longer description explains what exactly the numbers mean
      and supports markdown formatting. This means that we can do _this_:

      * Something important
      * Something else important
      * Best of all - some `code`

      Doesn't matter if this is copied from documentation - makes it
      easier for people to find quickly.
      ''',
      plot=html)


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
      description = description if description else '',
      plot = table.plot(data, headers, table_config)
    )