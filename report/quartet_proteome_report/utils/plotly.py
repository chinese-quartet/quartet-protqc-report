#!/usr/bin/env python
""" MultiQC functions to use plotly library """

import logging
import base64
from plotly.io import to_json
from multiqc.utils import report

logger = logging.getLogger(__name__)


def fig_to_json_html(fig, pconfig):
    if pconfig.get('auto_margin'):
        fig.update_layout(margin=dict(l=40, r=20, t=40, b=40))

    if pconfig.get('ylab'):
        fig.update_layout(yaxis=dict(title_text=pconfig['ylab']))

    if pconfig.get('xlab'):
        fig.update_layout(xaxis=dict(title_text=pconfig['xlab']))

    if pconfig.get('title'):
        fig.update_layout(title_text=pconfig['title'], title_x=0.5)

    json_str = to_json(fig)
    html = '<script id="{id}" type="text/json">{json}</script>'.format(
        id=pconfig['data_id'], json=json_str)
    return html


def plot(fig, pconfig):
    data_html = fig_to_json_html(fig, pconfig)
    html = '''
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <div class="hc-plot-wrapper">
    <div id="{id}" class="hc-plot not_rendered">
      <small>loading..</small>
    </div>
  </div>
  {data_html}
  <script type="text/javascript">
    var figure = JSON.parse($("#{data_id}").html());
    figure.layout.autosize = true
    Plotly.newPlot("{id}", figure.data, figure.layout);
    // When plotly is working, hide something
    $("#{id}").removeClass("not_rendered");
    $("#{id} small").hide();
    // update the layout to expand to the available size
    // when the window is resized
    window.onresize = function() {{
      $(".js-plotly-plot").each(function () {{
        const id = $(this).attr('id');
        Plotly.relayout(id, {{
          "xaxis.autorange": true,
          "yaxis.autorange": true
        }});
      }})
    }};
  </script>
  '''.format(id=pconfig['id'], data_id=pconfig['data_id'], data_html=data_html)
    return html