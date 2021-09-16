"""
OMBE Temperature Calibration
@author: Alexandra Tully
@date: July 2021
"""

import numpy as np
import plotly.graph_objects as go
from arpes_functions import fitting_functions as ff

DEFAULT_RENDERER = 'browser'

xdata = np.array([50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 275, 300, 310, 320, 330, 340, 350, 378, 400, 425, 450,
                  475, 500, 525, 550, 575, 598])
ydata = np.array([41.9, 60.5, 76.6, 110.4, 151.8, 171.3, 189.7, 217, 233, 256, 256, 280, 286, 295, 304, 311, 321, 344,
                  362, 382, 401, 420, 439, 427, 454, 474, 492])

fig = go.Figure()
line = go.Scatter(x=xdata, y=ydata, mode='markers+lines', name='data')
fig.add_trace(line)
fig.update_layout(title='Temperature Calibration (OMBE Heater)',
                  xaxis_title='Manipulator Temperature', yaxis_title='Sample Temperature')
fig.update_layout(template='plotly_dark')
fig.show(renderer=DEFAULT_RENDERER)

fit1 = ff.fit_linear_data(x=xdata, data=ydata,
                          num=1,
                          aes=1, bes=1,
                          offset_type=None)

fit2 = ff.offset_model(offset_type='quadratic')
fig.add_trace(go.Scatter(x=xdata,
                         y=fit2.eval(x=xdata),
                         mode='lines', name='fit', line=dict(color='turquoise', dash='dash')))

fig.add_trace(go.Scatter(x=xdata,
                         y=fit1.eval(x=xdata),
                         mode='lines', name='linear fit', line=dict(color='turquoise', dash='dash')))
fig.show()