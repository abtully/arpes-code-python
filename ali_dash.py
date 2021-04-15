import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, MATCH, ALL, ALLSMALLER
from typing import List, Union, Optional, Tuple
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
import os

from ali_functions import bin_2D
from ali_plotting_functions import filepath_plot2D, filepath_plot3D, plot2D, plot3D
from ali_classes import Data3D, Data2D

"""Filepath Definition (Desktop and Macbook)"""
Path = 'C:/Users/atully/Code/ARPES Code Python/data/'
# Path = os.path.abspath('/Users/alexandratully/Desktop/ARPES Code Python/data/')


"""Dash App Begins"""
# Setup some theme stuff
THEME = ''
if THEME == 'dark':
    DASH_THEME = dbc.themes.DARKLY
    plotly_theme = 'plotly_dark'
else:
    DASH_THEME = dbc.themes.BOOTSTRAP
    plotly_theme = 'none'
pio.templates.default = plotly_theme

# Individual layout components
top_bar = dbc.Col(html.Div(html.H1('OMBE ARPES Data')))
main_area = html.Div([dbc.Card(dcc.Graph(id='fig')), dbc.Card(dcc.Graph(id='fig2'))])
side_bar = html.Div([
    dbc.Toast(
        [html.P(id='fp')],
        header="Current Filepath", style={'maxWidth': '600px'}
    ),
    dbc.Button("Load", color="info", id='but-load'),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Filename Override:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-fn', placeholder='Type full filename',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Year:", addon_type="prepend"),
         dcc.Dropdown(id='drop-year',
                      options=[
                          {'label': '2020', 'value': '2020'},
                          {'label': '2021', 'value': '2021'}
                      ],
                      value='2020',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Month:", addon_type="prepend"),
         dcc.Dropdown(id='drop-month',
                      value='October',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Light Source:", addon_type="prepend"),
         dcc.Dropdown(id='drop-laser',
                      value='lamp',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Colorscale:", addon_type="prepend"),
         dcc.Dropdown(id='drop-color',
                      options=[
                          {'label': 'Plasma', 'value': 'Plasma'},
                          {'label': 'Teal', 'value': 'Teal'},
                          {'label': 'Gray', 'value': 'gray'},
                          {'label': 'Ice', 'value': 'ice'},
                          {'label': 'Greys', 'value': 'Greys'},
                      ],
                      value='Plasma',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Temperature:", addon_type="prepend"),
         dcc.Dropdown(id='drop-temp',
                      options=[
                          {'label': 'Room Temp', 'value': 'RT'},
                          {'label': '140K', 'value': '140K'},
                          {'label': '9K', 'value': '9K'}
                      ],
                      value='RT',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Scan Dimension:", addon_type="prepend"),
         dcc.Dropdown(id='drop-sdim',
                      options=[
                          {'label': '2D', 'value': '2D'},
                          {'label': '3D', 'value': '3D'}
                      ],
                      value='2D',
                      multi=False,
                      style={'width': '70%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Scan Type:", addon_type="prepend"),
         dcc.Dropdown(id='drop-stype',
                      options=[
                          {'label': 'High Symmetry', 'value': 'HS'},
                          {'label': 'Raster', 'value': 'Raster'},
                          {'label': 'UPS', 'value': 'UPS'},
                          {'label': 'Fermi Surface', 'value': 'FS'},
                      ],  # NOTE: this whole options chunk does nothing, it's called in toggle_subtype_active callback
                      value='HS',
                      multi=False,
                      style={'width': '75%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("K-Cut:", addon_type="prepend"),
         dcc.Dropdown(id='drop-kcut',
                      options=[
                          {'label': 'Gamma-M', 'value': 'GM'},
                          {'label': 'K-Gamma', 'value': 'KG'},
                      ],
                      value='GM',
                      multi=False,
                      style={'width': '75%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Scan Subtype:", addon_type="prepend"),
         dcc.Dropdown(id='drop-subtype',
                      options=[
                          {'label': 'Fermi Surface', 'value': 'Fermi_Surface'},
                          {'label': 'Fermi Edge', 'value': 'Fermi_Edge'},
                          {'label': 'Calibrate', 'value': 'calibrate'},
                          {'label': 'Coarse', 'value': 'coarse'},
                          {'label': 'Offset', 'value': 'offset'},
                          {'label': 'HOMO', 'value': 'HOMO'},
                          {'label': 'HOMO-1', 'value': 'HOMO-1'}
                      ],
                      value='Fermi_Surface',
                      multi=False,
                      style={'width': '75%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Slice Dimension (3D):", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-slice-dim', placeholder='x, y, or z (typ. z)',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Integration Range:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-int', placeholder='Choose Int (Float)',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Scan Number:", addon_type="prepend"),
         dbc.Input(type='number', bs_size='sm', autoFocus=True, id='inp-snum', placeholder='Choose Scan Num',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Fermi Energy:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-ef', placeholder='Choose E_f (Float)',
                   min=0, debounce=True)], size="sm"),
    'Constant Slice (3D):',
    dcc.Slider(id='sl-slice', updatemode='drag'),
    'Constant Slice (2D):',
    dcc.Slider(id='sl-slice2', updatemode='drag'),
    'Bin X Data:',
    dcc.Slider(id='sl-binx', min=0, max=10, marks={int(s): f'{s:.0f}' for s in np.linspace(0, 10, 11)}, value=5,
               updatemode='mouseup'),
    'Bin Y Data:',
    dcc.Slider(id='sl-biny', min=0, max=10, marks={int(s): f'{s:.0f}' for s in np.linspace(0, 10, 11)}, value=5,
               updatemode='mouseup'),
    dbc.Checklist(
        options=[
            {"label": "Symmetric Binning", "value": 1},
        ],
        value=[1],
        id="tog-bin",
        switch=True
    )
])

# Make overall layout
app_layout = html.Div(
    [
        dbc.Row(top_bar),
        dbc.Row([
            dbc.Col(main_area, width=8), dbc.Col(side_bar, width=4)  # Dash app width = 12 total
        ])
    ], style={'width': '90vw', 'margin': 'auto'}
)

# dbc.Container()
# style={'width': '90vw', 'position': 'absolute'}
# , fluid=False

# Make the dash app
app = dash.Dash(__name__, external_stylesheets=[DASH_THEME])
app.layout = app_layout


# Make callbacks
@app.callback(
    Output('fig', 'figure'),
    Output('fp', 'children'),
    Input('but-load', 'n_clicks'),
    Input('sl-slice', 'value'),
    Input('sl-slice2', 'value'),
    Input('sl-binx', 'value'),
    Input('sl-biny', 'value'),
    State('drop-month', 'value'),
    State('drop-year', 'value'),
    State('drop-laser', 'value'),
    State('drop-color', 'value'),
    State('drop-sdim', 'value'),
    State('inp-snum', 'value'),
    State('drop-stype', 'value'),
    State('drop-subtype', 'value'),
    State('drop-kcut', 'value'),
    State('inp-slice-dim', 'value'),
    State('inp-int', 'value'),
    State('drop-temp', 'value'),
    State('inp-ef', 'value'),
    State('inp-fn', 'value')
)
def get_figure(clicks, slice_val, slice_line, bin_x, bin_y, month, year, laser, colorscale, scan_dim, scan_num, scan_type,
               scan_subtype, k_cut, slice_dim, int_range, cryo_temp, E_f, filename):
    if not scan_num:
        scan_num = 1
    if all([year, month, colorscale, scan_dim, scan_num, scan_type, cryo_temp]):
        if scan_dim == '2D':
            try:
                data = Data2D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp, scan_type=scan_type,
                                          k_cut=k_cut, scan_number=scan_num, filename=filename, filepath=Path)
                x, y, d = data.xaxis, data.yaxis, data.data
                check_filepath = data.fp
                if bin_x is None or bin_x == 0:
                    bin_x = 1
                if bin_y is None or bin_y == 0:
                    bin_y = 1
                if bin_x > 1 or bin_y > 1:
                    x, y, d = bin_2D(data, bin_x, bin_y)
                if scan_type == 'HS':
                    title = f'High Symmetry Cut: {k_cut} ({cryo_temp})'
                else:
                    title = f'{scan_type}: Scan #{scan_num} ({cryo_temp})'
                fig = plot2D(x, y, d, show=False, colorscale=colorscale, E_f=E_f, xlabel='Theta', ylabel='Energy',
                             title=title)
                if slice_line:
                    fig.update_layout(shapes=[dict(type='line', y0=slice_line, y1=slice_line, xref='paper', x0=0,
                                                   x1=1)])
            except FileNotFoundError as e:
                fig = go.Figure()
                check_filepath = str(e)
            return fig, check_filepath
        elif scan_dim == '3D':
            if slice_dim is None or slice_dim == '':
                slice_dim = 'z'
            if int_range is None or int_range == '':
                int_range = 0.02
            try:
                # print(month, year, scan_type, scan_subtype, scan_num, slice_dim, slice_val, int_range)
                slice_val = float(slice_val)
                int_range = float(int_range)
                data = Data3D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp, scan_type=scan_type,
                                          scan_subtype=scan_subtype, scan_number=scan_num, filename=filename,
                                          filepath=Path)
                x, y, z, d = data.xaxis, data.yaxis, data.zaxis, data.data
                check_filepath = data.fp
                title = f'Fermi Surface: {scan_subtype}, Scan #{scan_num} (slice dimension: {slice_dim}, temp: ' \
                        f'{cryo_temp})'
                fig = plot3D(y, x, z, data=np.moveaxis(d, 2, 1),
                             slice_dim=slice_dim, slice_val=slice_val, int_range=int_range, show=False,
                             colorscale=colorscale, xlabel='Theta (Deg)', ylabel='Phi (Deg)', title=title)
                if slice_line:
                    fig.update_layout(shapes=[dict(type='line', y0=slice_line, y1=slice_line, xref='paper', x0=0,
                                                   x1=1)])
            except FileNotFoundError as e:
                fig = go.Figure()
                check_filepath = str(e)
            return fig, check_filepath
    else:
        return go.Figure()


@app.callback(
    Output('drop-month', 'options'),
    Input('drop-year', 'value')
)
def year_specifies_month(year):
    print(type(year))
    if year == '2020':
        opts = [
            {'label': 'October', 'value': 'October', 'disabled': False},
            {'label': 'December', 'value': 'December', 'disabled': False}
        ]
        return opts
    elif year == '2021':
        opts = [
            {'label': 'October', 'value': 'October', 'disabled': True},
            {'label': 'December', 'value': 'December', 'disabled': True}
        ]
        return opts
    else:
        raise ValueError(f'No year {year}')


@app.callback(
        Output('drop-laser', 'options'),
        Input('drop-year', 'value'),
        Input('drop-month', 'value')
)
def mmyy_specifies_laser(year, month):
    print(type(year))
    if year == '2020' and month == 'October':
        opts = [
            {'label': 'He Lamp', 'value': 'lamp', 'disabled': False},
            {'label': '6eV', 'value': '6eV', 'disabled': True},
            {'label': 'XUV', 'value': 'XUV', 'disabled': True}
        ]
        return opts
    elif year == '2020' and month == 'December':
        opts = [
            {'label': 'He lamp', 'value': 'lamp', 'disabled': False},
            {'label': '6eV', 'value': '6eV', 'disabled': True},
            {'label': 'XUV', 'value': 'XUV', 'disabled': False}
        ]
        return opts
    else:
        raise ValueError(f'No year {year}, month {month}')


@app.callback(
    Output('drop-kcut', 'disabled'),
    Input('drop-stype', 'value'),
    Input('drop-sdim', 'value')
)
def toggle_kcut_active(scan_type, scan_dim):
    if scan_type and scan_dim:
        if scan_type == 'Raster':
            return True
        if scan_type == 'UPS':
            return True
        if scan_dim == '3D':
            return True
    return False


@app.callback(
    Output('drop-subtype', 'disabled'),
    Output('sl-slice', 'disabled'),
    Output('inp-int', 'disabled'),
    Output('drop-stype', 'options'),
    Output('inp-slice-dim', 'disabled'),
    Input('drop-sdim', 'value')
)
def toggle_subtype_active(scan_dim):
    if scan_dim:
        if scan_dim == '2D':
            opts = [
                {'label': 'High Symmetry', 'value': 'HS'},
                {'label': 'Raster', 'value': 'Raster'},
                {'label': 'UPS', 'value': 'UPS'},
                {'label': 'Fermi Surface', 'value': 'FS', 'disabled': True},
            ]
            return [True, True, True, opts, True]
        elif scan_dim == '3D':
            opts = [
                {'label': 'High Symmetry', 'value': 'HS', 'disabled': True},
                {'label': 'Raster', 'value': 'Raster', 'disabled': True},
                {'label': 'UPS', 'value': 'UPS', 'disabled': True},
                {'label': 'Fermi Surface', 'value': 'FS'},
            ]
            return [False, False, False, opts, False]
        else:
            raise ValueError(f'{scan_dim} not 2D or 3D')


# # @lru_cache(maxsize=100)
# @app.callback(
#     Output('inp-snum', 'options'),
#     Input('but-load', 'n_clicks'),
#     State('drop-month', 'value'),
#     State('drop-year', 'value'),
#     State('drop-sdim', 'value'),
#     State('drop-stype', 'value'),
#     State('drop-subtype', 'value'),
#     State('drop-kcut', 'value'),
#     State('drop-temp', 'value')
# )
# def get_scan_num_opts(clicks, month, year, scan_dim, scan_type, scan_subtype, k_cut, cryo_temp):
#     ddir = os.path.normpath(Path)
#     if scan_dim == '2D':
#         fp = os.path.join(ddir, f'{month}_{year}', cryo_temp, '2D', scan_type, scan_subtype)
#     if scan_dim == '3D':
#         fp = os.path.join(ddir, f'{month}_{year}', cryo_temp, '3D')






@app.callback(
    Output('inp-ef', 'disabled'),
    Input('drop-sdim', 'value')
)
def toggle_ef_active(scan_dim):
    if scan_dim:
        if scan_dim == '2D':
            return False
        else:
            return True


@app.callback(
    Output('sl-slice', 'min'),
    Output('sl-slice', 'max'),
    Output('sl-slice', 'marks'),
    Output('sl-slice', 'value'),
    Output('sl-slice', 'step'),
    Input('but-load', 'n_clicks'),
    State('drop-month', 'value'),
    State('drop-year', 'value'),
    State('drop-sdim', 'value'),
    State('inp-snum', 'value'),
    State('drop-stype', 'value'),
    State('drop-subtype', 'value'),
    State('inp-slice-dim', 'value'),
    State('drop-temp', 'value')
)
def get_slice_values(clicks, month, year, scan_dim, scan_num, scan_type, scan_subtype,
                     slice_dim, cryo_temp):
    if not scan_num:
        scan_num = 1
    if slice_dim is None or slice_dim == '':
        slice_dim = 'z'
    slice_min, slice_max = (0, 1)
    value = 0
    if all([month, year, scan_dim, scan_num, scan_type, scan_subtype, slice_dim, cryo_temp]):
        if scan_dim and scan_dim == '3D':
            try:
                data = Data3D.single_load(month, year=year, cryo_temp=cryo_temp, scan_type=scan_type,
                                          scan_subtype=scan_subtype, scan_number=scan_num, filepath=Path)
                if slice_dim == 'x':
                    slice_min = data.xaxis[0]
                    slice_max = data.xaxis[-1]
                elif slice_dim == 'y':
                    slice_min = data.yaxis[0]
                    slice_max = data.yaxis[-1]
                elif slice_dim == 'z':
                    slice_min = data.zaxis[0]
                    slice_max = data.zaxis[-1]
                else:
                    raise ValueError(f'{slice_dim} not in x, y, z')
            except FileNotFoundError:
                fig = go.Figure()
    step = (slice_max - slice_min) / 15
    marks = {s: f'{s:.2f}' for s in np.linspace(slice_min, slice_max, 8)}
    value = np.mean([slice_min, slice_max])
    return slice_min, slice_max, marks, value, step


@app.callback(
    Output('sl-slice2', 'min'),
    Output('sl-slice2', 'max'),
    Output('sl-slice2', 'marks'),
    Output('sl-slice2', 'step'),
    Input('fig', 'figure')
)
def get_slice2_values(fig):
    data = fig.get('data')
    if data:
        y_len = len(data[0]['y'])
        min = data[0]['y'][0]
        max = np.array(data[0]['y'])[-1]
        # min, max= 0, y_len-1
        # step = 1
        step = abs(max - min) / 20
        marks = {i: f'{i:.2f}' for i in np.arange(min, max, abs(max - min) / 10)}
        marks = {k if not np.isclose(int(k), k) else int(k): v if not np.isclose(int(k), k) else str(int(k)) for k, v in
                 marks.items()}
        # mark = {int(i): f'{i:.0f}' for i in np.arange(0, 15, 1).astype(int)}
        # print(min, max, marks, step)
        return min, max, marks, step
    else:
        return 0, 1, {}, 1


@app.callback(
    Output('fig2', 'figure'),
    Input('sl-slice2', 'value'),
    Input('fig', 'figure')
)
def update_fig2(slice_line, fig1):
    if slice_line and fig1.get('data'):
        data = fig1.get('data')
        data = data[0]
        y_index = np.argmin(np.abs(np.array(data['y']) - slice_line)).astype(int)
        # print(y_index)
        one_d = data.get('z')[y_index]  # z = intensity in heatmaps
        x = data.get('x')
        fig = go.Figure(go.Scatter(mode='lines', x=x, y=one_d))
        fig.update_layout(title=f'1D Slice: {slice_line:.2f}')
        return fig
    return go.Figure()


@app.callback(
    Output('sl-biny', 'value'),
    Output('sl-biny', 'disabled'),
    Input('tog-bin', 'value'),
    Input('sl-binx', 'value'),
    State('sl-biny', 'value')
)
def sym_binning(bin_toggle, bin_x, bin_y):
    if bin_toggle == [1]:
        return bin_x, True
    else:
        return bin_y, False


if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload=True, port=8001)
