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

from misc_functions import bin_2D, bin_2D_array, sanitize
from plotting_functions import filepath_plot2D, filepath_plot3D, plot2D, plot3D
from arpes_dataclasses import Data3D, Data2D
from HDF5_loader import data_from_hdf
from k_correction import fix_EkatEF, kcorrect2D, kcorrect3D, kcorrect_phimotor
from polygons import gen_polygon, gen_tiled_hexagons, plot_polygon, plot_polygons

"""Filepath Definition (Desktop and Macbook)"""
PATH = 'C:/Users/atully/Code/ARPES Code Python/analysis_data/'
# PATH = os.path.abspath('/Users/alexandratully/Desktop/ARPES Data/')


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
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-fn', placeholder='Type full filename (+ parent'
                                                                                        ' folder for phims)',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Year:", addon_type="prepend"),
         dcc.Dropdown(id='drop-year',
                      options=[
                          {'label': '2020', 'value': '2020'},
                          {'label': '2021', 'value': '2021'}
                      ],
                      value='2021',
                      multi=False,
                      style={'width': '70%'},
                      persistence=True,
                      persistence_type='local'
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Month:", addon_type="prepend"),
         dcc.Dropdown(id='drop-month',
                      value='January',
                      multi=False,
                      style={'width': '70%'},
                      persistence=True,
                      persistence_type='local'
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Light Source:", addon_type="prepend"),
         dcc.Dropdown(id='drop-laser',
                      options=[
                          {'label': 'He Lamp', 'value': 'lamp'},
                          {'label': '6eV', 'value': '6eV'},
                          {'label': 'XUV', 'value': 'XUV'}
                      ],
                      value='Lamp',
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
                          {'label': 'Low Temp', 'value': 'LT'},
                      ],
                      value='LT',
                      multi=False,
                      style={'width': '70%'},
                      persistence=True,
                      persistence_type='local'
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
                          {'label': 'Deflector Voltages', 'value': 'deflector_scan'},
                          {'label': 'Phi Motor', 'value': 'phi_motor_scan'},
                      ],
                      value='deflector_scan',
                      multi=False,
                      style={'width': '75%'}
                      )
         ], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Slice Dimension (3D):", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-slice-dim', placeholder='x, y, or z (dvs: z, '
                                                                                               'phims: y)',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Integration Range:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-int', placeholder='Choose Int (Float, '
                                                                                         'default=0.02)',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Scan Number:", addon_type="prepend"),
         dbc.Input(type='number', bs_size='sm', autoFocus=True, id='inp-snum', placeholder='Choose Scan Num',
                   min=0, debounce=True)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Fermi Energy:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-ef', placeholder='Choose E_f (Float, '
                                                                                        'default=16.8)',
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
side_bar_2 = html.Div([
    dbc.Checklist(
        options=[
            {'label': 'K-Correction', 'value': 'k'}
        ],
        value=['k'],
        id="tog-k",
        # switch=True
    ),
    'Angles:',
    dbc.InputGroup(
        [dbc.InputGroupAddon("Theta0:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-theta',
                   min=0, debounce=True, value=0)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Phi0:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-phi',
                   min=0, debounce=True, value=0)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Theta_m:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-thetam',
                   min=0, debounce=True, value=0)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Phi_m:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-phim',
                   min=0, debounce=True, value=0)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Alpha:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-alpha',
                   min=0, debounce=True, value=0)], size="sm"),
    '',
    'Hexagon Overlay:',
    dbc.InputGroup(
        [dbc.InputGroupAddon("Radius:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-hexrad', placeholder='Default: 0.42 A-1',
                   min=0, debounce=True, value=0.42)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Rotation:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-hexrot', placeholder='Degrees',
                   min=0, debounce=True, value=30)], size="sm"),
    dbc.InputGroup(
        [dbc.InputGroupAddon("Center:", addon_type="prepend"),
         dbc.Input(type='value', bs_size='sm', autoFocus=True, id='inp-hexcent', placeholder='Default: (0,0)',
                   min=0, debounce=True, value=(0, 0))], size="sm"),
    # 'Hexagon Center X:',
    # dcc.Slider(id='sl-hexcent-x', updatemode='drag'),
    # 'Hexagon Center Y:',
    # dcc.Slider(id='sl-hexcent-y', updatemode='drag'),
    dbc.RadioItems(
        options=[
            {'label': 'Plot Hexagon', 'value': 1},
            {'label': 'Plot Hexagons', 'value': 2},
            {'label': 'None', 'value': 3}
        ],
        value=3,
        id="rad-hex",
        inline=True,
    )
])

# Make overall layout
app_layout = html.Div(
    [
        dbc.Row(top_bar),
        dbc.Row([
            dbc.Col(main_area, width=6), dbc.Col(side_bar, width=4), dbc.Col(side_bar_2, width=2)
            # Dash app width = 12 total
        ])
    ], style={'width': '90vw', 'margin': 'auto'}
)

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
    Input('tog-k', 'value'),
    Input('inp-theta', 'value'),
    Input('inp-phi', 'value'),
    Input('inp-thetam', 'value'),
    Input('inp-phim', 'value'),
    Input('inp-alpha', 'value'),
    Input('rad-hex', 'value'),
    Input('inp-hexrad', 'value'),
    Input('inp-hexrot', 'value'),
    Input('inp-hexcent', 'value'),
    # Input('sl-hexcent-x', 'value'),
    # Input('sl-hexcent-y', 'value'),
    State('drop-month', 'value'),
    State('drop-year', 'value'),
    State('drop-laser', 'value'),
    State('drop-color', 'value'),
    State('drop-sdim', 'value'),
    State('inp-snum', 'value'),
    State('drop-stype', 'value'),
    State('inp-slice-dim', 'value'),
    State('inp-int', 'value'),
    State('drop-temp', 'value'),
    State('inp-ef', 'value'),
    State('inp-fn', 'value')
)
def get_figure(clicks, slice_val, slice_line, bin_x, bin_y, kcorrect, theta0, phi0, theta_m, phi_m, alpha, hexagons,
               hex_radius, hex_rotation, hex_center, month, year, laser,
               colorscale, scan_dim, scan_num, scan_type, slice_dim, int_range, cryo_temp, E_f, filename):
    # hex_center_slice_x, hex_center_slice_y,
    theta0 = sanitize(theta0, 'float', default_val=0)
    theta_m = sanitize(theta_m, 'float', default_val=0)
    phi0 = sanitize(phi0, 'float', default_val=0)
    phi_m = sanitize(phi_m, 'float', default_val=0)
    alpha = sanitize(alpha, 'float', default_val=0)
    hex_radius = sanitize(hex_radius, 'float', default_val=0)
    hex_rotation = sanitize(hex_rotation, 'float', default_val=0)
    E_f = sanitize(E_f, 'float', default_val=16.8)
    if type(hex_center) == str:
        hex_center = sanitize(hex_center, 'csv-to-valuelist', default_val=0)
    if not scan_num:
        scan_num = 1
    if all([year, month, colorscale, scan_dim, scan_num, scan_type, cryo_temp]):
        if scan_dim == '2D':
            try:
                data = Data2D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp,
                                          scan_number=scan_num, filename=filename, filepath=PATH)
                check_filepath = data.fp
                if kcorrect:
                    if not E_f:
                        E_f = 16.8
                    fix_EkatEF(data, E_f)
                    x, y, d = kcorrect2D(data)
                    xlabel, ylabel = 'kx [A-1]', 'Energy'
                elif not kcorrect:
                    x, y, d = data.xaxis, data.yaxis, data.data
                    xlabel, ylabel = 'Theta', 'Energy'
                if bin_x is None or bin_x == 0:
                    bin_x = 1
                if bin_y is None or bin_y == 0:
                    bin_y = 1
                if bin_x > 1 or bin_y > 1:
                    x, y, d = bin_2D_array(d, x, y, bin_x, bin_y)
                title = f'{laser} 2D: Scan #{scan_num} ({cryo_temp})'
                fig = plot2D(x, y, d, show=False, colorscale=colorscale, E_f=E_f, xlabel=xlabel,
                             ylabel=ylabel, title=title)
                if slice_line:
                    fig.update_layout(shapes=[dict(type='line', y0=slice_line, y1=slice_line, xref='paper', x0=0,
                                                   x1=1)])
            except FileNotFoundError as e:
                fig = go.Figure()
                check_filepath = str(e)
            return fig, check_filepath
        elif scan_dim == '3D':
            if int_range is None or int_range == '':
                int_range = 0.02
            try:
                slice_val = float(slice_val)
                int_range = float(int_range)
                if scan_type == 'deflector_scan':
                    if slice_dim is None or slice_dim == '':
                        slice_dim = 'z'
                    data = Data3D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp,
                                              scan_type=scan_type, scan_number=scan_num, filename=filename,
                                              filepath=PATH)
                    check_filepath = data.fp
                    title = f'Angular Scan (Deflector Voltages): Scan #{scan_num} (slice dimension: {slice_dim}, ' \
                            f'temp: {cryo_temp})'
                    if not E_f:
                        E_f = 16.8
                    if kcorrect:
                        # if not E_f:
                        #     E_f = 16.8
                        # fix_EkatEF(data, E_f)
                        kx, ky, kd = kcorrect3D(data, val=slice_val, Eint=int_range, fixEF=True, EF=E_f,
                                                slice_dim=slice_dim, theta0=theta0, phi0=phi0, theta_m=theta_m,
                                                phi_m=phi_m, alpha0=alpha)
                        fig = plot2D(kx, ky, kd, show=False, xlabel='kx [A-1]', ylabel='ky [A-1]', title=title)
                        if hexagons == 3:
                            print('No hexagons plotted')
                        elif hexagons == 2 or 1:
                            print(f'hexagons: {hexagons}, {hex_radius}')
                            # hex_cent = (hex_center_slice_x, hex_center_slice_y)
                            coords = gen_polygon(6, radius=hex_radius, rotation=hex_rotation, translation=hex_center)
                            plot_polygon(coords, fig=fig, show=False)
                            if hexagons == 2:
                                new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(
                                    radius=hex_radius, rotation=hex_rotation, translation=hex_center)
                                plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r],
                                              fig=fig, show=False)
                        else:
                            return ValueError(f'{hexagons} must be 1, 2, or 3')
                    elif not kcorrect:
                        x, y, z, d = data.xaxis, data.yaxis, data.zaxis, data.data
                        fig = plot3D(y, x, z, data=np.moveaxis(d, 2, 1),
                                     slice_dim=slice_dim, slice_val=slice_val + E_f, int_range=int_range, show=False,
                                     colorscale=colorscale, xlabel='Theta (Deg)', ylabel='Phi (Deg)', title=title)
                elif scan_type == 'phi_motor_scan':
                    if slice_dim is None or slice_dim == '':
                        slice_dim = 'y'
                    fp = os.path.join(PATH, f'{month}_{year}', cryo_temp, laser, '3D', scan_type)
                    if filename is None:
                        raise ValueError('Filename (including parent folder) is required for phi motor scans')
                    check_filepath = os.path.join(fp, filename)
                    title = f'Angular Scan (Phi Motor): Scan #{scan_num} (slice dimension: {slice_dim}, temp: ' \
                            f'{cryo_temp})'
                    if not E_f:
                        E_f = 16.8
                    if kcorrect:
                        print(theta0, phi0, theta_m, phi_m, alpha)
                        kx, ky, kdata = kcorrect_phimotor(fp=fp, fn=filename, slice_dim=slice_dim, val=slice_val,
                                                          Eint=int_range, EF=E_f, theta0=theta0, phi0=phi0,
                                                          theta_m=theta_m, phi_m=phi_m, alpha0=alpha)
                        fig = plot2D(kx, ky, kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=title, show=False)
                        if hexagons == 3:
                            print('No hexagons plotted')
                        elif hexagons == 2 or 1:
                            print(f'hexagons: {hexagons}, {hex_radius}')
                            # hex_cent = (hex_center_slice_x, hex_center_slice_y)
                            coords = gen_polygon(6, radius=hex_radius, rotation=hex_rotation, translation=hex_center)
                            plot_polygon(coords, fig=fig, show=False)
                            if hexagons == 2:
                                new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(
                                    radius=hex_radius, rotation=hex_rotation, translation=hex_center)
                                plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r],
                                              fig=fig, show=False)
                        else:
                            return ValueError(f'{hexagons} must be 1, 2, or 3')
                    elif not kcorrect:
                        d, x, y, z = data_from_hdf(fp, f'{filename}')
                        if slice_dim == 'y':
                            z = np.flip(z)
                        else:
                            print(f'Slice dimension is {slice_dim} not y, might need to flip ky axis')
                        fig = plot3D(x, y, z, data=d, slice_dim=slice_dim, slice_val=slice_val + E_f,
                                     int_range=int_range, show=False, colorscale=colorscale, xlabel='Theta (Deg)',
                                     ylabel='Phi (Deg)', title=title)
                if slice_line:
                    fig.update_layout(shapes=[dict(type='line', y0=slice_line, y1=slice_line, xref='paper', x0=0,
                                                   x1=1)])
            except (FileNotFoundError, OSError) as e:
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
    # print(type(year))
    if year == '2020':
        opts = [
            {'label': 'October', 'value': 'October', 'disabled': False},
            {'label': 'December', 'value': 'December', 'disabled': False},
            {'label': 'January', 'value': 'January', 'disabled': True},
            {'label': 'April', 'value': 'April', 'disabled': True}
        ]
        return opts
    elif year == '2021':
        opts = [
            {'label': 'October', 'value': 'October', 'disabled': True},
            {'label': 'December', 'value': 'December', 'disabled': True},
            {'label': 'January', 'value': 'January', 'disabled': False},
            {'label': 'April', 'value': 'April', 'disabled': False}
        ]
        return opts
    else:
        raise ValueError(f'No year {year}')


@app.callback(
    # Output('drop-subtype', 'disabled'),
    Output('sl-slice', 'disabled'),
    Output('inp-int', 'disabled'),
    Output('drop-stype', 'disabled'),
    Output('inp-slice-dim', 'disabled'),
    Input('drop-sdim', 'value')
)
def toggle_scantype_active(scan_dim):
    if scan_dim:
        if scan_dim == '2D':
            return [True, True, True, True]
        elif scan_dim == '3D':
            return [False, False, False, False]
        else:
            raise ValueError(f'{scan_dim} not 2D or 3D')


@app.callback(
    Output('inp-ef', 'disabled'),
    Input('drop-sdim', 'value')
)
def toggle_ef_active(scan_dim):
    if scan_dim:
        if scan_dim == '2D':
            return False
        else:
            print('E_f enabled for all')
            return False


# TODO: have type in box for constant slice 3D in binding energy (so can type in -1.1 and it gives that slice, based on
#  the Ef)

@app.callback(
    Output('sl-slice', 'min'),
    Output('sl-slice', 'max'),
    Output('sl-slice', 'marks'),
    Output('sl-slice', 'value'),
    Output('sl-slice', 'step'),
    Input('but-load', 'n_clicks'),
    Input('inp-ef', 'value'),
    State('drop-month', 'value'),
    State('drop-year', 'value'),
    State('drop-sdim', 'value'),
    State('inp-snum', 'value'),
    State('drop-stype', 'value'),
    # State('drop-subtype', 'value'),
    State('inp-slice-dim', 'value'),
    State('drop-temp', 'value'),
    State('drop-laser', 'value'),
    State('inp-fn', 'value')
)
def get_slice_values(clicks, E_f, month, year, scan_dim, scan_num, scan_type,
                     slice_dim, cryo_temp, laser, filename):
    E_f = sanitize(E_f, 'float', default_val=16.8)
    if not E_f:
        E_f = 16.8
    if not scan_num:
        scan_num = 1
    if scan_type == 'deflector_scan':
        if slice_dim is None or slice_dim == '':
            slice_dim = 'z'
    elif scan_type == 'phi_motor_scan':
        if slice_dim is None or slice_dim == '':
            slice_dim = 'y'
    slice_min, slice_max = (0, 1)
    value = 0
    if all([month, year, scan_dim, scan_num, scan_type, slice_dim, cryo_temp]):
        if scan_dim and scan_dim == '3D':
            try:
                if scan_type == 'deflector_scan':
                    data = Data3D.single_load(month, year=year, cryo_temp=cryo_temp, light_source=laser,
                                              scan_type=scan_type, scan_number=scan_num, filepath=PATH)
                    x, y, z = data.xaxis, data.yaxis, data.zaxis
                elif scan_type == 'phi_motor_scan':
                    fp = os.path.join(PATH, f'{month}_{year}', cryo_temp, laser, '3D', scan_type)
                    d, x, y, z = data_from_hdf(fp, f'{filename}')
                if slice_dim == 'x':
                    slice_min = x[0]
                    slice_max = x[-1]
                elif slice_dim == 'y':
                    slice_min = y[0]
                    slice_max = y[-1]
                elif slice_dim == 'z':
                    slice_min = z[0]
                    slice_max = z[-1]
                else:
                    raise ValueError(f'{slice_dim} not in x, y, z')
            except (FileNotFoundError, OSError, ValueError):
                fig = go.Figure()
    slice_min -= E_f
    slice_max -= E_f
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
    if fig is None:
        return 0, 1, {}, 1
    data = fig.get('data')
    if data:
        y_len = len(data[0]['y'])
        min = data[0]['y'][0]
        max = np.array(data[0]['y'])[-1]
        if min > max:
            print(f'min: {min}, max: {max}')
            m1, m2 = min, max
            min, max = m2, m1
            print(f'min: {min}, max: {max}')
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
