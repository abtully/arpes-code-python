from __future__ import annotations
from typing import List, Tuple, Dict, Optional, Callable
import logging
import pandas as pd

from dash_dashboard.base_classes import PageInteractiveComponents, \
    CommonInputCallbacks, PendingCallbacks
from ali_dash.base_class_overrides import AliMain, AliSidebar, AliPageLayout
import dash_dashboard.component_defaults as c
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash import no_update

from leed_analysis.leed_import import get_filepath, load_image, show_image

logger = logging.getLogger(__name__)

NAME = 'LEED'
URL_ID = 'LEED'
page_collection = None  # Gets set when running in multipage mode


class Components(PageInteractiveComponents):
    def __init__(self, pending_callbacks: Optional[PendingCallbacks] = None):
        super().__init__(pending_callbacks)

        # sidebar components
        self.txt_filepath = html.P(id='txt-filepath')
        self.inp_example = c.input_box(id_name='inp-exmaple', val_type='text', debounce=False,
                                       placeholder='Example Input', persistence=False)
        self.inp_filepath = c.input_box(id_name='inp-filepath', val_type='text', debounce=True,
                                        placeholder='PATH (excl. filename)', persistence=True)
        self.inp_filename = c.input_box(id_name='inp-filename', val_type='text', debounce=True,
                                        placeholder='filename (.tiff)', persistence=True)
        self.ref_table = c.table(id_name='tab-ref', editable=True)
        self.dat_table = c.table(id_name='tab-dat', editable=True)
        self.tog_ref = c.toggle(id_name='tog-ref', label='Reference Image')
        self.txt_printout = html.P(id='txt-printout')

        # main area components
        self.div_example = c.div(id_name='div-example')
        self.fig_leed = c.graph_area(id_name='fig-leed', graph_header='LEED Data')
        self.clicks = html.Pre(id='click-data')

        self.initial_setup()

    def initial_setup(self):
        self.ref_table.columns = [{'name': n, 'id': n} for n in ['x', 'y', 'Intensity', 'n', 'eV', 'a']]
        self.dat_table.columns = [{'name': n, 'id': n} for n in ['x', 'y', 'Intensity', 'n', 'eV', 'a']]


class LEEDLayout(AliPageLayout):
    top_bar_title = 'Title -- May want to override in PageLayout override'

    # Defining __init__ only for typing purposes (i.e. to specify page specific Components as type for self.components)
    def __init__(self, components: Components):
        super().__init__(page_components=components)
        self.components = components

    def get_mains(self) -> List[AliMain]:
        return [LEEDMain(self.components), ]

    def get_sidebar(self) -> AliSidebar:
        return LEEDSidebar(self.components)


class LEEDMain(AliMain):
    name = 'LEED'

    # Defining __init__ only for typing purposes (i.e. to specify page specific Components as type for self.components)
    def __init__(self, components: Components):
        super().__init__(page_components=components)
        self.components = components

    def layout(self):
        lyt = html.Div([
            self.components.div_example,
            self.components.fig_leed,
            self.components.clicks,
        ])
        return lyt

    def set_callbacks(self):
        self.make_callback(outputs=(self.components.div_example.id, 'children'), inputs=[
            (self.components.inp_example.id, 'value'),
        ], func=lambda text: text)

        self.make_callback(outputs=(self.components.fig_leed.graph_id, 'figure'),
                           inputs=LeedFigureCallbacks.get_inputs(),
                           states=LeedFigureCallbacks.get_states(),
                           func=LeedFigureCallbacks.get_callback_func('LEED Image'))

        self.make_callback(outputs=(self.components.clicks.id, 'children'),
                           inputs=(self.components.fig_leed.graph_id, 'clickData'),
                           # func=lambda clickData: json.dumps(clickData, indent=2))
                           func=get_click_data)


def get_click_data(clickData):
    if clickData:
        point = clickData['points'][0]
        return f'x: {point["x"]}, y: {point["y"]}'
    else:
        return 'no clicks'


class LEEDSidebar(AliSidebar):
    id_prefix = 'LEEDSidebar'

    # Defining __init__ only for typing purposes (i.e. to specify page specific Components as type for self.components)
    def __init__(self, components: Components):
        super().__init__(page_components=components)
        self.components = components

    def layout(self):
        lyt = html.Div([
            dbc.Toast(
                self.components.txt_filepath,
                header="Current Filepath", style={'maxWidth': '600px'}
            ),
            self.components.dd_main,
            self.components.inp_example,
            self.input_wrapper('Filepath:', self.components.inp_filepath),
            self.components.inp_filename,
            c.space(height='10px'),
            html.Hr(),
            self.components.tog_ref,
            html.H5('Reference Table'),
            self.components.ref_table,
            c.space(height='10px'),
            html.H5('Data Table'),
            self.components.dat_table,
            self.components.txt_printout
        ])
        return lyt

    def set_callbacks(self):
        self.make_callback(outputs=(self.components.txt_filepath.id, 'children'),
                           inputs=[(self.components.inp_filepath.id, 'value'),
                                   (self.components.inp_filename.id, 'value')],
                           func=lambda path, name: get_filepath(path, name))

        # self.make_callback(outputs=(self.components.ref_table.id, 'data'),
        #                    inputs=(self.components.fig_leed.graph_id, 'clickData'),
        #                    states=[(self.components.ref_table.id, 'data'),
        #                            (self.components.tog_ref.id, 'value')],
        #                    func=ref_table)
        #
        # self.make_callback(outputs=(self.components.dat_table.id, 'data'),
        #                    inputs=(self.components.fig_leed.graph_id, 'clickData'),
        #                    states=[(self.components.dat_table.id, 'data'),
        #                            (self.components.tog_ref.id, 'value')],
        #                    func=dat_table)
        #
        self.make_callback(outputs=(self.components.ref_table.id, 'data'),
                           inputs=TableCallbacks.get_inputs(),
                           states=TableCallbacks.get_states(),
                           func=TableCallbacks.get_callback_func('ref_table'))

        self.make_callback(outputs=(self.components.dat_table.id, 'data'),
                           inputs=TableCallbacks.get_inputs(),
                           states=TableCallbacks.get_states(),
                           func=TableCallbacks.get_callback_func('dat_table'))

        self.make_callback(outputs=(self.components.txt_printout.id, 'children'),
                           inputs=TableCallbacks.get_inputs(),
                           states=TableCallbacks.get_states(),
                           func=TableCallbacks.get_callback_func('table_printout'))


# def ref_table(clickData, data, value):
#     if value == [1]:
#         return click_to_table(clickData, data)
#     else:
#         raise ValueError(f'Value {value} is not 1')
#
#
# def dat_table(clickData, data, value):
#     if not value:
#         return click_to_table(clickData, data)
#     else:
#         raise ValueError(f'Value {value} is not None')
#
#
# def click_to_table(clickData, data):
#     if clickData:
#         x, y = clickData['points'][0]['x'], clickData['points'][0]['y']
#         df_coords = []
#         I, n, eV, a = 0, 1, 0, 0
#         if data:
#             df_coords = [[d['x'], d['y'], d['Intensity'], d['n'], d['eV'], d['a']] for d in data][-1:]
#             I, n, eV, a = data[-1]['Intensity'], data[-1]['n'], data[-1]['eV'], data[-1]['a']
#         df = pd.DataFrame([*df_coords, [x, y, I, n, eV, a]],
#                           columns=['x', 'y', 'Intensity', 'n', 'eV', 'a'])  # * means unpack outer list
#         return df.to_dict('records')


class TableCallbacks(CommonInputCallbacks):
    components = Components()

    def __init__(self, clickData, ref_data, dat_data, tog_value):
        super().__init__()
        self.clickData = clickData
        self.ref_data = ref_data
        self.dat_data = dat_data
        self.tog_value = True if tog_value == [1] else False

    @classmethod
    def get_inputs(cls) -> List[Tuple[str, str]]:
        return [
            (cls.components.fig_leed.graph_id, 'clickData')
        ]

    @classmethod
    def get_states(cls) -> List[Tuple[str, str]]:
        return [(cls.components.ref_table.id, 'data'),
                (cls.components.dat_table.id, 'data'),
                (cls.components.tog_ref.id, 'value')]

    def callback_names_funcs(self) -> Dict[str, Callable]:
        return {
            'ref_table': self.ref_table,
            'dat_table': self.dat_table,
            'get_distance': self.get_distance,
            'table_printout': self.table_printout
        }

    def ref_table(self):
        if self.tog_value is True:
            return self.click_to_table(self.ref_data)
        else:
            return no_update

    def dat_table(self):
        if self.tog_value is False:
            return self.click_to_table(self.dat_data)
        else:
            return no_update

    def click_to_table(self, data):
        if self.clickData:
            x, y = self.clickData['points'][0]['x'], self.clickData['points'][0]['y']
            df_coords = []
            I, n, eV, a = 0, 1, 0, 0
            if data:
                df_coords = [[d['x'], d['y'], d['Intensity'], d['n'], d['eV'], d['a']] for d in data][-1:]
                I, n, eV, a = data[-1]['Intensity'], data[-1]['n'], data[-1]['eV'], data[-1]['a']
            df = pd.DataFrame([*df_coords, [x, y, I, n, eV, a]],
                              columns=['x', 'y', 'Intensity', 'n', 'eV', 'a'])  # * means unpack outer list
            return df.to_dict('records')

    def get_distance(self, data):
        if data is not None:

            x = [d['x'] for d in data][-1:][0]
            y = [d['y'] for d in data][-1:][0]
            return [d['x'] for d in data]

    def table_printout(self, variable=None):
        data = self.click_to_table(self.ref_data)
        var = self.get_distance(data)
        return f'{var}'

    # def table_printout(self):
    #     return f'{self.ref_data}'


class LeedFigureCallbacks(CommonInputCallbacks):  # useful for figure that uses a bunch of the same inputs
    components = Components()  # gets ids from above

    def __init__(self, path, name):  # order matters!!
        super(LeedFigureCallbacks, self).__init__()
        self.path = path
        self.name = name

    @classmethod
    def get_inputs(cls) -> List[Tuple[str, str]]:
        return [
            (cls.components.inp_filepath.id, 'value'),  # order matters!!
            (cls.components.inp_filename.id, 'value')
        ]

    @classmethod
    def get_states(cls) -> List[Tuple[str, str]]:  # order matters!!
        pass

    def callback_names_funcs(self) -> Dict[str, Callable]:
        return {
            'LEED Image': self.get_leed_figure
        }

    def get_leed_figure(self):
        # full_fp = get_filepath(self.path, self.name)
        data = load_image(self.name, self.path)
        fig = show_image(data, show=False)
        return fig


"""Multipage layout requirements"""


def layout(*args):  # *args only because dash_extensions passes in the page name for some reason
    inst = LEEDLayout(Components())
    inst.page_collection = page_collection
    return inst.layout()


def callbacks(app):
    inst = LEEDLayout(Components(pending_callbacks=PendingCallbacks()))
    inst.page_collection = page_collection
    inst.layout()  # Most callbacks are generated while running layout
    return inst.run_all_callbacks(app)


if __name__ == '__main__':
    from dash_dashboard.app import test_page

    test_page(layout=layout, callbacks=callbacks, port=8062)
