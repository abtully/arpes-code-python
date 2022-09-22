"""
LEED analysis functions
@author: Alexandra Tully
@date: April 2021
"""

from PIL import Image
import numpy as np
from scipy.constants import h, m_e
from typing import List, Optional

DEFAULT_RENDERER = 'browser'  # this is a constant of this file

PATH = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\October_2020\LEED\Spinview'


def intensity(array):
    def pixel_intensity(pixel):  # define pixel intensity subfunction
        return np.sqrt(pixel[0] ** 2 + pixel[1] ** 2 + pixel[2] ** 2)  # adding RBG values in quadrature

    array = array.astype(np.float32)  # change dtype from 8 byte integer to 32 byte float
    if array.ndim == 1 and len(array) == 4:  # if array is a pixel, utilize pixel_intensity subfunction
        return pixel_intensity(array)
    elif array.shape[-1] == 4:  # [-1] is the last value
        sq_array = array[:, :, 0:3] ** 2  # all x, all y, and 0-3 not inclusive, squared
        sum_sq = sq_array[:, :, 0] + sq_array[:, :, 1] + sq_array[:, :, 2]
        return np.sqrt(sum_sq)


def get_z_data(array, x, y):
    return array[:][y][x]


def get_line_distance(coord1: tuple, coord2: tuple) -> float:
    x1, y1 = coord1[0], coord1[1]
    x2, y2 = coord2[0], coord2[1]
    x = np.abs(x1 - x2)
    y = np.abs(y1 - y2)
    return np.sqrt(x ** 2 + y ** 2)


def atomic_row_spacing(x: float, R: float, E: int, n=1):  # E in eV, R and x in same units
    return (n * R * 12.3) / (x * np.sqrt(E))  # angstroms


def distance_to_screen(a: float, x: float, E: int, n=1):  # E in eV, a and x n same units
    return (a * x * np.sqrt(E)) / (12.3 * n)  # in angstroms


if __name__ == '__main__':
    from leed_import import load_image, show_image
    import os
    import plotly.graph_objects as go

    """Use LEED Analysis Dash Page to generate (x,y) coordinates for spots. Then proceed as follows:
    1. d_material = get_line_distance((x1, y1), (x2, y2))  # specify order of spots: n=1, and beam energy: 17eV
    2. d_reference = get_line_distance((x1, y1), (x2, y2))  # specify order of spots: n=1, and beam energy: 90eV
    2. r = distance_to_screen(a=2.88, x=d_reference, E=90, n=1)  # you know a from literature, x is calculated from 
                                                                image, E=beam energy, n=order of spots
    3. finally: atomic_row_spacing(x=d_material, R=r, E=17, n=1)  # you know x from image, r from reference calculation 
                                                                (distance to screen), E=beam energy, n=order of spots"""

    ### C60/Au(111) Growth Paper v2 (thorough) ###

    """
    Take au, c60 bilayer, and c60 5-6 ML LEED images from August 2022 (all same height). Get LEED spot coords
    from all images, in order. I know that differences in distances between coords is much greater source of error
    than small uncertainty in pixel location of spot (+/- 5 pixels) -- see robustness check below. Find all distances,
    calculate mean of each. Use au mean with a=2.88 as reference. Calculate a for c60 films. 
    a_c60_bilayer = 9.11 A
    a_c60_5/6ML = 10.26 A
    """

    au_spots = [(776, 518), (646, 682), (430, 651), (347, 447), (488, 276), (707, 317)]
    c60_thin_spots = [(745, 638), (522, 709), (338, 551), (392, 308), (627, 240), (799, 407)]
    c60_thick_spots = [(723, 620), (526, 688), (365, 547), (405, 331), (613, 268), (770, 415)]

    def get_line_distances(coords):
        prev_coord = coords[-1]
        distances = []
        for coord in coords:
            distances.append(get_line_distance(coord, prev_coord))
            prev_coord = coord
        return distances

    d_au = np.mean(get_line_distances(au_spots))
    d_c60_bilayer = np.mean(get_line_distances(c60_thin_spots))
    d_c60_thick = np.mean(get_line_distances(c60_thick_spots))

    r = distance_to_screen(a=2.88, x=d_au, E=123, n=1)
    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.11
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.26


    ### C60/Au(111) Growth Paper ###
    ddir = r'C:\Users\atully\OneDrive\Physics.UBC\Papers\C60 Growth\Figs\LEED'
    fn = r'au_for_paper.tiff'
    fp = os.path.join(ddir, fn)

    d_au = get_line_distance((707, 317), (775, 517))  # n=1, 123eV
    d_c60_bilayer = get_line_distance((799, 407), (745, 637))  # n=1, 10eV
    d_c60_thick = get_line_distance((773, 417), (723, 621))  # n=1, 10eV

    r = distance_to_screen(a=2.88, x=d_au, E=123, n=1)
    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.03
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.16

    "different spots"
    d_au = get_line_distance((430, 652), (348, 447))  # n=1, 123eV
    d_c60_bilayer = get_line_distance((522, 710), (745, 637))  # n=1, 10eV
    d_c60_thick = get_line_distance((527, 688), (720, 620))  # n=1, 10eV

    r = distance_to_screen(a=2.88, x=d_au, E=123, n=1)
    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.5
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.9

    "robustness check"
    b = 3
    c = 2
    d = -3
    g = -2
    d_c60_bilayer = get_line_distance((799 + b, 407 + c), (745 + d, 637 + g))  # n=1, 10eV
    d_c60_thick = get_line_distance((773 + b, 417 + c), (723 + d, 621 + g))  # n=1, 10eV

    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.12
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.27

    b = -5
    c = 5
    d = 4
    g = -4
    d_c60_bilayer = get_line_distance((799 + b, 407 + c), (745 + d, 637 + g))  # n=1, 10eV
    d_c60_thick = get_line_distance((773 + b, 417 + c), (723 + d, 621 + g))  # n=1, 10eV

    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.46
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.7

    d_au = get_line_distance((707 + 5, 317 - 5), (775 + 3, 517 - 4))  # n=1, 123eV
    d_c60_bilayer = get_line_distance((799, 407), (745, 637))  # n=1, 10eV
    d_c60_thick = get_line_distance((773, 417), (723, 621))  # n=1, 10eV

    r = distance_to_screen(a=2.88, x=d_au, E=123, n=1)
    atomic_row_spacing(x=d_c60_bilayer, R=r, E=10, n=1)  # 9.045
    atomic_row_spacing(x=d_c60_thick, R=r, E=10, n=1)  # 10.17

    # ### C60/Au(111) June 2020 ###
    # ddir = r'E:\atully\LEED data\2020_June'
    # fn = r'30_06_20_Au111_C60_1158.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((696, 451), (629, 333))  # n=1, 25eV
    # d_c60_n3 = get_line_distance((868, 431), (730, 695))  # n=3, 25eV
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=25, n=1)
    # atomic_row_spacing(x=d_c60_n3, R=r, E=25, n=3)  # n=3 -> 13.8; n=2 -> 9.2
    #
    # ### C60/Au(111) October 2020 ###
    # ddir = r'C:\\Users\\atully\\Code\\ARPES Code Python\\analysis_data\\October_2020\\LEED\\Spinview'
    # fn = r'14_10_20_C60_Au111_1654.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((699, 443), (612, 300))  # n=1, 17eV
    # d_c60_n2 = get_line_distance((813, 590), (559, 762))  # n=2, 17eV
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=17, n=1)
    # atomic_row_spacing(x=d_c60_n2, R=r, E=17, n=2)  # 11.0
    #
    #
    # ### C60/Au(111) April 2021 ###
    # ddir = r'E:\atully\LEED data\2021_April\April26'
    # fn = r'enhanced_04_26_21_Au111_C60_1130am_15eV.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((621, 362), (454, 414))  # n=1, 15eV
    # d_c60_n2 = get_line_distance((822, 281), (478, 212))  # n=2, 15eV
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=15, n=1)
    # atomic_row_spacing(x=d_c60_n2, R=r, E=15, n=2)  # 10.05
    #
    # ddir = r'E:\atully\LEED data\2021_April\April26'
    # fn = r'enhanced_04_26_21_Au111_C60_1130am_25eV.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((607, 381), (484, 423))  # n=1, 25eV
    # d_c60_n3 = get_line_distance((870, 415), (813, 715))  # n=3, 25eV
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=25, n=1)
    # atomic_row_spacing(x=d_c60_n3, R=r, E=25, n=2)  # n=3 -> 12.9, n=2 -> 8.6
    #
    # # ARPES LEED Image
    #
    # ddir = r'E:\atully\analysis_data\April_2021\LEED_ARPES\best spot'
    # fn = r'C60_30eV.bmp'
    # fp = os.path.join(ddir, fn)
    # #
    # # d_c60_n1 = get_line_distance((, ), (, ))  # n=1, 30eV
    # # d_c60_n3 = get_line_distance((, ), (, ))  # n=3, 30eV
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=30, n=1)
    # atomic_row_spacing(x=d_c60_n3, R=r, E=30, n=3)  # n=3 -> 12.9, n=2 -> 8.6
    #
    # ### C60/Graphene/Cu/SiC ###
    #
    # # C60/Graphene/Cu/SiC - 15eV
    # ddir = r'E:\atully\LEED data\2022_February\28_Feb_2022'
    # fn = r'enhanced_02_28_22_C60_graphene_SiC_1100am_15eV.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((486, 628), (433, 474))  # n=1, 15eV
    #
    # # C60/Graphene/Cu/SiC - 24eV
    # ddir = r'E:\atully\LEED data\2022_February\28_Feb_2022'
    # fn = r'enhanced_02_28_22_C60_graphene_SiC_1100am_24eV.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # d_c60_n1 = get_line_distance((509, 587), (461, 471))  # n=1, 24eV -> reference
    # d_c60_unknown = get_line_distance((815, 408), (640, 249))  # n=unknown, 24eV -> mystery spots
    # d_c60_n3 = get_line_distance((879, 545), (782, 268))  # n=3, 24eV -> double check (second order spots)
    #
    # r = distance_to_screen(a=10.08, x=d_c60_n1, E=24, n=1)
    # r2 = distance_to_screen(a=10.08, x=d_c60_unknown, E=24, n=2)
    #
    # atomic_row_spacing(x=d_c60_n3, R=r, E=24, n=3)  # n=2 --> a=8.623; n=3 --> 12.9
    # atomic_row_spacing(x=d_c60_n1, R=r, E=24, n=1)  # n=1 --> a=10.08  # duh, because I used it as my reference!
    # atomic_row_spacing(x=d_c60_unknown, R=r, E=24, n=2)  # n=2 --> a=10.7; n=3 --> a=16.1
    #
    # atomic_row_spacing(x=d_c60_n3, R=r2, E=24, n=2)  # n=2 --> a=8.1; n=3 --> 12.2
    #
    # # Double-check using C60/Au(111)
    # fn = r'14_10_20_C60_Au111_1654.tiff'
    # fp = os.path.join(PATH, fn)
    #
    # d_c60_n1_ref = get_line_distance((699, 398), (534, 358))  # n=1, 17eV
    # d_c60_n2_ref = get_line_distance((553, 764), (816, 589))  # n=2, 17eV
    #
    # r2 = distance_to_screen(a=10.08, x=d_c60_n1_ref, E=17, n=1)
    #
    # atomic_row_spacing(x=d_c60_n2_ref, R=r2, E=17, n=2)  # n=2 --> a=10.835, so that "reconstruction" is just 2nd order spots
    #
    # ### Au(111) ###
    #
    # ddir = r'E:\atully\LEED data\2020_January'
    # fn = r'17_01_20_Au111_1819.tiff'
    # fp = os.path.join(ddir, fn)
    #
    # # d_au_n1 = get_line_distance((, ), (, ))  # n=1, 190eV
    # # d_au_n2 = get_line_distance((, ), (, ))  # n=2, 190eV
    #
    # # r = distance_to_screen(a=2.88, x=d_au_n1, E=190, n=1)
    # # atomic_row_spacing(x=d_au_n2, R=r, E=190, n=2)  #

    """Want to click on dot. Have figure update with a point marked and x,y coordinates appear at the bottom"""

##########

# fn = r'14_10_20_C60_Au111_1654.tiff'
#
# d_c60 = get_line_distance((699, 398), (534, 358))  # n=1, 17eV
# d_au = get_line_distance((652, 271), (385, 339))  # n=1, 90eV
#
# r = distance_to_screen(a=2.88, x=d_au, E=90, n=1)
# atomic_row_spacing(x=d_c60, R=r, E=17, n=1)  # a=10.75 (vs 10.07)
#
# r = distance_to_screen(a=10, x=d_c60, E=17, n=1)
# atomic_row_spacing(x=d_au, R=r, E=90, n=1)  # a=2.68 (vs 2.88)


# data = load_image(fn, path=PATH)
# fig = show_image(data)
# data[:][0][0]  # z=[2,4,0], y=0, x=0
# np.array(fig.data[0]).shape
# fig.data

# fig.layout.hovermode = 'closest'  # default
# print(fig.layout.hoverlabel)
# fig.layout.hoverlabel.
# scatter = fig.data
# scatter.on_click()
# print(fig)
# fig.data.
# fig.add_annotation(print(hovertext))
# print(fig.layout.hover)
#
# from plotly.callbacks import Points, InputDeviceState
# points, state = Points(), InputDeviceState()
# f = go.FigureWidget(fig)
# f.on_click()
#
# test = go.FigureWidget([go.Scatter()])
# help(test.data[0].on_click)
#
# def click_fn(trace, points, state):
#     inds = points.point_inds
#     print(inds)
#     print(scatter.selectedpoints)
#     return inds
#
# fig2=go.FigureWidget(go.Scatter(x=[1, 2], y=[3, 0]))
# # fig3=fig.add_trace()
# scatter = fig2.data[0]
# scatter.on_click(click_fn)
# fig2.show(renderer=DEFAULT_RENDERER)
