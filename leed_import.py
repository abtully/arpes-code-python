"""
Low energy electron diffraction (LEED) images import functions (tested on .tiff files)
@author: Alexandra Tully
@date: April 2021
"""

from PIL import Image
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from typing import List, Optional

DEFAULT_RENDERER = 'browser'  # this is a constant of this file

PATH = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\October_2020\LEED\Spinview'


def load_image(filename, path=PATH) -> np.ndarray:  # this function returns a numpy array
    """
    This is to load an image from PATH and convert to numpy array

    example fp: /Volumes/AT/LEED data/March_2020/05_03_20_Au111_C60_1705.tiff

    Args:
        filename (str): full name of image in month folder, necessary parameter because no default value
        path (str): global variable for filepath

    Returns:
         (np.ndarray) a numpy array, type ndarray

    """
    fp = get_filepath(path, filename)  # full path
    img_data = plt.imread(fp)  # take tiff file and converts to numpy array because imported with mpl
    return img_data


def get_filepath(path, filename):
    return os.path.join(path, filename)


def show_image(img_data: np.ndarray, show=True):
    """
    Plots numpy array from load_image function on given axes
    Args:
        img_data: numpy array (a matrix)
        show: optional argument to create new figure

    Returns: plots the numpy array

    """
    fig = go.Figure(go.Image(z=img_data))
    if show:
        fig.show(renderer=DEFAULT_RENDERER)
    return fig


def show_image_mpl(img_data, ax=None):
    """
    Plots numpy array from load_image function on given axes
    Args:
        ax (plt.Axes): optional axes
        img_data: numpy array (a matrix)

    Returns: plots the numpy array

    """
    if ax is None:
        fig, ax = plt.subplots(1)  # 1 row 1 column, could do (2,2) etc. Then call it fig, axs.
    ax.imshow(img_data)  # imshow is what causes the image file to pop up
    ax.axes.set_axis_off()
    return ax


if __name__ == '__main__':
    fn = r'14_10_20_C60_Au111_1654.tiff'

    # fp = os.path.join(PATH + fn)
    # image = Image.open(fp)
    # image.show()
    # imarray = np.array(image)

    data = load_image(fn)
    show_image(data)
