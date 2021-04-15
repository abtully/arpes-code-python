import numpy as np


def fourier_2d(data: np.ndarray, xaxis: np.ndarray, yaxis: np.ndarray):
    fft = np.fft.fft2(data)
    lx, ly = np.shape(xaxis)[0], np.shape(yaxis)[0]
    dx = (xaxis[-1] - xaxis[0]) / lx
    dy = (yaxis[-1] - yaxis[0]) / ly
    x_fft = 1 / dx * np.linspace(-0.5, 0.5, lx)
    y_fft = 1 / dy * np.linspace(-0.5, 0.5, ly)
    dat_fft = 2 * np.abs(np.fft.fftshift(fft))
    return dat_fft, x_fft, y_fft