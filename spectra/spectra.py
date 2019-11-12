from scipy import signal as sg
from scipy import optimize as opt
from astropy.io import fits
import numpy as np
import os

# IO functions for fits files

def read_obsfile(file_name):
    """Return a 2 numpy arrays

    Reads a fits file of observed flux, generates the corresponding wavelength from file header and returns
    (wavelength, flux) data.
    """
    hdu = fits.open(file_name)[0]
    data = hdu.data
    wl = np.arange(data.size) * hdu.header["cdelt1"] + hdu.header["crval1"]
    return wl, data

def read_sintfile(wl_file, flux_file):
    """Return 2 numpy arrays

    Reads a simulated wavelength file and a simulated flux file and returns (wavelength, flux) data.
    """
    hdu_wl = fits.open(wl_file)[0]
    hdu_flux = fits.open(flux_file)[0]
    wl = hdu_wl.data
    flux = hdu_flux.data
    return wl, flux


# Data analysis

def limit_spec(data, lambda_min, lambda_max):
    """Returns 2 numpy arrays

    Reads (wavelength, flux) data and limits it's range by wavelength.
    """
    data1, data2 = data
    wl = data1[(data1>lambda_min)&(data1<lambda_max)]
    flux = data2[(data1>lambda_min)&(data1<lambda_max)]
    return wl, flux

def gauss(x, a, b, c,d):
    """Returns numpy array

    Gaussian generator.
    """
    return a*np.exp( - (x - b)**2 / 2 / c**2 ) + d


def W_gauss(a, b, c, d):
    """Returns float

    Calculates equivalent width of line W_lambda from analytical formula (gaussian approximation).
    """
    return np.absolute(np.sqrt(2 * np.pi) * a * c / d)


def ajust_gauss(data):
    """Returns list

    Takes spectra data, minimum and maximum wavelengths and returns parameters of gaussian fit
    """
    params0=[-1, data[0][data[1].argmin()], .1, 1]
    popt, r = (0, 0, 0, 0), 0
    while r < .9 and params0[2] > 0:
        print("iteration")
        try:
            popt, pcov = opt.curve_fit(gauss, *data, params0)
            y = gauss(data[0], *popt)
            s = (data[1] - y)**2/data[1]
            r = 1 - s.sum()
            params0[2] -= .02
        except RuntimeError:
            print("Error")
            break
    return popt, r



def get_W(data):
    """Returns float

    Takes spectra data, minimum and maximum wavelengths of line limit and return equivalent width W_lambda.
    """
    ps, r = ajust_gauss(data)
    print("r = ",r)
    return W_gauss(*ps)



def get_temp_estimate(data, EP1, EP2, wls1, wls2, lgf1, lgf2, k=.4, N = 1000):
    """Returns float

    Takes observed spectra data, EP ranges, corresponding wavelength and log(gf) lists and estimates temperature.
    Takes optional parameters k and N for line width estimation and number of points in linear regression, respectively.
    """
    W1, W2 = np.empty(wls1.size), np.empty(wls2.size)
    for i, wl in enumerate(wls1):
        line = limit_spec(data, wl - k, wl + k)
        W1[i] = get_W(line)
    for i, wl in enumerate(wls2):
        line = limit_spec(data, wl - k, wl + k)
        W2[i] = get_W(line)
    x1 = np.array(lgf1 + np.log10(wls1))
    y1 = np.absolute(np.array(W1*1000/wls1))
    x2 = np.array(lgf2 + np.log10(wls2))
    y2 = np.absolute(np.array(W2*1000/wls2))
    p1, o1 = opt.curve_fit(lambda x, a, b: a*x+b, x1, y1)
    p2, o2 = opt.curve_fit(lambda x, a, b: a*x+b, x2, y2)
    ymin, ymax = max(y1[0], y2[0]), min(y1[-1], y2[-1])
    ys= np.array([ymin, ymax])
    d = np.absolute(np.mean((ys * (p2[0] - p1[0]) + p1[0]*p2[1] - p1[1]*p2[0]) / p1[0] / p2[0]))
    return np.absolute(5040 * (np.mean(EP2) - np.mean(EP1))) / d



def get_temp_range(T, delta_T=400):
    Ts = np.arange(4000, 7250, 250)
    return Ts[(Ts > T - delta_T) & ( Ts < T + delta_T)]


def read_library(library_path, Ts, logg_min = 3.5):
    data = []
    for file in os.listdir(library_path):
        for T in Ts:
            if file.startswith("p{}".format(T)):
                data.append(file)
    return data


