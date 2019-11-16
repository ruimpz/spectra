from scipy import signal as sg
from scipy import optimize as opt
from astropy.io import fits
import matplotlib.pyplot as plt
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
    return np.array([wl, flux])

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
    params0=[-1, data[0][data[1].argmin()], .1, data[0][data[1].argmax()]]
    popt, r = (0, 0, 0, 0), 0
    while r < .9 and params0[2] > 0:
        try:
            popt, pcov = opt.curve_fit(gauss, *data, params0)
            y = gauss(data[0], *popt)
            s = (1 - y/data[1])**2
            r = 1 - s.sum()
            params0[2] -= .02
        except RuntimeError:
            print("Error")
            break
    return popt, r



def get_W(data, r_tol = .994):
    """Returns float

    Takes spectra data, minimum and maximum wavelengths of line limit and return equivalent width W_lambda.
    """
    ps, r = ajust_gauss(data)
    if r > r_tol:
        return W_gauss(*ps), ps, r
    else:
        return 0, 0, 0

def get_line_Ws(data, ws, k = .2, get_inds = False, get_zeros=False):
    W = np.zeros(len(ws))
    for i, wl in enumerate(ws):
        line = limit_spec(data, wl-k, wl+k)
        W[i], p, r = get_W(line)
#        if r > .98:
#            print(r)
#           plt.plot(line[0], line[1])
#            plt.plot(line[0], gauss(line[0], *p))
#            plt.show()
    inds = np.where(W != 0)
    if get_zeros:
        return W
    elif get_inds:
        return ws[inds], W[inds], inds
    else:
        return ws[inds], W[inds]


def get_temp_estimate(data, EP1, EP2, wls1, wls2, lgf1, lgf2):
    """Returns float

    Takes observed spectra data, EP ranges, corresponding wavelength and log(gf) lists and estimates temperature.
    Takes optional parameters k and N for line width estimation and number of points in linear regression, respectively.
    """
    wls1, W1, inds1 = get_line_Ws(data, wls1, get_inds=True)
    wls2, W2, inds2 = get_line_Ws(data, wls2, get_inds=True)
    lgf1 = lgf1[inds1]
    lgf2 = lgf2[inds2]
    x1 = np.array(lgf1 + np.log10(wls1))
    y1 = np.array(W1*1000/wls1)
    x2 = np.array(lgf2 + np.log10(wls2))
    y2 = np.array(W2*1000/wls2)
#    plt.plot(x1, y1, "+")
#    plt.plot(x2, y2, "+")
#    plt.show()
    p1, o1 = opt.curve_fit(lambda x, a, b: a*x+b, x1, y1)
    p2, o2 = opt.curve_fit(lambda x, a, b: a*x+b, x2, y2)
    ymin, ymax = max(y1[0], y2[0]), min(y1[-1], y2[-1])
    ys= np.array([ymin, ymax])
    d = np.absolute(np.mean((ys * (p2[0] - p1[0]) + p1[0]*p2[1] - p1[1]*p2[0]) / p1[0] / p2[0]))
    return np.absolute(5040 * (np.mean(EP2) - np.mean(EP1))) / d



def get_temp_range(T, delta_T=400):
    """Returns array

    Takes a temperature estimate and generates list of possible temperature values within an error detla_T to look up on the database.
    """
    Ts = np.arange(4000, 7250, 250)
    return Ts[(Ts > T - delta_T) & ( Ts < T + delta_T)]


def read_library(library_path, Ts, logg_min = 3.5):
    """Returns list
    """
    data = []
    for file in os.listdir(library_path):
        for T in Ts:
            if file.startswith("p{}".format(T)):
                data.append(file)
    return data

def get_spec_fit(library, wl_obs, W_obs, \
                 library_path="../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/", \
                 lambda_file="../imagens_teste/GES_UVESRed580_Lambda.fits"):
    """Returns string.

    Compares equivalent width of lines of observed spectra with database of synthetic spectra and returns name of best fitting synthetic spectra.
    """
    diffs = np.zeros(len(library))
    for i, synt_spec in enumerate(library):
        synt_data = np.array(read_sintfile(lambda_file, library_path + synt_spec))
        W_synt = get_line_Ws(synt_data, wl_obs, get_zeros=True)
        diff = (W_obs*1000 - W_synt*1000)**2
        diffs[i] = diff.sum()
    return library[np.argmin(diffs)]

