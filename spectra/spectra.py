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
    wl = data1[(data1>=lambda_min)&(data1<=lambda_max)]
    flux = data2[(data1>=lambda_min)&(data1<=lambda_max)]

    return np.array([wl, flux])

def gauss(x, a, b, c, d):
    """Returns numpy array

    Gaussian generator.
    """
    return a*np.exp( - (x - b)**2 / 2 / c**2 ) + d

def W_gauss(a, b, c, d, normed = True):
    """Returns float

    Calculates equivalent width of line W_lambda from analytical formula (gaussian approximation).
    """
    if normed:
        return np.absolute(np.sqrt(2 * np.pi) * a * c)
    else:
        return np.absolute(np.sqrt(2 * np.pi) * a * c) / d

def get_exp_params(wl, R):

    """Returns float.

    Takes wavelength array from synthetic spectra and returns parameters for simulating experimental resolution on said spectra
    """
    wl_med = (wl[0] + wl[-1]) / 2
    dl = wl_med / R
    return dl / (2 * np.sqrt(2 * np.log(2)))


def get_conv_gauss(wl, R):
    """Returns array.

    Takes wavelength array from synthetic spectra and returns gaussian for convolution.
    """
    c = get_exp_params(wl, R)
    dx = wl[1] - wl[0]
    x = np.arange(-5*c, 5*c + dx, dx)
    y = gauss(x, 1, 0, c, 0)
    return y / y.sum()


def conv(data, g):
    """Returns array.

    Takes synthetic spectra and returns the result of it's convolution with function g.
    """
    data[1] =  sg.convolve(data[1], g, mode="same")
    return data


def interp(obs_data, synt_data):
    """Returns array.

    Takes observed and synthetic spectra and returns flux of synthetic spectra interpolated with observed data.
    """
    return np.interp(obs_data[0], synt_data[0], synt_data[1])

def get_rot_G(obs_data, v, mins, e = .6):
    x = np.linspace(-.5, .5, 1000)
    sun_limbo = np.array([.66, 1.162, 1.661])
    lm = np.mean((sun_limbo / mins))
    lm /= 1-2
    G = ( 2*(1-e)*np.sqrt(1-(x/lm)**2)+np.pi*e/2*(1-(x/lm)**2) )/( np.pi*lm*(1-e/3) )
    G[np.isnan(G)] = 0
    return G/np.sum(G)


def process_spectra(l1, l2, obs_data, synt_data, wls, R, rotation=False):

    """Returns array.

    Takes observed and synthetic spectra and returns synthetic spectra after application of convolution with gaussian and interpolation with observed spectra data.
    """
    if rotation:
        v, mins = estimate_vsinI(obs_data, wls)
        G = get_rot_G(synt_data, v, mins)
    obs_data = limit_spec(obs_data, l1, l2)
    synt_data = limit_spec(synt_data, l1, l2) 
    x = np.copy(synt_data)
    g = get_conv_gauss(x[0], R)
    x = conv(x, g)
    if rotation:
        x = conv(x, G)
    return np.array([obs_data[0], interp(obs_data, x)]), obs_data


def limit_line(line, wl, k = .2):
    try:
        ind_min = line[1].argmin()
        wl_min = line[0, ind_min]
        step1 = line[1, :ind_min-5] - line[1, 1:ind_min-4]
        ind1 = np.where(step1 < 0)
        if ind1[0].size>0:
            k1 = wl_min - line[0, ind1[0][-1]]
        else:
            k1 = k
        step2 = line[1, ind_min:-5] - line[1, ind_min+1:-6]
        ind2 = np.where(step2 > 0)
        if ind2[0].size>0:
            k2 = line[0, ind2[0][0] + ind_min] - wl_min
        else:
            k2 = k
        return limit_spec(line, wl_min - k1, wl_min + k2), wl_min
    except ValueError:
        return limit_spec(line, wl_min - k, wl_min + k), wl_min


def ajust_gauss(data, wl_min, std_dev=.15):
    """Returns list

    Takes spectra data, minimum and maximum wavelengths and returns parameters of gaussian fit
    """
    params0=[-1, wl_min, .1, 1]

    popt, r = (0, 0, 0, 0), 0
    while r < .9 and params0[2] > 0:
        try:
            popt, pcov = opt.curve_fit(gauss, *data, params0)
            y = gauss(data[0], *popt)
            s = (1 - y/data[1])**2
            r = 1 - s.sum()
            params0[2] -= .01
        except RuntimeError or TypeError:
            print("Error. Line ignored.")
            params0[2] -= .01
    return popt, r


def get_W(data, wl_min, r_tol = 0.96, normed = True):
    """Returns float

    Takes spectra data, minimum and maximum wavelengths of line limit and return equivalent width W_lambda.
    """
    ps, r = ajust_gauss(data, wl_min)
    if r > r_tol:
        return W_gauss(*ps, normed = normed), ps, r
    else:
#        print("Bad line fit -> ignored")
        return 0, [0, 0, 0, 0], 0


def get_line_Ws(data, ws, k = .2, limit = False, get_inds = False, get_zeros=False, plot = False, count=False, normed=True):
    W = np.zeros(len(ws))
    if plot:
        for i, wl in enumerate(ws):
            line = limit_spec(data, wl-k, wl+k)
            if limit:
                line, wl_min = limit_line(line, wl, k = k)
            else:
                wl_min = line[0, line[1].argmin()]
            W[i], ps, r = get_W(line, wl_min, normed=normed)
            plt.plot(line[0], line[1], label="Observed line")
            plt.plot(line[0], gauss(line[0], *ps), label="Gaussian fit")
            plt.xlabel(r"$\lambda$ $(\AA)$")
            plt.ylabel("Flux (normalized)")
            plt.legend()
            plt.show()
    else:
        for i, wl in enumerate(ws):
            line = limit_spec(data, wl-k, wl+k)
            if limit:
                line, wl_min = limit_line(line, wl, k = k)
            else:
                wl_min = line[0, line[1].argmin()]
            W[i], ps, r = get_W(line, wl_min, normed=normed)
    inds = np.where(W != 0)
    if count:
        print(len(W) - len(inds[0]))
    if get_zeros:
        return W
    elif get_inds:
        return ws[inds], W[inds], inds
    else:
        return ws[inds], W[inds]


def get_temp_estimate(data, EP1, EP2, wls1, wls2, lgf1, lgf2, k = .2):
    """Returns float

    Takes observed spectra data, EP ranges, corresponding wavelength and log(gf) lists and estimates temperature.
    """
    wls1, W1, inds1 = get_line_Ws(data, wls1, get_inds=True, k = k, normed = False)
    wls2, W2, inds2 = get_line_Ws(data, wls2, get_inds=True, k = k, normed = False)
    lgf1 = lgf1[inds1]
    lgf2 = lgf2[inds2]
    x1 = np.array(lgf1 + np.log10(wls1))
    y1 = np.array(np.log10(W1*1000/wls1))
    x2 = np.array(lgf2 + np.log10(wls2))
    y2 = np.array(np.log10(W2*1000/wls2))
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

#def get_spec_fit(library, wl_obs, W_obs, \
#                 library_path="../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/", \
#                 lambda_file="../imagens_teste/GES_UVESRed580_Lambda.fits", k = .2):

#
#    Compares equivalent width of lines of observed spectra with database of synthetic spectra and returns name of best fitting synthetic spectra.
#    """
#    diffs = np.zeros(len(library))
#    for i, synt_spec in enumerate(library):
#        synt_data = np.array(read_sintfile(lambda_file, library_path + synt_spec))
#        W_synt = get_line_Ws(synt_data, wl_obs, get_zeros=True, k = k, limit=True, count=True)
#        diff = (W_obs- W_synt)**2
#        diffs[i] = diff.sum()
#    return library[np.argmin(diffs)]


def get_spec_fit(library, wl_obs, W_obs, \
                 library_path="../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/", \
                 lambda_file="../imagens_teste/GES_UVESRed580_Lambda.fits", k = .2):
    """Returns string.

    Compares equivalent width of lines of observed spectra with database of synthetic spectra and returns name of best fitting synthetic spectra.
    """
    fits = np.zeros(len(library))
    for i, synt_spec in enumerate(library):
        try:
            synt_data = np.array(read_sintfile(lambda_file, library_path + synt_spec))
            wl_calculated, W_synt, inds = get_line_Ws(synt_data, wl_obs, get_inds=True, k = k, limit=True)
            if W_synt.size > 2:
                p, o = opt.curve_fit(lambda x, a, b: a*x+b, W_synt, W_obs[inds])
            else:
                p = [0, 0]
            fits[i] = p[0]
        except TypeError:
            fits[i] = 0
    return fits


def get_fft_min(line, wl):
    c = 2.99e3
    min_freqs = np.zeros(3)
    sun_limbo = np.array([.66, 1.162, 1.661]) 
    j = 0
    line[1] /= np.max(line[1])
    ext_line = np.concatenate((line[1], np.zeros(10000)))
    f_line = np.fft.fft(ext_line)
    f_freq = np.fft.fftfreq(ext_line.size)
    f_line = np.absolute(f_line)
    inds = np.where(f_freq > 1e-6)

    f_freq, f_line = f_freq[inds], f_line[inds]
    f_line /= np.max(f_line)
    for i in range(1,f_line.size-1):
        if j < min_freqs.size and f_line[i-1] > f_line[i] and f_line[i+1] > f_line[i]:
            min_freqs[j] = f_freq[i]
            j += 1
    vsini = sun_limbo * c / min_freqs / wl
    return np.mean(vsini[1:]), min_freqs



def estimate_vsinI(obs_data, wls, k = .2):
    vs = np.zeros(wls.size)
    minima = np.zeros((wls.size, 3))
    for i, wl in enumerate(wls):
        line = limit_spec(obs_data, wl-k, wl+k)
        vs[i], minima[i, :] = get_fft_min(line, wl)
    minima = np.sum(minima, axis= 0) / wls.size
    return np.mean(vs), minima
