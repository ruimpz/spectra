import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import spectra as spt

#obs_data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits"))
obs_data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela1.fits"))
synt_data = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/p4250:g+4.0:m0.0:t01:z-1.00:a+0.40.GES4750.fits"))

# Read line library
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]
wls, EP, lgf, EW = np.array(df["lambda"]), np.array(df["EP"]), np.array(df["loggf"]), np.array(df["EW"])

EP1 = [2.1, 2.3]
EP2 = [4.6, 4.7]

# Get lines for temperature estimation
wl = wls[(wls > obs_data[0][0]) & (wls < obs_data[0][-1])]
valid_lines1 = np.where((EP > EP1[0]) & (EP < EP1[1]) & (wls > obs_data[0][0]) & (wls < obs_data[0][-1]))
valid_lines2 = np.where((EP > EP2[0]) & (EP < EP2[1]) & (wls > obs_data[0][0]) & (wls < obs_data[0][-1]))
wl1, lgf1 = wls[valid_lines1], lgf[valid_lines1]
wl2, lgf2 = wls[valid_lines2], lgf[valid_lines2]


print("Estimating temperature....")
T = spt.get_temp_estimate(obs_data, EP1, EP2, wl1, wl2, lgf1, lgf2)
print("Done. T = {}K".format(T))
t_range = spt.get_temp_range(T)

# Read synth spectra library
lib = spt.read_library("../library/GES_UVESRed580_deltaAlphaFe+0.0_fits", t_range)

print("Calculating observed Ws...")
obs_wl, obs_W, inds = spt.get_line_Ws(obs_data, wl, get_inds=True)
comparison_wls = np.where((obs_wl > synt_data[0][0]) & (obs_wl < synt_data[0][-1]))
obs_wl, obs_W = obs_wl[comparison_wls], obs_W[comparison_wls]
print("Done.")

def get_spec_fit(library, wl_obs, W_obs, \
                 library_path="../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/", \
                 lambda_file="../imagens_teste/GES_UVESRed580_Lambda.fits", k = .2):
    """Returns string.

    Compares equivalent width of lines of observed spectra with database of synthetic spectra and returns name of best fitting synthetic spectra.
    """
    fits = np.zeros(len(library))
    for i, synt_spec in enumerate(library):
        synt_data = np.array(spt.read_sintfile(lambda_file, library_path + synt_spec))
        wl_calculated, W_synt, inds = spt.get_line_Ws(synt_data, wl_obs, get_inds=True, k = k)
        p, o = opt.curve_fit(lambda x, a, b: a*x+b, W_synt, W_obs[inds])
        fits[i] = p[0]
    return fits

print("Fitting spectra...")
f = get_spec_fit(lib, obs_wl, obs_W, k = .16)
print("Done.")
print(np.min(np.absolute(f-1)))
best_fit = lib[np.argmin(np.absolute(f-1))]
print(best_fit)

best_data = spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/{}".format(best_fit))

l1, l2 = 5700, 5720
rbest_data = spt.limit_spec(best_data, l1, l2)
robs_data = spt.limit_spec(obs_data, l1, l2)
plt.plot(robs_data[0], robs_data[1]/1.4)
plt.plot(rbest_data[0], rbest_data[1])
plt.show()


#l1, l2 = 6000, 6030
#rbest_data = spt.limit_spec(best_data, l1, l2)
#robs_data = spt.limit_spec(obs_data, l1, l2)
#plt.plot(robs_data[0], robs_data[1])
#plt.plot(rbest_data[0], rbest_data[1])
#plt.show()

#l1, l2 =6300, 6330
#rbest_data = spt.limit_spec(best_data, l1, l2)
#robs_data = spt.limit_spec(obs_data, l1, l2)
#plt.plot(robs_data[0], robs_data[1])
#plt.plot(rbest_data[0], rbest_data[1])
#plt.show()
