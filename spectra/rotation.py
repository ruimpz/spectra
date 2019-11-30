import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import spectra as spt

plt.style.available
plt.style.use('seaborn-paper')

obs_data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits"))
sint_data = spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/p6000:g+5.0:m0.0:t01:z-0.25:a+0.10.GES4750.fits")
obs_data=spt.limit_spec(obs_data, 5832, 6817)

# Read line library
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]
wls, EP, lgf, EW = np.array(df["lambda"]), np.array(df["EP"]), np.array(df["loggf"]), np.array(df["EW"])

wl = wls[(wls > obs_data[0][0]) & (wls < obs_data[0][-1])]

sint = spt.process_synt(obs_data, sint_data, wl, 50000, rotation=True)

plt.plot(sint[0],sint[1])
plt.title("Rotational profile")
plt.xlabel(r"$\lambda - \lambda_0$ $(\AA)$")
plt.ylabel("Flux (normalized)")
plt.show()
