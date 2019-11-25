import spectra as spt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize as opt

# Read synthetic spectra
data = spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../imagens_teste/p6250-g+5.0-m0.0-t01-z+0.75-a+0.00.GES4750.fits")
#data = spt.read_obsfile("../estrelas_a_analisar/estrela1.fits")

# Read FeI line library
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]

# Read relevant wavelengths and limit them to spectra range
wls = np.array(df["lambda"])
wls = wls[(wls > data[0][0]) & (wls < data[0][-1])]

k = .2

for i, wl in enumerate(wls):
    line = spt.limit_spec(data, wl - k, wl + k)
    plt.plot(line[0], line[1])
    plt.xlabel(r"$\lambda$ $(\AA)$")
    plt.ylabel("Flux (normalized)")
    plt.show()
