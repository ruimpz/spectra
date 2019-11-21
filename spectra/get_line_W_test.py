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
for i in range(len(wls)):
    line = spt.limit_spec(data, wls[i]-k, wls[i]+k)
    plt.plot(line[0], line[1])
    plt.show()
    ind_min = line[1].argmin()
    step1 = line[1, :ind_min-1] - line[1, 1:ind_min]
    ind1 = np.where(step1 < 0)
    if ind1[0].size>0:
        k1 = line[0, ind1[0][-1]]
    else:
        k1 = wls[i] - k
    step2 = line[1, ind_min:-1] - line[1, ind_min +1:]
    ind2 = np.where(step2 > 0)
    if ind2[0].size>0:
        k2 = line[0, ind2[0][0] + ind_min]
    else:
        k2 = wls[i] + k
    line = spt.limit_spec(line, k1, k2)
    ps, r = spt.ajust_gauss(line, line[0, ind_min])
    print(r)
    plt.plot(line[0], line[1], label="Observed spectra")
    plt.plot(line[0], spt.gauss(line[0], *ps), label="Gaussian fit")
    plt.ylabel("Flux (normalized)")
    plt.xlabel(r"$\lambda$ ($\AA$) ")
    plt.legend()
    plt.show()
