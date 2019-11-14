import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import spectra as spt

data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits"))
#data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela1.fits"))

line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]

EP1 = [2.1, 2.3]
EP2 = [4.5, 4.6]
wls, EP, lgf, EW = df["lambda"], df["EP"], df["loggf"], df["EW"]

wl1 = wls[(EP > EP1[0]) & (EP < EP1[1]) & (wls > data[0,0]) & (wls < data[0, -1])]
wl2 = wls[(EP > EP2[0]) & (EP < EP2[1]) & (wls > data[0,0]) & (wls < data[0, -1])]
lgf1 = lgf[(EP > EP1[0]) & (EP < EP1[1]) & (wls > data[0,0]) & (wls < data[0, -1])]
lgf2 = lgf[(EP > EP2[0]) & (EP < EP2[1]) & (wls > data[0,0]) & (wls < data[0, -1])]
EW1 = EW[(EP > EP1[0]) & (EP < EP1[1]) & (wls > data[0,0]) & (wls < data[0, -1])]
EW2 = EW[(EP > EP2[0]) & (EP < EP2[1]) & (wls > data[0,0]) & (wls < data[0, -1])]


k = .2
wls = wls[(wls>data[0][0]) & (wls<data[0][-1])]
W = np.zeros(wls.size)
for i, wl in enumerate(wls):
    line = spt.limit_spec(data, wl - k, wl + k)
    ps, r = spt.ajust_gauss(line)
    if r > .8:
        W[i] = spt.get_W(line)

good_wls = wls[W > 0.]
good_W = W[W != 0.]


#T = spt.get_temp_estimate(data, EP1, EP2, wl1, wl2, lgf1, lgf2)
#print(T)
#t_range = spt.get_temp_range(T)
#print(t_range)

