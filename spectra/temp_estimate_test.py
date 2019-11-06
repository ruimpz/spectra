import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import spectra as spt

line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
wl, EP, lgf, EW = line_df["lambda"], line_df["EP"], line_df["loggf"], line_df["EW"]
data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela1.fits"))

wl1 = wl[(EP > EP1) & (EP < EP2) & (wl > data[0,0]) & (wl < data[0, -1])]
wl2 = wl[(EP > EP3) & (EP < EP4) & (wl > data[0,0]) & (wl < data[0, -1])]
lgf1 = lgf[(EP > EP1) & (EP < EP2) & (wl > data[0,0]) & (wl < data[0, -1])]
lgf2 = lgf[(EP > EP3) & (EP < EP4) & (wl > data[0,0]) & (wl < data[0, -1])]
EW1 = EW[(EP > EP1) & (EP < EP2) & (wl > data[0,0]) & (wl < data[0, -1])]
EW2 = EW[(EP > EP3) & (EP < EP4) & (wl > data[0,0]) & (wl < data[0, -1])]

ep1 = [2.1, 2.3]
ep2 = [4.5, 4.6]
T = spt.get_temp_estimate(data, ep1, ep2, wl1, wl2, lgf1, lgf2)
print(T)

W1 = np.empty(wl1.size)
W2 = np.empty(wl2.size)
k = .2
for i, wls in enumerate(wl1):
    data_l = spt.limit_spec(data, wls - k, wls + k)
#    ps = spt.ajust_gauss(data_l)
#    x = np.linspace(wls - k, wls + k, 1000)
#    plt.plot(data_l[0], data_l[1])
#    y = spt.gauss(x, *ps)
#    plt.plot(x, y)
    W1[i] = spt.get_W(data_l)
for i, wls in enumerate(wl2):
    data_l = spt.limit_spec(data, wls - .2, wls + .2)
    W2[i] = spt.get_W(data_l)


x1 = np.array(lgf1 + np.log10(wl1))
y1 = np.absolute(np.array(W1*1000/wl1))
x2 = np.array(lgf2 + np.log10(wl2))
y2= np.absolute(np.array(W2*1000/wl2))

p1, o = opt.curve_fit(lambda x, a, b: a*x + b, x1, y1)
p2, o = opt.curve_fit(lambda x, a, b: a*x + b, x2, y2)

ys = np.linspace(max(y1[0], y2[0]), min(y1[-1], y2[-1]), 100)

d = np.absolute(np.mean((ys * (p2[0] - p1[0]) + p1[0]*p2[1] - p1[1]*p2[0]) / p1[0] / p2[0]))
T = np.absolute(5040 * (4.55 - 2.2)) / d







