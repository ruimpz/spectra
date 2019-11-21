import spectra as spt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as opt

obs_data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela1.fits"))
#obs_data = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits"))

# Read FeI line library
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]
wls, EP, lgf, EW = np.array(df["lambda"]), np.array(df["EP"]), np.array(df["loggf"]), np.array(df["EW"])
EP1 = [2.1, 2.3]
EP2 = [4.6, 4.7]

# Limit wavelength lists to selected multiplets
valid_lines1 = np.where((EP > EP1[0]) & (EP < EP1[1]) & (wls > obs_data[0,0]) & (wls < obs_data[0, -1]))
valid_lines2 = np.where((EP > EP2[0]) & (EP < EP2[1]) & (wls > obs_data[0,0]) & (wls < obs_data[0, -1]))
#valid_lines1 = np.where((EP > EP1[0]) & (EP < EP1[1]))
#valid_lines2 = np.where((EP > EP2[0]) & (EP < EP2[1]))
wl1 = wls[valid_lines1]
wl2 = wls[valid_lines2]
lgf1 = lgf[valid_lines1]
lgf2 = lgf[valid_lines2]
EW1 = EW[valid_lines1]
EW2 = EW[valid_lines2]

# get Ws
fitted_wl1, W1 = spt.get_line_Ws(obs_data, wl1)
fitted_wl2, W2 = spt.get_line_Ws(obs_data, wl2)

# Get temperature estimation
T0 = spt.get_temp_estimate(obs_data, EP1, EP2, wl1, wl2, lgf1, lgf2)
print("T = {} K".format(T0))

# We replicate the get_temp_estimate method for representation
x1 = np.array(lgf1 + np.log10(wl1))
y1 = np.array(np.log10(EW1/wl1))
x2 = np.array(lgf2 + np.log10(wl2))
y2 = np.array(np.log10(EW2/wl2))
p1, o1 = opt.curve_fit(lambda x, a, b: a*x+b, x1, y1)
p2, o2 = opt.curve_fit(lambda x, a, b: a*x+b, x2, y2)
ymin, ymax = max(y1[0], y2[0]), min(y1[-1], y2[-1])
ys= np.array([ymin, ymax])
d = np.absolute(np.mean((ys * (p2[0] - p1[0]) + p1[0]*p2[1] - p1[1]*p2[0]) / p1[0] / p2[0]))
print("T_sun = ", np.absolute(5040 * (np.mean(EP2) - np.mean(EP1))) / d)



plt.plot(x1, y1, "+", label="EP1 = [2.1, 2.3]")
plt.plot(x2, y2, "+", label="EP2 = [4.5, 4.6]")
plt.plot(x1, p1[0]*x1 + p1[1], label="EP1 fit")
plt.plot(x2, p2[0]*x2 + p2[1], label="EP2 fit")
plt.legend()
plt.title("Sun (T = 5769K)")
plt.xlabel(r"$\log \left( gf \lambda \right)$")
plt.ylabel(r"$\log \left( \frac{W}{\lambda} \right)$")
plt.show()

