import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import spectra as spt

#data0 = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits"))
data0 = np.array(spt.read_obsfile("../estrelas_a_analisar/estrela1.fits"))
data1 = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/p4250:g+4.0:m0.0:t01:z-1.00:a+0.40.GES4750.fits"))
#data1 = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/p5750:g+3.5:m0.0:t01:z-1.00:a+0.40.GES4750.fits"))
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]
wls, EP, lgf, EW = np.array(df["lambda"]), np.array(df["EP"]), np.array(df["loggf"]), np.array(df["EW"])
EP1 = [2.1, 2.3]
EP2 = [4.5, 4.6]

#limit spectra to usable part
#data1 = spt.limit_spec(data1, 4800, 6800)
#g = spt.get_conv_gauss(wls, 50000)
#data1 = spt.conv(data1, g)
#plt.plot(data1[0], data1[1])
#plt.show()
#wls = wls[(wls > data0[0,0]) & (wls < data0[0, -1])]

#get line data and limit to obs spectra
valid_lines1 = np.where((EP > EP1[0]) & (EP < EP1[1]) & (wls > data0[0,0]) & (wls < data0[0, -1]))
valid_lines2 = np.where((EP > EP2[0]) & (EP < EP2[1]) & (wls > data0[0,0]) & (wls < data0[0, -1]))
wl1 = wls[valid_lines1]
wl2 = wls[valid_lines2]
lgf1 = lgf[valid_lines1]
lgf2 = lgf[valid_lines2]
EW1 = EW[valid_lines1]
EW2 = EW[valid_lines2]


wls = wls[(wls > data1[0, 0]) & (wls < data1[0, -1])]



print("Estimating temperature...")
T0 = spt.get_temp_estimate(data0, EP1, EP2, wl1, wl2, lgf1, lgf2)
t_range = spt.get_temp_range(T0)
print("Done")
print("T = ", T0)


lib = spt.read_library("../library/GES_UVESRed580_deltaAlphaFe+0.0_fits", t_range)

print("Calculating observed Ws...")
w0, W0 = spt.get_line_Ws(data0, wls)

print("Done")

data1 = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/{}".format(lib[0])))
w0, W0 = spt.get_line_Ws(data1, wls)
plt.plot(data1[0], data1[1])
plt.show()


#diffs = np.zeros(len(lib))
#for i, fit in enumerate(lib):
#    data1 = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/{}".format(fit)))
#    W1 = spt.get_line_Ws(data1, w0, k = .15, get_zeros=True, get_inds=True)
#    diff = (W0*1000 - W1*1000)**2
#    diffs[i] = diff.sum()

print("Calculating best fit...")
best_fit = spt.get_spec_fit(lib, w0, W0, k = .2)
best_fit = "p5500:g+5.0:m0.0:t01:z-1.00:a+0.40.GES4750.fits"
print("Done")
#
print("T = ", T0)
print("Best fit: ", best_fit)
datab = np.array(spt.read_sintfile("../imagens_teste/GES_UVESRed580_Lambda.fits", "../library/GES_UVESRed580_deltaAlphaFe+0.0_fits/{}".format(best_fit)))
#
#representation range
l1, l2 = 5200, 5240
#
rdata0 = spt.limit_spec(data0, l1, l2)
rdatab = spt.limit_spec(datab, l1, l2)

rdata = spt.process_synt(rdata0, rdatab, 50000)

plt.plot(rdata0[0], rdata0[1]/1.4, label="Observado")
plt.plot(rdata[0], rdata[1], label = "Sintético")
plt.legend()
plt.show()
#
#l1, l2 = 5600, 5640
#
#rdata0 = spt.limit_spec(data0, l1, l2)
#rdatab = spt.limit_spec(datab, l1, l2)
#
#rdata = spt.process_synt(rdata0, rdatab, 50000)
#
#plt.plot(rdata0[0], rdata0[1]/1.4, label="Observado")
#plt.plot(rdata[0], rdata[1], label = "Sintético")
#plt.legend()
#plt.show()
#
#l1, l2 = 6200, 6240
#
#rdata0 = spt.limit_spec(data0, l1, l2)
#rdatab = spt.limit_spec(datab, l1, l2)
#
#rdata = spt.process_synt(rdata0, rdatab, 50000)
#
#plt.plot(rdata0[0], rdata0[1]/1.4, label="Observado")
#plt.plot(rdata[0], rdata[1], label = "Sintético")
#plt.legend()
#plt.show()
