import matplotlib.pyplot as plt
import spectra as spt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#data = spt.read_obsfile("../estrelas_a_analisar/estrela2_vcor.fits")
data = spt.read_obsfile("../estrelas_a_analisar/estrela1.fits")
line_df = pd.read_csv("~/University/Comp_Astro/Project_1/line_list_tsantaki.dat", sep="\t")
df = line_df[line_df.El.str.contains("FeI ")]

# Read relevant wavelengths and limit them to spectra range
wls = np.array(df["lambda"])
wls = wls[(wls > data[0][0]) & (wls < data[0][-1])]

k = .25
n = 31
line = spt.limit_spec(data, wls[n]-k, wls[n]+k)

plt.plot(line[0], line[1])
plt.show()

#nline = line[1] / np.max(line[1])

#pad_line = np.concatenate((nline, np.zeros(10000)))
#fline = np.fft.fft(pad_line)
#freq = np.fft.fftfreq(pad_line.size)
#inds = np.where(freq > 1e-8)
#gline = np.absolute(fline[inds])
#
#n_max = 4
#minima = np.zeros(n_max)
#N = 0
#for i in range(1, gline.size-1):
#    if N < n_max:
#        if gline[i-1] > gline[i] and gline[i+1] > gline[i]:
#            minima[N] = i
#            N += 1
#
#dl = wls[1] - wls[0]
#e = .6
#K = (line[0] - wls[n])/dl
#s = 2*np.pi*q
#g = 2*((1-e)*)
#
#plt.plot(gfreq[indG], np.absolute(gg))
#plt.plot(freq[inds], gline)
#plt.show()

