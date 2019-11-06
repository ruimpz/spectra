import spectra as sp
import numpy as np
import matplotlib.pyplot as plt

data = sp.read_sintfile("/home/flower/University/Comp_Astro/Project_1/imagens_teste/GES_UVESRed580_Lambda.fits",
"/home/flower/University/Comp_Astro/Project_1/imagens_teste/p6250-g+5.0-m0.0-t01-z+0.75-a+0.00.GES4750.fits")
w1, w2 = 5007.8, 5008.2
data = sp.limit_spec(data, w1, w2)
ps = sp.ajust_gauss(data)
x = np.linspace(w1, w2, 1000)
W = sp.get_W(data)
print(W)
a, b, c, d = ps
print(a, b, c, d)
plt.plot(data[0], data[1])
plt.plot(x, sp.gauss(x, a, b, c, d))
plt.show()
