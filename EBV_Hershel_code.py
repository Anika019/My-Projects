from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from decimal import *
import pandas as pd

#Load in the AGN obserbed spectrum
from astropy.table import Table
EBV = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\EBV_tablefit.fits')
EBVdata=EBV[1].data
EBV.close()
gordon=EBVdata.EBV_gordon
lyu=EBVdata.EBV_lyu
lyu0=np.array(lyu)
spectrum=EBVdata.spectrum
spectrum0=np.array(spectrum)

#load w3 from full catalog
W3 = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\W3table.fits')
W3data=W3[1].data
W3.close()
w3values=W3data.W3
w3values0=np.array(w3values)
rec=W3data.REC_NO
rec0=np.array(rec)
df=pd.DataFrame({'w3':w3values0,'rec':rec0})
w3range=[]
for x in spectrum0:
    w3range.append(df.loc[df['rec'] == x, 'w3'].iloc[0])

coldquasar=np.array([475,507,2435,2480,2551,3819,3871,4252,4285,4336,4472,4979,5074,5097,5122])

w3coldquasar=[]
for x in coldquasar:
    w3coldquasar.append(df.loc[df['rec'] == x, 'w3'].iloc[0])

df2=pd.DataFrame({'spec':spectrum0,'lyu':lyu0})
lyurange=[]
for x in coldquasar:
    lyurange.append(df2.loc[df2['spec'] == x, 'lyu'].iloc[0])


# plt.plot(gordon,w3range,'.', color='red')
# plt.title('Fig 1: Reddening_Gordon vs W3 values')
# plt.xlabel('Gordon Template EB-V values')
# plt.ylabel('W3 values (AB magnitude)')
# plt.legend(['Gordon EB-V'],shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.ylim(10,20)
# ax = plt.axes()
# ax.set_facecolor('black') #[0.05,0.4,0.4]
# plt.show()

plt.plot(lyu,w3range,'.', color='red')
plt.plot(lyurange,w3coldquasar,'.', color='white')
plt.title('Fig 1: Reddening_Lyu vs W3 values')
plt.xlabel('Lyu Template E(B-V) values')
plt.ylabel('W3 values (AB magnitude)')
plt.legend(['Lyu E(B-V)'],shadow=True, handlelength=1.5,fontsize=10)
plt.grid(linestyle='dashed')
plt.ylim(10,20)
ax = plt.axes()
ax.set_facecolor('black')
ax.annotate('Cold Quasars', xy=(0.9, 16), xytext=(2, 19),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='white'), color='white')
plt.show()


#
# f3=interpolate.interp1d(spectrum,w3range, bounds_error=False)
# gnew=f3(gordon)
# lnew=f3(lyu)
#
# plt.plot(gnew,w3range,'.', color='red')
# plt.title('Fig: interpolated Reddening Gordon vs W3 values')
# plt.xlabel('Gordon EBV values ()')
# plt.ylabel('W3 values ()')
# plt.legend(('Gordon EBV','W3'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.ylim(10,25)
# plt.show()
#
# plt.plot(lnew,w3range,'.', color='red')
# plt.title('Fig: interpoolated Reddening Lyu vs W3 values')
# plt.xlabel('Lyu EBV values ()')
# plt.ylabel('W3 values ()')
# plt.legend(('Lyu EBV','W3'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.ylim(10,25)
# plt.show()



# plt.plot(spectrum,gordon,'^', color='blue')
# plt.xlabel('Hershel Spectrum')
# plt.ylabel('EB-V values')
# plt.plot(spectrum,lyu,'+', color='red')
# plt.title('Fig 3:Spectrum Vs Reddening Lyu and reddening gordon')
# plt.xlabel('Hershel Spectrum')
# plt.ylabel('EB-V values')
# # plt.text(4000, 15, 'Fit parameters:normalizaton:3.36517990e+02, ebv:-2.75477242e-01', fontsize=8)
# plt.legend(('Gordon','Lyu'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.show()




# f=interpolate.interp1d(spectrum,lyu, bounds_error=False)
# ynew=f(gordon)
#
# plt.plot(ynew,lyu,'.', color='red')
# plt.title('Fig: Reddening Gordon vs Interpolated Reddening Lyu')
# plt.xlabel('Gordon EBV values ()')
# plt.ylabel('Lyu EBV values ()')
# plt.legend(('Reddening curve'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.show()
#
# ff=interpolate.interp1d(spectrum,gordon, bounds_error=False)
# ynew2=f(lyu)
#
# plt.plot(gordon, ynew2,'.', color='red')
# plt.title('Fig: interpolated Reddening Gordon vs Reddening Lyu')
# plt.xlabel('Gordon EBV values ()')
# plt.ylabel('Lyu EBV values ()')
# plt.legend(('Reddening curve'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.show()


# keep this plot for poster
plt.plot(gordon,lyu,'.', color='white')
plt.title('Fig 2: Reddening Gordon Vs Reddening Lyu')
plt.xlabel('E(B-V) Gordon')
plt.ylabel('E(B-V) Lyu')
# plt.text(4000, 15, 'Fit parameters:normalizaton:3.36517990e+02, ebv:-2.75477242e-01', fontsize=8)
plt.legend(['Reddening Curve'],shadow=True, handlelength=1.5,fontsize=10)
plt.grid(linestyle='dashed')
ax = plt.axes()
ax.set_facecolor('black')
plt.show()



#load luminosity from full catalog
lum = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\luminositytable.fits')
lumdata=lum[1].data
lum.close()
specvalues=lumdata.LUMINOSITY_SPEC
specvalues0=np.array(specvalues)
phot=lumdata.LUMINOSITY_PHOT
phot0=np.array(phot)
no0=lumdata.REC_NO
no=np.array(no0)
df2=pd.DataFrame({'spec':specvalues0,'no':no0, 'phot':phot0})
specrange=[]
photrange=[]
for x in spectrum0:
    specrange.append(df2.loc[df2['no'] == x, 'spec'].iloc[0])
    photrange.append(df2.loc[df2['no'] == x, 'phot'].iloc[0])

#LUMINOSITY SPEC
# plt.plot(gordon,specrange,'.', color='red')
# plt.title('Fig 4: Reddening Gordon vs luminosity spec values')
# plt.xlabel('Gordon EB-V values')
# plt.ylabel('Luminosity spec values (Ergs/sec) in Log Scale')
# plt.legend(('Gordon EB-V','luminosity Spec'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.yscale('log')
# #plt.ylim(0,3)
# plt.show()

plt.plot(lyu,specrange,'.', color='red')
plt.title('Fig 3: Reddening Lyu vs Luminosity spec values')
plt.xlabel('Lyu E(B-V) values')
plt.ylabel('Luminosity spec values (Ergs/sec) in Log Scale')
plt.legend(['Lyu E(B-V)'],shadow=True, handlelength=1.5,fontsize=10)
plt.grid(linestyle='dashed')
plt.yscale('log')
ax = plt.axes()
ax.set_facecolor('black')
# plt.ylim(10,25)
plt.show()



#LUMINOSITY PHOT
# plt.plot(gordon,photrange,'.', color='red')
# plt.title('Fig 5: Reddening Gordon vs luminosity phot values')
# plt.xlabel('Gordon EB-V values')
# plt.ylabel('luminosity phot values (Ergs/sec) in Log Scale')
# plt.legend(('Gordon EB-V','luminosity phot'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.yscale('log')
# # plt.ylim(10,25)
# plt.show()
#
# plt.plot(lyu,photrange,'.', color='red')
# plt.title('Fig 4: Reddening Lyu vs luminosity phot values')
# plt.xlabel('Lyu EB-V values')
# plt.ylabel('luminosity phot values (Ergs/sec) in Log Scale')
# plt.legend(('Lyu EB-V','luminosity phot'),shadow=True, handlelength=1.5,fontsize=10)
# plt.grid(linestyle='dashed')
# plt.yscale('log')
# ax = plt.axes()
# ax.set_facecolor('black')
# # plt.ylim(10,25)
# plt.show()
