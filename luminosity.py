from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from decimal import *
import decimal
#Load in the data
from astropy.table import Table
import pandas as pd
from numpy import arange
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.cosmology import WMAP9 as cosmo
import scipy.constants
import math

table = fits.open(r'C:\Users\anika\Downloads\catalog_file_updated2.fits')
tabledata = table[1].data
table.close()


rec=tabledata.REC_NO
mbh=tabledata.M_BH
kmass=tabledata.Kevin_Mass
kerr=tabledata.Kevin_Error
w1=tabledata.W1
w2=tabledata.W2
w3=tabledata.W3
w4=tabledata.W4
redshift=tabledata.SPEC_Z

rec0=[float(i) for i in rec]
mbh0=[float(i) for i in mbh]
kmass0=[float(i) for i in kmass]
kerr0=[float(i) for i in kerr]
w10=[float(i) for i in w1]
w20=[float(i) for i in w2]
w30=[float(i) for i in w3]
w40=[float(i) for i in w4]
redshift0=[float(i) for i in redshift]

kmass0=[10**i for i in kmass0]



mbhfinal=[]
for item in mbh0:
    if item>0:
        mbhfinal.append(item)

df=pd.DataFrame({'r':rec0, 'm':mbh0, 'w1': w10, 'w2': w20,'w3': w30,'w4': w40, 's':redshift0})

recrange=[]
massrange=[]
w1range=[]
w2range=[]
w3range=[]
w4range=[]
specrange=[]


for x in mbhfinal:
    recrange.append(df.loc[df['m'] == x, 'r'].iloc[0])
    massrange.append(df.loc[df['m'] == x, 'm'].iloc[0])
    w1range.append(df.loc[df['m'] == x, 'w1'].iloc[0])
    w2range.append(df.loc[df['m'] == x, 'w2'].iloc[0])
    w3range.append(df.loc[df['m'] == x, 'w3'].iloc[0])
    w4range.append(df.loc[df['m'] == x, 'w4'].iloc[0])
    specrange.append(df.loc[df['m'] == x, 's'].iloc[0])

value=[]
for i in range(len(mbhfinal)):
    fr=recrange[i]
    fm=massrange[i]
    fw1=w1range[i]
    fw2=w2range[i]
    fw3=w3range[i]
    fw4=w4range[i]
    fs=specrange[i]
    #
    # # de redshift WISE values
    # fw1=fw1/(1+fs)
    # fw2=fw2/(1+fs)
    # fw3=fw3/(1+fs)
    # fw4=fw4/(1+fs)

    # Convert the Wise values from AB magnitude to fluxes. axis is ergs/s/Hz/cm^2
    w1f=3631/(100**(fw1/5)) * 10**(-23)
    w2f=3631/(100**(fw2/5)) * 10**(-23)
    w3f=3631/(100**(fw3/5)) * 10**(-23)
    w4f=3631/(100**(fw4/5)) * 10**(-23)

    # convert WISE values to luminosity but first get distance
    d=cosmo.luminosity_distance(fs).value
    d1=d*3.086*(10**(24))
    # print(type(d1))
    p=np.pi

    # Now get lum values for all WISE VALUES
    lum1=4*p*(d1**2)*(w1f/(1+fs))
    lum2=4*p*(d1**2)*(w2f/(1+fs))
    lum3=4*p*(d1**2)*(w2f/(1+fs))
    lum4=4*p*(d1**2)*(w4f/(1+fs))
    #
    # print("redshift is = " + str(fs))
    # print("distance d is = " + str(d1))
    #
    # print("lum wise = " + str(lum1))
    # print("lum wise = " + str(lum2))
    # print("lum wise = " + str(lum3))
    # print("lum wise = " + str(lum4))

    # Get the template we will fit the WISE values
    template = fits.open(r'C:\Users\anika\Downloads\featureless_AGN_SED.fits')
    templatedata = template[1].data
    template.close()

    # Convert the lnu to smaller values by dividing by 10**24
    Wavelength=templatedata.Wavelength
    # wavelength=[i for i in Wavelength if i <600]

    Lnu=templatedata.Lnu
    lnu=Lnu/10**24
    # lnu=lnu[:5888]


    wise=np.array([[3.4, lum1], [4.6, lum2], [12, lum3], [22, lum4]])
    list=[3.4,4.6,12,22]
    listf=[lum1,lum2,lum3,lum4]

    # convert list wavelengths to nu and c is in microns
    list2=[(scipy.constants.c *(10**6))/i for i in list]
    l6x=(scipy.constants.c *(10**6))/6
    # interpolate the wise values to template

    ff=interpolate.interp1d(list,listf, bounds_error=False)
    ynew=ff(Wavelength)

    Wavelength1=[(scipy.constants.c*(10**6))/i for i in Wavelength]
     #c is in 10**8 and it need to be in 10**14
    # integrate the lnu vakues for given wavelength VALUES


    wave=[i for i in Wavelength if (i >5.7 and i<6.3)]
    df=pd.DataFrame({'w':Wavelength,'l':ynew})
    lnurange=[]
    for x in wave:
        lnurange.append(df.loc[df['w'] == x, 'l'].iloc[0])

    # print("lnurange = " + str(lnurange))
    # print("wave = " + str(wave))
    df=pd.DataFrame({'w':Wavelength,'l':Wavelength1})
    nu=[]
    for x in wave:
        nu.append(df.loc[df['w'] == x, 'l'].iloc[0])
    nu2=nu[::-1]
    # print("nu values = " + str(nu))
    import scipy.integrate as it
    integral = it.cumtrapz(lnurange, x=nu2)

    value.append(integral[0]) #THIS IS IN UNITS OF ERG/S/HZ#

#
# print(value)
# f=[np.abs(x) for x in value]
L6 = fits.BinTableHDU.from_columns(
        [fits.Column(name='REC_NO', format='J', array=recrange),
        fits.Column(name='L6', format='D', array=value)])

L6.writeto(r'C:\Users\anika\Downloads\L6table2.fits')

# plt.plot(Wavelength,ynew)
# plt.plot(list[0],w1f, 'ro', color='r')
# plt.plot(list[1],w2f, 'ro', color='b')
# plt.plot(6,value, 'ro', color='k')
# plt.plot(list[2],w3f, 'ro', color='g')
# plt.plot(list[3],w4f, 'ro', color='y')
# plt.xlabel('Wavelength (Microns)')
# plt.ylabel('Flux (W/Hz)')
# plt.yscale("log")
# plt.xscale("log")
# plt.show()
