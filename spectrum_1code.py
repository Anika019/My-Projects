from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from decimal import *

#Load in the AGN obserbed spectrum
from astropy.table import Table
AGN = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\Hershel spectra\sp_1\spectrum_1fit.fits')
print(AGN)
AGNdata = AGN[1].data
AGN.close()

#set redshift  from the catalog
z=0.989

#set up our variables
AGNwavelength=AGNdata.Wavelength/(1+z)   #I change the wavelength to rest-frame by dividing by 1+z  Reminder:  This is in angstroms
AGNflux=AGNdata.Flux




#Load in the reddening template
template = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\kext_albedo_WD_SMCbar_0.fits')
templatedata = template[1].data
template.close()

draineind0=np.where((templatedata.micron > 0.05) & (templatedata.micron < 2))
draineind=draineind0[0]

originaltemplatewavelength=templatedata.micron[draineind]*10**4 #wavelength converted to angstroms
originaltemplatekappa=templatedata.kappa[draineind]*0.0001

#Now we want the template to have the same number of data points, and spaced equally as, your data, so let's use your interp1d code

ff=interpolate.interp1d(originaltemplatewavelength,originaltemplatekappa, bounds_error=False)
templatekappa=ff(AGNwavelength)
#now we have a set of kappa values for each of our wavelength points of the spectrum,


#Now let's load in the Gordon template that provides the f_0 values for the function
gordon = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\gordon06_templates.fits')
gordondata = gordon[1].data
gordon.close()

#Since Gordon uses log(Hz), we need to change frequency to wavelength in angstroms
gordonwavelength=(2.998*10**18)/(10**gordondata.hz) #now we have wavelength in angstroms
gordonwavelength=gordonwavelength

#Now we need f_0 values also in terms of per wavelength instead of per Hz, with an additional division by 10^44 to make it smaller
gordonf0=(10**(gordondata.qsolum-44))/gordonwavelength

#also we need to make sure f0 is interpolated to the same wavelength scale as we used earlier...
ffgordon=interpolate.interp1d(gordonwavelength,gordonf0, bounds_error=False)
f0=ffgordon(AGNwavelength)

#Now let's build the reddened spectrum function
#This function is modifying y input (the spectrum flux) and a set of parameters that will be fit
#the kappas from the template are fixed, so I'm just using them and not fitting them here.
def func(X,norm, ebv):
    f0,kappa = X
    f_lambda=norm*f0*np.exp(-(kappa*ebv)/1.086)
    return(f_lambda)

#I personally use curvefit to do fitting, but there's probably a similar call for lmfit like you've been using.
popt, pcov = curve_fit(func,(f0,templatekappa), AGNflux)#,bounds=[[0,0],[20,10]])
print(popt)
#popt gives us the fit parameters, which in this case is the normalization which controls the vertical placement, and ebv, which controls the curvature
x=popt[0]
y=popt[1]
textstr= '\n'.join(('Fit Parameters Gordon:',
                    'Normalization Gordon:'+ str(x),
                    'EBV Gordon:'+ str(y)))

fig = plt.figure()
plt.step(AGNwavelength, AGNflux,linewidth=0.5,color='blue')
plt.step(AGNwavelength, func((f0,templatekappa), *popt),linewidth=2,color='red')
plt.title('Fig: Reddening code for spectrum 1 with Gordon template')
plt.xlabel('AGN Wavelength (Angstrom)')
plt.ylabel('AGN Flux (10^(-17) ergs/s/cm^2/Å))')
# plt.text(4000, 15, 'Fit parameters:normalizaton:3.36517990e+02, ebv:-2.75477242e-01', fontsize=8)
plt.text(3000, 15, textstr, fontsize=8,
        verticalalignment='top',bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),))
plt.legend(('Spectrum','Reddening Curve'),shadow=True, handlelength=1.5,fontsize=10)
plt.show()





#Now let's load in the Lyu template that provides the f_0 values for the function
lyu = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\FALL 2020 RESEARCH\Lyu_qso_luminosity1.fits')
lyudata = lyu[1].data
lyu.close()

#Since Lyu uses log(Hz), we need to change frequency to wavelength in angstroms
lyuwavelength=(2.998*10**18)/(10**lyudata.frequency) #now we have wavelength in angstroms
lyuwavelength=lyuwavelength

#Now we need f_0 values also in terms of per wavelength instead of per Hz, with an additional division by 10^44 to make it smaller
#lyuf0=(10**(lyudata.qsolum-11))/lyuwavelength
lyuf0=lyudata.qsolum/lyuwavelength


#also we need to make sure f0 is interpolated to the same wavelength scale as we used earlier...
fflyu=interpolate.interp1d(lyuwavelength,lyuf0, bounds_error=False)
f1=fflyu(AGNwavelength)



#Now let's build the reddened spectrum function
#This function is modifying y input (the spectrum flux) and a set of parameters that will be fit
#the kappas from the template are fixed, so I'm just using them and not fitting them here.
def func1(X,norm, ebv):
    f1,kappa = X
    f1_lambda=norm*f1*np.exp(-(kappa*ebv)/1.086)
    return(f1_lambda)

#I personally use curvefit to do fitting, but there's probably a similar call for lmfit like you've been using.
popt1, pcov1 = curve_fit(func1,(f1,templatekappa), AGNflux)#,bounds=[[0,0],[20,10]])
print(popt1)
#popt gives us the fit parameters, which in this case is the normalization which controls the vertical placement, and ebv, which controls the curvature
x1=popt1[0]
y1=popt1[1]
textstr1= '\n'.join(('Fit Parameters Lyu:',
                    'Normalization Lyu:'+ str(x1),
                    'EBV Lyu:'+ str(y1)))

fig1 = plt.figure()
plt.step(AGNwavelength, AGNflux,linewidth=0.5,color='blue')
plt.step(AGNwavelength, func1((f1,templatekappa), *popt1),linewidth=2,color='red')
plt.title('Fig: Reddening code for spectrum 1 Lyu Template')
plt.xlabel('AGN Wavelength (Angstrom)')
plt.ylabel('AGN Flux (10^(-17) ergs/s/cm^2/Å))')
# plt.text(4000, 15, 'Fit parameters:normalizaton:3.36517990e+02, ebv:-2.75477242e-01', fontsize=8)
plt.text(3000, 15, textstr1, fontsize=8,
        verticalalignment='top',bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),))
plt.legend(('Spectrum','Reddening Curve'),shadow=True, handlelength=1.5,fontsize=10)
plt.show()





import pandas as pd

AGN_w0=[float(x) for x in AGNwavelength]
AGN_f0=[float(x) for x in AGNflux]
df=pd.DataFrame({'agn_w':AGN_w0,'agn_f':AGN_f0})
AGN_wl = [i for i in AGN_w0 if i > 2500 and i<4500]
AGN_fx=[]
for x in AGN_wl:
    AGN_fx.append(df.loc[df['agn_w'] == x, 'agn_f'].iloc[0])
template_kappa=ff(AGN_wl)

#Now let's load in the Glikman qso template that provides the f_1 values for the function
glik = fits.open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\SPRING 2021 documents\Glikman template.fits')
glikdata = glik[1].data
glik.close()

#Since glik uses log(Hz), we need to change frequency to wavelength in angstroms
#glikwavelength=(2.998*10**18)/(glikdata.frequency) #now we have wavelength in angstroms
glikwavelength=glikdata.wavelength

#Now we need f_0 values also in terms of per wavelength instead of per Hz, with an additional division by 10^44 to make it smaller
#lyuf0=(10**(lyudata.qsolum-11))/lyuwavelength
glikf0=(glikdata.qsolum*10**14)/glikwavelength


#also we need to make sure f0 is interpolated to the same wavelength scale as we used earlier...
ffglik=interpolate.interp1d(glikwavelength,glikf0, bounds_error=False)
f2=ffglik(AGN_wl)


#Now let's build the reddened spectrum function
#This function is modifying y input (the spectrum flux) and a set of parameters that will be fit
#the kappas from the template are fixed, so I'm just using them and not fitting them here.
def func2(X,norm, ebv):
    f2,kappa = X
    f2_lambda=norm*f2*np.exp(-(kappa*ebv)/1.086)
    return(f2_lambda)

#I personally use curvefit to do fitting, but there's probably a similar call for lmfit like you've been using.
popt2, pcov2 = curve_fit(func2,(f2,template_kappa), AGN_fx)#,bounds=[[0,0],[20,10]])
print(popt2)
#popt gives us the fit parameters, which in this case is the normalization which controls the vertical placement, and ebv, which controls the curvature
x2=popt2[0]
y2=popt2[1]
textstr2= ('E(B-V) ='+ str(round(y2,3)))




fig2 = plt.figure()
plt.step(AGN_wl, AGN_fx,linewidth=0.5,color='blue')
plt.step(AGN_wl, func2((f2,template_kappa), *popt2),linewidth=2,color='red')
plt.title("Spectrum 1")
plt.xlabel('Rest Wavelength (Angstrom)')
plt.ylabel('AGN Flux (10^(-17) ergs/s/cm^2/Å))')
# plt.text(4000, 15, 'Fit parameters:normalizaton:3.36517990e+02, ebv:-2.75477242e-01', fontsize=8)
# plt.text(textstr2, fontsize=8,
#         verticalalignment='top',bbox=dict(boxstyle="round",
#                    ec=(1., 0.5, 0.5),
#                    fc=(1., 0.8, 0.8),))
plt.legend([textstr2],ncol=1,handlelength=0,fontsize=10)
# plt.ylim(1,4)
#plt.legend(('Spectrum','Reddening Curve'),shadow=True, handlelength=1.5,fontsize=10)
plt.show()

#

file2write=open(r'C:\Users\anika\OneDrive - The University of Kansas\Desktop\ASTRO\SPRING 2021 documents\Glikman EBV.txt', "a")
file2write.writelines(["Spectrum 1","\t", str(round(y2,3)),"\n"])
file2write.close()
