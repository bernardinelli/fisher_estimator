'''
Author: Pedro Bernardinelli (pedro.h.bernardinelli@gmail.com)
Date: Mar. 25, 2016

This short code gives an estimation of the survey dependent part of the Fisher matrix in a given supernovae Ia simulation in bins of redshift and
signal-to-noise ratio. 

Usage: call the analyse_and_plot(data_filte, plotformat, scaling) function. 
data_file is a file with the redshift z, the signal-to-noise ratio snr and the error in the distance modulus mu for each object. 
A complementary shell script file is provided for transforming SALT2mu's fitres files into the required format for this code. 
plotformat defines how the style of the lines of the plot. pylab/matplotlib accepts "--", ".-", "-", ":", ".", among others. 
Check http://matplotlib.org/1.3.1/examples/pylab_examples/line_styles.html for more information.
scaling is a parameter for scaling the amplitude of the Fisher matrix. 
A simple way to do use this parameter would be the number of objects per season divided by the number of simulated objects, or the purity and contamination
of the sample coming from photometric classification.

The function also returns dN/dz*\sigma_\mu^{-2} marginalized over the SNR distribution. This is interesting for comparing two 
simulations with very different ranges of redshifts (for example, SDSS and J-PAS). 
The marginalized_plot(marginalization, plotformat, curve_label) function plots this for the average z in each z bin. 


Warning: the values of each bin are only to be compared within other values for the same bin. Don't compare one bin to another. This happens because 
the full Fisher matrix has a cosmology dependent parameter that is not included here.

'''



import numpy as np
import pylab as pl

def redshiftbin(zmin,zmax, z, snr, sigma):
		i = 0
		zbin=[]
		snrbin=[]
		sigmabin=[]
		while i < len(z):
			if z[i] < zmax and z[i] > zmin:
				zbin.append(z[i])
				snrbin.append(snr[i])
				sigmabin.append(sigma[i])
			i = i + 1
		return zbin, snrbin, sigmabin

def SNRbin(SNRmin,SNRmax,SNR,sigma):
		i = 0
		sigmabin=[]
		while i < len(SNR):
			if SNR[i] < SNRmax and SNR[i] > SNRmin:
				sigmabin.append(sigma[i]**-2.)
			i = i + 1

		SNRav = (SNRmin+SNRmax)*0.5
		sig = np.average(sigmabin)
		A = len(sigmabin)*(sig)/0.05
		return A,SNRav

def findtheparameter(SNRmin,SNRmax,zmin,zmax,z,snr, sigma):
		zb,snrb,sigmab = redshiftbin(zmin,zmax,z, snr, sigma)
		A,AvSNR = SNRbin(SNRmin,SNRmax,snrb,sigmab)
		return A,AvSNR


def analyse_and_plot(data_file, plotformat, scaling):
	z,snr,sigma = np.loadtxt(data_file, unpack=True)

	marginalization = []
	SNRbinsize=1

	maxSNR = np.amax(snr)
	
	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.1,0.15,z,snr,sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(221)
	pl.ylabel(r"$\frac{\mathrm{d} N_z}{\mathrm{d} z} \sigma_\mu^{-2} $",fontsize=15)
	pl.plot(SNRdat,Adat,'b' + plotformat,label="$0.10<z<0.15$")
	marginalization.append(np.nansum(Adat))


	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.15,0.20,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(221)
	pl.ylabel(r"$\frac{\mathrm{d} N_z}{\mathrm{d} z} \sigma_\mu^{-2} $",fontsize=15)
	pl.plot(SNRdat,Adat,'r' + plotformat,label="$0.15<z<0.20$")
	marginalization.append(np.nansum(Adat))


	legend = pl.legend(loc=0, shadow=True,fontsize=15)


	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.20,0.25,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(222)
	pl.plot(SNRdat,Adat,'g' + plotformat,label='$0.20<z<0.25$')
	marginalization.append(np.nansum(Adat))



	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.25,0.30,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(222)
	pl.plot(SNRdat,Adat,'m' + plotformat,label='$0.25<z<0.30$')
	legend = pl.legend(loc=0, shadow=True,fontsize=15)
	marginalization.append(np.nansum(Adat))


	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.3,0.35,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(223)
	pl.xlabel("$SNR$",fontsize=15)
	pl.ylabel(r"$\frac{\mathrm{d} N_z}{\mathrm{d} z} \sigma_\mu^{-2} $",fontsize=15)
	pl.plot(SNRdat, Adat,'c' + plotformat,label='$0.30<z<0.35$')
	marginalization.append(np.nansum(Adat))



	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.35,0.40,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(223)
	pl.xlabel("$SNR$")
	pl.plot(SNRdat,Adat,'k' + plotformat,label="$0.35<z<0.40$")
	legend = pl.legend(loc=0, shadow=True,fontsize=15)
	marginalization.append(np.nansum(Adat))



	SNRi=0
	Adat=[]
	SNRdat=[]
	while SNRi<maxSNR:
		A,AvSNR = findtheparameter(SNRi,SNRi+1,0.4,0.45,z,snr, sigma)
		Adat.append(scaling*A)
		SNRdat.append(AvSNR)
		SNRi = SNRi + SNRbinsize
	pl.subplot(224)
	pl.plot(SNRdat,Adat,'m' + plotformat,label='$0.40<z<0.45$')
	pl.xlabel("$SNR$",fontsize=15)
	marginalization.append(np.nansum(Adat))

	return marginalization

'''
-------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

def marginalized_plot(marginalization, plotformat, curve_label):
	zbins = [0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425]

	pl.plot(zbins, marginalization, plotformat, label=curve_label)
	pl.xlabel("$z$",fontsize=15)
	pl.ylabel(r"$\sum _i\frac{\mathrm{d} N_z}{\mathrm{d} z} \sigma_\mu^{-2} SNR_i$",fontsize=15)

'''
-------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
jpas_mar = analyse_and_plot("jpas.log", "-", 227./10000.)
jpas_double_mar = analyse_and_plot("jpas.log", ":", 2*227./10000.)
jpas_3k_mar = analyse_and_plot("jpas.log", ":", 908./10000.)


sdss_mar = analyse_and_plot("SDSS.log", "--", 373./10000.)


pl.suptitle(r'Information', fontsize=15)
legend = pl.legend(loc=0, shadow=True)
frame = legend.get_frame()
fig = pl.gcf()
fig.set_size_inches(15,10)
pl.savefig("jpas.pdf")

pl.clf()

marginalized_plot(jpas_mar, '.-', curve_label="J-PAS 750 sq deg")
marginalized_plot(jpas_double_mar, '.-', curve_label="J-PAS 1500 sq deg")
marginalized_plot(jpas_3k_mar, '.-', curve_label="J-PAS 3000 sq deg")
marginalized_plot(sdss_mar, '.-', curve_label="SDSS")

pl.suptitle(r'Marginalized Information', fontsize=15)
legend = pl.legend(loc=0, shadow=True)
frame = legend.get_frame()
fig = pl.gcf()
fig.set_size_inches(15,10)
pl.savefig("marginal.pdf")
