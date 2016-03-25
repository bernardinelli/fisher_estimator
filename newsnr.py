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
				sigmabin.append(sigma[i])
			i = i + 1

		SNRav = (SNRmin+SNRmax)*0.5
		sig = np.average(sigmabin)
		A = len(sigmabin)*(sig**-2)/0.05
		return A,SNRav

def findtheparameter(SNRmin,SNRmax,zmin,zmax,z,snr, sigma):
		zb,snrb,sigmab = redshiftbin(zmin,zmax,z, snr, sigma)
		A,AvSNR = SNRbin(SNRmin,SNRmax,snrb,sigmab)
		return A,AvSNR


def analyse_and_plot(data_file, plotformat, scaling):
	z,snr,sigma = np.loadtxt(data_file, unpack=True)

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

'''
-------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

analyse_and_plot("jplus.log", "--", 1)

analyse_and_plot("jpas.log", "-", 2)

pl.suptitle(r'Information', fontsize=15)
legend = pl.legend(loc=0, shadow=True)
frame = legend.get_frame()
fig = pl.gcf()
fig.set_size_inches(15,10)
pl.savefig("jpas.pdf")