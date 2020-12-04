# Microwave coaxial cable loss calculator
# Adrian Sinclair
# Modified and updated with thermal calculations by Jeremy Meinke
# References:
#   coax co,ltd SC-086/50-SS-SS
#   http://www.coax.co.jp/en/product/sc/086-50-ss-ss.html
#   http://www.gb.nrao.edu/electronics/edir/edir223.pdf
#   http://www.coax.co.jp/en/product/sc/219-50-ss-ss.html
#   http://www.coax.co.jp/en/product/sc/086-50-nbti-nbti.html
#   https://rfcoax.com/technical-drawing/S086MMHF-20.5.html
#   https://www.tek-stock.com/ut-085-ss-ss/

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# plotting = False
# input parameters for testing
#ctype = "SC-219/50-SS-SS"
#f = 0.500 # GHz
#length = 1000.0 # mm
#Temp = 300. # K

####setting up thermal calculations and functions
def func_nist(x, a, b, c, d, e, f, g, h, i):
	logT = np.log10(x)
	return pow(10,a + b*logT + c*pow(logT,2) + d*pow(logT,3) + e*pow(logT,4) + f*pow(logT,5) + g*pow(logT,6) + h*pow(logT,7) + i*pow(logT,8))
def func_nist_2(x, a, b, c, d, e, f, g, h, i):
	return pow(10, (a + c * pow(x, 0.5) + e * x + g * pow(x, 1.5) + i * x * x)/(1 + b * pow(x, 0.5) + d * x + f * pow(x, 1.5) + h * x * x))
def func_exp_simple(x,a,b):
	return (a * pow(x,b))

def func_nist_304_SS_thermal_conductivity(T):   ##Same as 316 SS
	return np.piecewise(T/1,[T<1,T>=1],[lambda x: pow(x,1.2)*func_nist(1,-1.4087,1.3982,0.2543,-0.626,0.2334,0.4256,-0.4658,0.1650,-0.0199),lambda x: func_nist(x,-1.4087,1.3982,0.2543,-0.626,0.2334,0.4256,-0.4658,0.1650,-0.0199)])
def func_NbTi_thermal_conductivity(T):  ###From (J.R. Olson) Thermal conductivity of some common cryostat materials between 0.05 and 2 K, likely good to ~4K....
	return func_exp_simple(T,0.015,2)
def func_nist_PTFE_thermal_conductivity(T):  ##AKA teflon, fixed so anything below 4K fit to 1.88
	return np.piecewise(T/1,[T<4,T>=4],[lambda x: pow(x/4,1.88)*func_nist(4, 2.738, -30.677, 89.43, -136.99, 124.69, -69.556, 23.32, -4.3135, 0.33829 ),lambda x: func_nist(x, 2.738, -30.677, 89.43, -136.99, 124.69, -69.556, 23.32, -4.3135, 0.33829 )])
# def func_RRR_50_copper_thermal_conductivity(T):
# 	return func_nist_2(T, 1.8743, -0.41538, -0.6018, 0.13294, 0.26426, -0.0219, -0.051276, 0.0014871, 0.003723 )
# def func_RRR_100_copper_thermal_conductivity(T):
# 	return func_nist_2(T, 2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281, 0.003207)
# def func_RRR_150_copper_thermal_conductivity(T):
# 	return func_nist_2(T, 2.3797, -0.4918, -0.98615, 0.13942, 0.30475, -0.019713, -0.046897, 0.0011969, 0.0029988)
# def func_RRR_300_copper_thermal_conductivity(T):
# 	return func_nist_2(T, 1.357, 0.3981, 2.669, -0.1346, -0.6683, 0.01342, 0.05773, 0.0002147, 0)
# def func_RRR_500_copper_thermal_conductivity(T):
# 	return func_nist_2(T, 2.8075, -0.54074, -1.2777, 0.15362, 0.36444, -0.02105, -0.051727, 0.0012226, 0.0030964)


# Cable Parameters
# f - frequency data points in GHz
# LT1, LT2 - loss in dB/m at two different temperatures
# T1, T2 - Temperatures at which loss was measured
def cableParams(ctype):
	global f
	global LT1
	global LT2
	global T1
	global T2
	if ctype == "SC-086/50-SS-SS":
		# SC-086/50-SS-SS
		# http://www.coax.co.jp/en/product/sc/086-50-ss-ss.html  
		f = np.array([0.5,1,5,10,20]) # GHz
		LT1 = np.array([7.3,10.3,23.,32.7,46.4]) # dB/m
		LT2 = np.array([4.7,6.6,14.8,20.9,29.5]) # dB/m
		T1 = 300 # K
		T2 = 4 # K
	elif ctype == "SC-086/50-CN-CN":
		# SC-086/50-CN-CN
		# http://www.coax.co.jp/en/product/sc/086-50-cn-cn.html  
		f = np.array([0.5,1,5,10,20]) # GHz
		LT1 = np.array([5.4,7.7,17.1,24.3,34.6]) # dB/m
		LT2 = np.array([4.1,5.7,12.8,18.1,25.7]) # dB/m
		T1 = 300 # K
		T2 = 4 # K
	elif ctype == "SC-219/50-SS-SS":
		# SC-219/50-SS-SS
		# http://www.coax.co.jp/en/product/sc/219-50-ss-ss.html
		f = np.array([0.5,1,5,10,20]) # GHz
		LT1 = np.array([3.0,4.2,9.4,13.5,19.2]) # dB/m
		LT2 = np.array([1.9,2.6,5.9,8.3,11.7]) # dB/m
		T1 = 300. # K
		T2 = 4. # K
	elif ctype == "SC-219/50-CN-CN":
		# SC-219/50-CN-CN
		# http://www.coax.co.jp/en/product/sc/219-50-cn-cn.html
		f = np.array([0.5,1,5,10,20]) # GHz
		LT1 = np.array([2.4,4.3,7.6,10.8,15.5]) # dB/m
		LT2 = np.array([1.6,2.3,5.1,7.2,10.2]) # dB/m
		T1 = 300. # K
		T2 = 4. # K
	elif ctype == "UT-085-SS-SS":
		# UT-085-SS-SS
		# https://www.tek-stock.com/ut-085-ss-ss/
		f = np.array([0.5,1,5,10,18]) # GHz
		LT1 = np.array([2.92,4.13,9.32,13.25,18.90])  # dB/m
		LT2 = np.array([2.92,4.12,9.31,13.24,18.89])  # dB/m
		T1 = 300. # K
		T2 = 299. # K
	elif ctype == "SC-086/50-NbTi-NbTi":  
		# SC-086/50-NbTi-NbTi
		# http://www.coax.co.jp/en/product/sc/086-50-nbti-nbti.html
		f = np.array([0.5,1,5,10,20]) # GHz
		LT1 = np.array([6.8,9.6,21.6,30.5,43.1]) # dB/m
		LT2 = np.array([6.7,9.5,21.5,30.4,43.0]) # dB/m
		T1 = 300. # K
		T2 = 299. # K
	elif ctype == "S086MMHF":  
		# S086MMHF
		# https://rfcoax.com/technical-drawing/S086MMHF-20.5.html
		f = np.array([18.,27.,40.,65.]) # GHz
		LT1 = np.array([3.6,4.6,6.6,7.9]) # dB/m
		LT2 = np.array([3.59,4.59,6.59,7.89]) # dB/m
		T1 = 300. # K
		T2 = 299. # K
	return f, LT1, LT2, T1, T2

def cableCond(ctype,Tin,Tout,L):
	Tmin=min(Tin,Tout)
	Tmax=max(Tin,Tout)
	if ctype == "SC-086/50-SS-SS":
		# SC-086/50-SS-SS
		# http://www.coax.co.jp/en/product/sc/086-50-ss-ss.html  
		# CryoCoax has nearly identical
		Cond=True
		OD=.00086
		DD=.00066
		ID=.000203  ##All in [m]
		conductor_x_sect=(OD**2-DD**2+ID**2)*np.pi/4   ##m
		PTFE_x_sect=(DD**2-ID**2)*np.pi/4  ##m
		def cond_func(T):
			return conductor_x_sect*func_nist_304_SS_thermal_conductivity(T)+PTFE_x_sect*func_nist_PTFE_thermal_conductivity(T)
	elif ctype == "SC-086/50-CN-CN":
		# SC-086/50-CN-CN
		# http://www.coax.co.jp/en/product/sc/086-50-cn-cn.html  
		Cond=False
	elif ctype == "SC-219/50-SS-SS":
		# SC-219/50-SS-SS
		# http://www.coax.co.jp/en/product/sc/219-50-ss-ss.html
		Cond=True
		OD=.00219
		DD=.00167
		ID=.00051
		conductor_x_sect=(OD**2-DD**2+ID**2)*np.pi/4   ##m
		PTFE_x_sect=(DD**2-ID**2)*np.pi/4  ##m
		def cond_func(T):
			return conductor_x_sect*func_nist_304_SS_thermal_conductivity(T)+PTFE_x_sect*func_nist_PTFE_thermal_conductivity(T)
	elif ctype == "SC-219/50-CN-CN":
		# SC-219/50-CN-CN
		# http://www.coax.co.jp/en/product/sc/219-50-cn-cn.html
		Cond=False
	elif ctype == "UT-085-SS-SS":
		# UT-085-SS-SS
		# https://www.tek-stock.com/ut-085-ss-ss/
		Cond=True	
	elif ctype == "SC-086/50-NbTi-NbTi":  
		# SC-086/50-NbTi-NbTi
		# http://www.coax.co.jp/en/product/sc/086-50-nbti-nbti.html
		print('NbTi thermal values only valid below 4K')
		Cond=True
		OD=.00090
		DD=.00066
		ID=.000203
		conductor_x_sect=(OD**2-DD**2+ID**2)*np.pi/4   ##m
		PTFE_x_sect=(DD**2-ID**2)*np.pi/4  ##m
		def cond_func(T):
			return conductor_x_sect*func_NbTi_thermal_conductivity(T)+PTFE_x_sect*func_nist_PTFE_thermal_conductivity(T)
	elif ctype == "S086MMHF":  
		# S086MMHF
		# https://rfcoax.com/technical-drawing/S086MMHF-20.5.html
		Cond=False
	if Cond:
		ratio=0
		tot=integrate.quad(cond_func,Tmin,Tmax)[0]
		Tavg=Tmin
		while ratio<0.5:
			Tavg+=(Tmax-Tmin)/1000
			ratio=integrate.quad(cond_func,Tmin,Tavg)[0]/tot
		return Tavg, tot/(L/1000)		###as L is in mm
	else:
		return  (Tmax+Tmin)/2,'Undefined'			###just to keep shape same



def lsFit(f,L):
	logf = np.log(f)
	logL = np.log(L)
	m,b = np.polyfit(logf,logL,1)
	Afit = np.exp(b)
	return Afit

def LT1fit(fin):
	return Afit_T1*np.sqrt(fin)

def LT2fit(fin):
	return Afit_T2*np.sqrt(fin)

def LperMeter(T,f): # frequency in GHz, Temperature in K
	M = ((LT1fit(f)-LT2fit(f)) / (T1-T2))
	b = LT1fit(f) - LT1fit(f)*T1/(T1-T2) + LT2fit(f)*T1/(T1-T2) 
	Lm = M*T + b # loss in [dB/m]
	return Lm

def L(Lm,length):
	L = length*Lm/1.0e3 # loss in [dB]
	return L
 
def getLoss(f_0,Temp,length,ctype):
	global Afit_T1
	global Afit_T2
	f_0 = f_0/1.e3 # convert from MHz to GHz
	if (ctype == "SC-086/50-NbTi-NbTi") and (Temp < 10.0):
		LtoReturn = 0.5*length/1.e3 # superconducting state loss max of 0.5 dB/m
	else:
		f , LT1, LT2, T1, T2 = cableParams(ctype)
		Afit_T1 = lsFit(f,LT1)
		Afit_T2 = lsFit(f,LT2)
		Lm = LperMeter(Temp,f_0)
		LtoReturn = L(Lm,length)
	#print "Loss at freq "+str(f_0)+"GHz at temp "+str(Temp)+"K with cable "+ctype+" :"+str(LtoReturn)+"[dB]"
	return LtoReturn # Loss in [dB]

# if plotting == True:
# 	# for plotting above fitting equations
# 	f_fit = np.linspace(0,20,2001)
# 	f , LT1, LT2, T1, T2 = cableParams(ctype)
# 	Afit_T1 = lsFit(f,LT1)
# 	Afit_T2 = lsFit(f,LT2)
# 	L_T1_fit = Afit_T1*np.sqrt(f_fit)
# 	L_T2_fit = Afit_T2*np.sqrt(f_fit)
# 	Lfit = np.zeros(len(f_fit))
# 	for i in range(len(f_fit)):
# 		Lfit[i] = LperMeter(Temp,f_fit[i])
# 	plt.ion()
# 	plt.scatter(f,LT1,label="meas 300K",c='red')
# 	plt.scatter(f,LT2,label="meas 4K",c='orange')
# 	plt.plot(f_fit,L_T1_fit,label="fit 300K")
# 	plt.plot(f_fit,Lfit,label="Lfit "+str(Temp)+"K",c='black')
# 	plt.plot(f_fit,L_T2_fit,label="fit 4K")
# 	plt.ylabel("Loss [dB/m]",size=14)
# 	plt.xlabel("Frequency [GHz]",size=14)
# 	plt.title("Frequency dependent loss of "+ctype)
# 	plt.legend()
# 	plt.tight_layout()
# 	plt.show()
