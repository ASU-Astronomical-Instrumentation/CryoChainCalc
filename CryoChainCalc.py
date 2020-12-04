# Cryogenic chain calculator for cascaded microwave systems
# Adrian Sinclair
# Updated and modified by Jeremy Meinke

import numpy as np
import SchemDraw as schem
import customElements as e
import matplotlib.pyplot as plt
from cable import *
# manual for SchemDraw
# https://cdelker.bitbucket.io/SchemDraw/SchemDraw.html

# initialize
#S = -10 # dBm input signal power
d = schem.Drawing() # initialize schematic
nParams = 5 # number of parameters for each component stored in array p
p = np.zeros(nParams) # initialize component parameter array p
plt.ion()
##########################
# Amplifier function
##########################
def addAmp(T,G,S=0.,p=0):
		"""
		Params: 
		T : noise temperature of the attenuator [K]
		G: gain value in [dB]
		S: signal power in [dBm]
		p : parameter array for components
		"""
		
		NSD = 0
		Tf = 0
		Ti = 0
		P = 0
		Gprod = 0
		
		if np.size(np.shape(p)) == 0:
			Tf = T
			# output signal power
			P = S+G # [dBm]
			NSD =10.0*np.log10( 1.38e-23*T*1000.0) # noise spectral density [dBm/Hz]
			p = np.array([S, T, G, P, NSD])
		
		else:
			#Noise calculation
			p = np.vstack((p,np.array([S,T,G,P,0])))      
			# output signal power
			S,p[-1][0] = p[-2][3], p[-2][3]
			P = p[-1][0]+G # [dBm]
			p[-1][3] = P  # assign output power to new array      
			Gc = 10.0**(p[:,2]/10.0) # component gain array and convert to linear
			p_flat = p.flatten()
			for i in range(int(len(p_flat)/nParams)): # Friis cascade noise loop
					Ti = p[i][1] # grab specific components noise, or thermal temp
					if i == 0:
						Tf += Ti
					else:
						Gprod = np.prod(Gc[:i])
						Tf += Ti / Gprod
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			# add NSD for component
			p[-1][-1] = NSD    
			Gprod = np.prod(Gc[:i+1])
		# draw component and label with params S, T, G, P, NSD
		a1 = d.add(e.AMP,d = 'right')
		a1.add_label("Amp",loc='top')
		a1.add_label("S [dBm]: "+str(round(S,4)),loc='bot')
		a1.add_label("T [K]: "+str(round(T,4)),loc='bot',ofst=1.0)
		a1.add_label("G [dB]: "+str(round(G,4)),loc='bot',ofst=2.0)
		a1.add_label("P [dBm]: "+str(round(P,4)),loc='bot',ofst=3.0)
		a1.add_label("Tcas [K]: "+str(round(Tf,4)),loc='bot',ofst=4.0)
		a1.add_label("T_N [K]: "+str(round(Tf*Gprod,4)),loc='bot',ofst=5.0)
		#a1.add_label("NSD [dBm/Hz]: "+str(round(NSD,2)),loc='bot',ofst=5.0)
		d.add(e.LINE, d='right', l=5)
		return p

##########################
# Attenuator function
##########################
def addAtten(T,A,ntones=1000,S=0.,p=0):
		"""
		Params:
		T : physical temperature of the attenuator [K]
		A: attenuation value in [dB]
		S: signal power in [dBm]
		p : parameter array for components
		"""
		NSD = 0
		Tf = 0
		Ti = 0
		P = 0
		Gprod = 0
		
		if np.size(np.shape(p)) == 0:
			Tf = T*(10**(abs(A)/10.) - 1.) # effective noise temp of atten
			# output signal power
			P = S-abs(A) # [dBm]
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			p = np.array([S, T*(10.**(abs(A)/10.0)-1.), -abs(A), P, NSD])
		
		else:
			#Noise calculation with temperature modified for attenuator T_equiv = T_phys*(Atten-1)
			p = np.vstack((p,np.array([S,T*(10.**(abs(A)/10.0)-1.),-abs(A),P,0])))
			# output signal power
			S,p[-1][0] = p[-2][3], p[-2][3]
			P = p[-1][0]-abs(A) # [dBm]
			p[-1][3] = P  # assign output power to new array      
			Gc = 10.0**(p[:,2]/10.0) # component gain array and convert to linear
			p_flat = p.flatten()
			for i in range(int(len(p_flat)/nParams)): # Friis cascade noise loop
					Ti = p[i][1] # grab specific components noise, or thermal temp
					if i == 0:
						Tf += Ti
					else:
						Gprod = np.prod(Gc[:i])
						Tf += Ti / Gprod
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			# add NSD for component
			p[-1][4] = NSD
			Gprod = np.prod(Gc[:i+1])
		# draw component and label with params S, T, G, P, NSD
		a1 = d.add(e.RES,d = 'right')
		a1.add_label("Atten",loc='top')
		a1.add_label("S [dBm]: "+str(round(S,4)),loc='bot',ofst=0.7)
		a1.add_label("T [K]: "+str(round(T,4)),loc='bot',ofst=1.7)
		a1.add_label("A [dB]: "+str(round(-abs(A),4)),loc='bot',ofst=2.7)
		a1.add_label("P [dBm]: "+str(round(P,4)),loc='bot',ofst=3.7)
		a1.add_label("Tcas [K]: "+str(round(Tf,4)),loc='bot',ofst=4.7)
		a1.add_label("T_N [K]: "+str(round(Tf*Gprod,4)),loc='bot',ofst=5.7)
		Pd=1e-3*10.0**(S/10.0)*(1-10**(-A/10))*ntones
		if Pd>.1:
			a1.add_label("%i tone P_diss [W]: "%ntones+str(round(Pd,2)),loc='bot',ofst=6.7)
		elif Pd>.001:
			a1.add_label("%i tone P_diss [mW]: "%ntones+str(round(Pd*1e3,2)),loc='bot',ofst=6.7)
		else:
			a1.add_label("%i tone P_diss [uW]: "%ntones+str(round(Pd*1e6,2)),loc='bot',ofst=6.7)
		
		#a1.add_label("NSD [dBm/Hz]: "+str(round(NSD,2)),loc='bot',ofst=5.0)
		d.add(e.LINE, d='right', l=5)
		return p

##########################
# Cable function
##########################
def addCable(Tin,Tout,L,fin, ntones=1000, ctype="SC-086/50-SS-SS",S=0.0,p=0):
		"""
		Params:
		Tin : temperature of stage of input [K]
		Tout : temperature of stage of output [K]
		T : physical temperature of the attenuation in [K]
		L : length of cable segment in [mm]
		fin : center frequency of input signal in [MHz]
		ctype: cable type as a string, see cable.py
		S: signal power in [dBm]
		p : parameter array for components
		"""
		T,therm=cableCond(ctype,Tin,Tout,L)
		NSD = 0
		Tf = 0
		Ti = 0
		P = 0
		Gprod = 0
		#Call cable temperature and frequency attenuation function from cable.py
		A = getLoss(fin,T,L,ctype)
		if np.size(np.shape(p)) == 0:
			Tf = T*(10**(abs(A)/10.) - 1.)
			# output signal power
			P = S-abs(A) # [dBm]
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			p = np.array([S, T*(10.**(abs(A)/10.0)-1.), -abs(A), P, NSD])
		
		else:
			#Noise calculation with temperature modified for attenuator T_equiv = T_phys*(Atten-1)
			p = np.vstack((p,np.array([S,T*(10.**(abs(A)/10.0)-1.),-abs(A),P,0])))
			# output signal power
			S,p[-1][0] = p[-2][3], p[-2][3]
			P = p[-1][0]-abs(A) # [dBm]
			p[-1][3] = P  # assign output power to new array      
			Gc = 10.0**(p[:,2]/10.0) # component gain array and convert to linear
			#print(Gc)
			p_flat = p.flatten()
			for i in range(int(len(p_flat)/nParams)): # Friis cascade noise loop
					Ti = p[i][1] # grab specific components noise, or thermal temp
					if i == 0:
						Tf += Ti
					else:
						Gprod = np.prod(Gc[:i])
						Tf += Ti / Gprod
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			# add NSD for component
			p[-1][4] = NSD
			Gprod = np.prod(Gc[:i+1])
		# draw component and label with params S, T, G, P, Tcas, NSD
		a1 = d.add(e.CABLE, d='right')
		a1.add_label(ctype+" Cable",loc='top')
		a1.add_label("S [dBm]: "+str(round(S,4)),loc='bot')
		a1.add_label("(%.2fK-%.2fK), Tavg [K]: "%(Tin,Tout)+str(round(T,2)),loc='bot',ofst=1.0)
		a1.add_label("A [dB]: "+str(round(-abs(A),4)),loc='bot',ofst=2.0)
		a1.add_label("P [dBm]: "+str(round(P,4)),loc='bot',ofst=3.0)
		a1.add_label("Tcas [K]: "+str(round(Tf,4)),loc='bot',ofst=4.0)
		a1.add_label("T_N [K]: "+str(round(Tf*Gprod,4)),loc='bot',ofst=5.0)
		#a1.add_label("NSD [dBm/Hz]: "+str(round(NSD,2)),loc='bot',ofst=5.0)
		Pd=1e-3*10.0**(S/10.0)*(1-10**(-A/10))*ntones
		if Pd>.1:
			a1.add_label("%i tone P_diss [W]: "%ntones+str(round(Pd,2)),loc='bot',ofst=6.0)
		elif Pd>.001:
			a1.add_label("%i tone P_diss [mW]: "%ntones+str(round(Pd*1e3,2)),loc='bot',ofst=6.0)
		elif Pd>1e-6:
			a1.add_label("%i tone P_diss [uW]: "%ntones+str(round(Pd*1e6,2)),loc='bot',ofst=6.0)
		else:
			a1.add_label("%i tone P_diss [nW]: "%ntones+str(round(Pd*1e9,2)),loc='bot',ofst=6.0)
		a1.add_label("Length [mm]: "+str(round(L,2)),loc='bot',ofst=7.0)
		if therm=='Undefined':
			a1.add_label("Thermal @ %.2f [W]: "%(min(Tin,Tout)) + therm,loc='bot',ofst=8.0)
		elif therm>0.1:
			a1.add_label("Thermal @ %.2f [W]: "%(min(Tin,Tout))+str(round(therm,2)),loc='bot',ofst=8.0)
		elif therm>1e-3:
			a1.add_label("Thermal @ %.2f [mW]: "%(min(Tin,Tout))+str(round(therm*1e3,2)),loc='bot',ofst=8.0)
		elif therm>1e-6:
			a1.add_label("Thermal @ %.2f [uW]: "%(min(Tin,Tout))+str(round(therm*1e6,2)),loc='bot',ofst=8.0)
		else:
			a1.add_label("Thermal @ %.2f [nW]: "%(min(Tin,Tout))+str(round(therm*1e9,2)),loc='bot',ofst=8.0)
		d.add(e.LINE, d='right', l=5)
		return p

##########################
# Optimized cable length function ####Requires thermal model of cable
##########################
def addCableOptimized(Tin, Tout, Lmin, Lmax, fin, ntones=1000, ctype="SC-086/50-SS-SS",S=0.0,p=0):
		"""
		Params:
		Tin : temperature of stage of input [K]
		Tout : temperature of stage of output [K]
		Lmax : Maximum physical length of cable in the event power dissipated by attenuation is super small compared to thermal [mm]
		T : physical temperature of the attenuation in [K]
		fin : center frequency of input signal in [MHz]
		ctype: cable type as a string, see cable.py
		S: signal power in [dBm]
		p : parameter array for components
		"""
		###For initial set-up
		L=1
		T,therm=cableCond(ctype,Tin,Tout,L)
		A = getLoss(fin,T,L,ctype)
		if therm=='Undefined':
			return '%.2f-%.2f: No thermal model to optimize, removing this cable from line....'%(Tin,Tout)
		if np.size(np.shape(p)) != 0:
			S = p[-1][3]
		###Can calculate once then scale for optimal length
		W_m = therm*(L/1000)	###Thermal power dissipation W*m
		dBpm = A/(L/1000)		###Attenuation dB per m
		def Net_P(L):	###L in mm
			Pd=1e-3*10.0**(S/10.0)*(1-10**(-dBpm*(L/1000)/10))*ntones
			return W_m/(L/1000)+Pd
		div=[10,1,0.1]
		Li=0.1
		for i in div:
			# print(i)
			L_searching=True
			j=0
			while L_searching:
				L_set=np.arange(0,3)*i+Li
				P_set=Net_P(L_set)
				# print(L_set,P_set)		###Check test
				if P_set[1]<P_set[0] and P_set[1]<P_set[2]:		###indicating the minimum power is within this range
					Li=L_set[0]
					L_searching=False
				else:
					Li=L_set[1]
				j+=1
				if Li>Lmax:										###Capping the max length
					Li=Lmax
					L_set=np.arange(0,3)*i+Li
					break
			if Li==Lmax:
				print('%.2f-%.2f: Maximum cable length of %.1f [mm] reached without optimization'%(Tin,Tout,Lmax))
				break
		L=L_set[0]
		if L<Lmin:
			print('%.2f-%.2f: Optimized length (%.1f [mm]) below minimum length of %.1f [mm]'%(Tin,Tout,Li,Lmin))
			L=Lmin
		
		T,therm=cableCond(ctype,Tin,Tout,L)
		A = getLoss(fin,T,L,ctype)
		NSD = 0
		Tf = 0
		Ti = 0
		P = 0
		Gprod = 0	
		if np.size(np.shape(p)) == 0:
			Tf = T*(10**(abs(A)/10.) - 1.)
			# output signal power
			P = S-abs(A) # [dBm]
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			p = np.array([S, T*(10.**(abs(A)/10.0)-1.), -abs(A), P, NSD])
		
		else:
			#Noise calculation with temperature modified for attenuator T_equiv = T_phys*(Atten-1)
			p = np.vstack((p,np.array([S,T*(10.**(abs(A)/10.0)-1.),-abs(A),P,0])))
			# output signal power
			S,p[-1][0] = p[-2][3], p[-2][3]
			P = p[-1][0]-abs(A) # [dBm]
			p[-1][3] = P  # assign output power to new array      
			Gc = 10.0**(p[:,2]/10.0) # component gain array and convert to linear
			#print(Gc)
			p_flat = p.flatten()
			for i in range(int(len(p_flat)/nParams)): # Friis cascade noise loop
					Ti = p[i][1] # grab specific components noise, or thermal temp
					if i == 0:
						Tf += Ti
					else:
						Gprod = np.prod(Gc[:i])
						Tf += Ti / Gprod
			NSD =10.0*np.log10( 1.38e-23*Tf*1000.0) # noise spectral density [dBm/Hz]
			# add NSD for component
			p[-1][4] = NSD
			Gprod = np.prod(Gc[:i+1])
		# draw component and label with params S, T, G, P, Tcas, NSD
		a1 = d.add(e.CABLE, d='right')
		a1.add_label(ctype+" Cable",loc='top')
		a1.add_label("S [dBm]: "+str(round(S,4)),loc='bot')
		a1.add_label("(%.2fK-%.2fK), Tavg [K]: "%(Tin,Tout)+str(round(T,2)),loc='bot',ofst=1.0)
		a1.add_label("A [dB]: "+str(round(-abs(A),4)),loc='bot',ofst=2.0)
		a1.add_label("P [dBm]: "+str(round(P,4)),loc='bot',ofst=3.0)
		a1.add_label("Tcas [K]: "+str(round(Tf,4)),loc='bot',ofst=4.0)
		a1.add_label("T_N [K]: "+str(round(Tf*Gprod,4)),loc='bot',ofst=5.0)
		#a1.add_label("NSD [dBm/Hz]: "+str(round(NSD,2)),loc='bot',ofst=5.0)
		Pd=1e-3*10.0**(S/10.0)*(1-10**(-A/10))*ntones
		if Pd>.1:
			a1.add_label("%i tone P_diss [W]: "%ntones+str(round(Pd,2)),loc='bot',ofst=6.0)
		elif Pd>.001:
			a1.add_label("%i tone P_diss [mW]: "%ntones+str(round(Pd*1e3,2)),loc='bot',ofst=6.0)
		elif Pd>1e-6:
			a1.add_label("%i tone P_diss [uW]: "%ntones+str(round(Pd*1e6,2)),loc='bot',ofst=6.0)
		else:
			a1.add_label("%i tone P_diss [nW]: "%ntones+str(round(Pd*1e9,2)),loc='bot',ofst=6.0)
		if Li==Lmax:
			a1.add_label("Length [mm]: "+str(round(L,2)),loc='bot',ofst=7.0)
		elif Li<Lmin:
			a1.add_label("Length [mm]: "+str(round(L,2)),loc='bot',ofst=7.0)
		else:
			a1.add_label("OPTIMIZED Length [mm]: "+str(round(L,2)),loc='bot',ofst=7.0)
		if therm=='Undefined':
			a1.add_label("Thermal @ %.2f [W]: "%(min(Tin,Tout)) + therm,loc='bot',ofst=8.0)
		elif therm>0.1:
			a1.add_label("Thermal @ %.2f [W]: "%(min(Tin,Tout))+str(round(therm,2)),loc='bot',ofst=8.0)
		elif therm>1e-3:
			a1.add_label("Thermal @ %.2f [mW]: "%(min(Tin,Tout))+str(round(therm*1e3,2)),loc='bot',ofst=8.0)
		elif therm>1e-6:
			a1.add_label("Thermal @ %.2f [uW]: "%(min(Tin,Tout))+str(round(therm*1e6,2)),loc='bot',ofst=8.0)
		else:
			a1.add_label("Thermal @ %.2f [nW]: "%(min(Tin,Tout))+str(round(therm*1e9,2)),loc='bot',ofst=8.0)
		d.add(e.LINE, d='right', l=5)
		return p


def addLEKID(arrayName):
		"""
		Params:
		arrayName: String of LEKID array
		"""
		# draw component and label
		a1 = d.add(e.LEKID,d = 'right')
		a1.add_label(arrayName,loc='top')
		d.add(e.LINE, d='right', l=4)
		return

def addSource(T,S=0.):		###Place first
		"""
		Params:
		T : physical temperature of source[K]
		A: attenuation value in [dB]
		S: signal power in [dBm]
		p : parameter array for components
		"""
		Tf = T	# effective noise temp
		# output signal power
		P = S # [dBm]
		p = np.array([S, Tf, 0, P, 0])
		
		# draw component and label with params S, T, G, P, NSD
		a1 = d.add(e.SOURCE,d = 'right')
		a1.add_label("Source",loc='top')
		a1.add_label("S [dBm]: "+str(round(S,4)),loc='bot')
		a1.add_label("T [K]: "+str(round(T,4)),loc='bot',ofst=1.0)
		a1.add_label("P [dBm]: "+str(round(P,4)),loc='bot',ofst=3.0)
		a1.add_label("Tcas [K]: "+str(round(Tf,4)),loc='bot',ofst=4.0)
		d.add(e.LINE, d='right', l=5)
		return p

