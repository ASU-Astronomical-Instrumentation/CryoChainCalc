# Requires python packages: numpy, scipy, and SchemDraw (with capitalized S and D, current version 0.6.0, lowercase one didn't work with this on python3.8.5)
# and associated .py files: linkBudgetCalc.py, cable.py, and customElements.py
# cryostat detailed microwave link budget
from CryoChainCalc import *
# Create cascade of components: addSource, addAmp, addAtten, addCable, addCableOptimized, addLEKID
# first component does not pass parameters array p, recommend addSource first to include input source noise temperature into calculation
# addCable ctype choices "SC-086/50-SS-SS" "SC-219/50-SS-SS" "SC-086/50-NbTi-NbTi" "S086MMHF"
# thermal currently exists for SS and NbTi (up to 4K for NbTi)
# addCableOptimized optimizes cable lengths to minimize power dissipated from attenuation + thermal (only really works for SS where attenuation ~= thermal)

cable_type_1 = "SC-219/50-SS-SS" # cable for lower attenuation, higher thermal (mostly higher temps)
cable_type_2 = "SC-086/50-SS-SS" # cable for higher attenuation, lower thermal (mostly lower temps)
cable_type_3 = "SC-086/50-NbTi-NbTi" # cable for near-zero attenuation, lower thermal (mostly lower temps from detectors to first amplifier)
cable_type_4 = "SC-086/50-CN-CN" # cable similar to SS version

# output_pdf_filename = input("output filename: ")	# For if you want to name in terminal
output_pdf_filename="example_cryostat"

ntones=2000.
Freq=4000.

##########################################
# RF input path
###########################################
p=addSource(300,-46)
p = addCable( 300, 50, 100, Freq, ntones, ctype=cable_type_1,p=p) # 300K to 50K 
p = addCable( 50, 4, 100, Freq, ntones, ctype=cable_type_1, p=p) # 50K to 4K
p = addAtten( 4, 10.0, ntones=ntones, p=p) # attenuator on 4K
p = addCableOptimized(4, 1, 100, 1000., Freq, ntones, ctype=cable_type_2, p=p) # 4K to 1K Optimized Length
p = addAtten( 1, 6, ntones=ntones, p=p) # attenuator on 1K ring
p = addCableOptimized(1, 0.35, 50, 500, Freq, ntones, ctype=cable_type_2, p=p) # 1K to 350mK Optimized Length
p = addAtten( 0.35, 3, ntones=ntones, p=p) # attenuator on 350mK
p = addCableOptimized(0.35, 0.25, 50, 300., Freq, ntones, ctype=cable_type_2, p=p) # 350mK to 250mK Optimized Length
p = addAtten( 0.25, 6,  ntones=ntones, p=p) # attenuator on 250mK
# ARRAY 
addLEKID("Detector Array") # function does not take or return parameter array p
###########################################
# RF output path
###########################################
p = addCable( 0.25, 0.35, 100, Freq, ntones, ctype=cable_type_3,p=p) # 250mK to 350mK
p = addCable( 0.35, 1, 100, Freq, ntones, ctype=cable_type_3,p=p) # 350mK to 1K
p = addCable( 1, 4, 100, Freq, ntones, ctype=cable_type_3,p=p) # 1K to LNA(4K)
p = addAmp( 6, 25, p=p) # 4K amplifier
p = addCable( 4, 50, 100, Freq, ntones, ctype=cable_type_1, p=p) #4K to 50K 
p = addAmp( 30, 15, p=p) # 50K amplifier
p = addCable( 50, 300, 100, Freq, ntones, ctype=cable_type_1, p=p) #50K to 300K

# DRAW DIAGRAM AND PLOT
d.draw()
d.save(output_pdf_filename+".pdf")
print('saved as: '+output_pdf_filename+".pdf")
