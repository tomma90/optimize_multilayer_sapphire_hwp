import numpy as np
import scipy
import matplotlib.pylab as plt
import scipy.constants as cost

#Rotation Matrix from Equation \ref{muller_hwp_rot}
def rotation(x):
	return np.array([[1., 0., 0., 0.], [0., np.cos(2.*x), -1.*np.sin(2.*x), 0.], [0., np.sin(2.*x), np.cos(2.*x), 0.], [0., 0., 0., 1.]])

#Birifrangent Matrix from Equation \ref{muller_hwp_rot}
def birefringent(d):
	return np.array([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., np.cos(d), -1.*np.sin(d)], [0., 0., np.sin(d), np.cos(d)]])

#Retardance from Equation \ref{eq:retardance}
def delta(nu,  no,  ne,  t):
	return 2.*np.pi*nu*abs(no-ne)*t/cost.c

df = 0.1 # Frequency step in GHz
f0 = 34. # LFT nu_0 in GHz
f1 = 161. # LFT nu_f in GHz
steps = int((f1-f0)/df+1) # Steps
f = np.linspace(f0, f1, steps)*1.e9 # Frequency range in Hz
fc = f[0]+(f[-1]-f[0])/2 # Central frequency
no = 3.047 # Refractive index ordinary axis
ne = 3.361 # Refractive index extra-ordinary axis
t0 = cost.c/abs(no-ne)/fc/2. # thickness at central frequency for retardance = np.pi
EFF = 0.98 # Target efficiency

# Define LFT bands lower and upper bounds 
i40_i = np.where(f==34.e9)[0][0]
i40_f = np.where(f==46.e9)[0][0]
i50_i = np.where(f==42.5e9)[0][0]
i50_f = np.where(f==57.5e9)[0][0]
i60_i = np.where(f==53.e9)[0][0]
i60_f = np.where(f==67.e9)[0][0]
i68_i = np.where(f==60.e9)[0][0]
i68_f = np.where(f==76.e9)[0][0]
i78_i = np.where(f==69.e9)[0][0]
i78_f = np.where(f==87.e9)[0][0]
i89_i = np.where(f==79.e9)[0][0]
i89_f = np.where(f==99.e9)[0][0]
i100_i = np.where(f==88.5e9)[0][0]
i100_f = np.where(f==101.5e9)[0][0]
i119_i = np.where(f==101.e9)[0][0]
i119_f = np.where(f==137.e9)[0][0]
i140_i = np.where(f==119.e9)[0][0]
i140_f = np.where(f==161.e9)[0][0]

eff = np.zeros(steps) # initialize efficiency vector 
phi = np.zeros(steps) # initialize phase vector
avg_eff = np.zeros(9) # initialize band-averaged efficiency vector

# Keep looping till band-averaged efficiency vector satisfies requirement (> 0.98)
while avg[0] <= EFF or avg[1] <= EFF or avg[2] <= EFF or avg[3] <= EFF or avg[4] <= EFF or avg[5] <= EFF or avg[6] <= EFF or avg[7] <= EFF or avg[8] <= EFF:
    # Random pair of angles 0 < a < 180
	a = np.random.randint(0, 1800, 2) / 10. / 180. * np.pi
	a1 = a[1] # First plate angle
	a2 = a[0] # Second plate angle
	a3 = np.deg2rad(0.) # Central plate (third) fixed
	a4 = a[1] # Fourth plate angle
	a5 = a[0] # Fifth plate angle
	
	# Loop over all frequency steps
	for i in range(steps):
		delta_0 = delta(f[i],  no,  ne,  t0) # Single pplate retardance
		biref_0 = birefringent(delta_0) # Single plate birefringent matrix
		
		# Single plate HWP matrix for different orientation
		gamma_5=np.matmul(rotation(-1.*a5), np.matmul(biref_0, rotation(1.*a5)))
		gamma_4=np.matmul(rotation(-1.*a4), np.matmul(biref_0, rotation(1.*a4)))
		gamma_3=np.matmul(rotation(-1.*a3), np.matmul(biref_0, rotation(1.*a3)))
		gamma_2=np.matmul(rotation(-1.*a2), np.matmul(biref_0, rotation(1.*a2)))
		gamma_1=np.matmul(rotation(-1.*a1), np.matmul(biref_0, rotation(1.*a1)))
        
        # Combining plates:
		gamma_tot = np.matmul(gamma_5, gamma_4)
		gamma_tot = np.matmul(gamma_tot, gamma_3)
		gamma_tot = np.matmul(gamma_tot, gamma_2)
		gamma_tot = np.matmul(gamma_tot, gamma_1)
		
		# AHWP modulation efficiency
		eff_num = np.sqrt((gamma_tot[1, 1] - gamma_tot[2, 2])**2. + (gamma_tot[1, 2] + gamma_tot[2, 1])**2.) / 4.
		eff_den = gamma_tot[0, 0]/2. + (gamma_tot[1, 1] + gamma_tot[2, 2]) / 4.
		eff[i] = eff_num / eff_den
		
		# AHWP phase
		phi[i] = np.arctan((gamma_tot[1, 2] + gamma_tot[2, 1]) / (gamma_tot[1, 1] - gamma_tot[2, 2])) / 4.
		
		# LFT band-averaged modulation efficiency
		avg[0]=sum(eff[i40_i:i40_f+1])/(i40_f+1-i40_i)
		avg[1]=sum(eff[i50_i:i50_f+1])/(i50_f+1-i50_i)
		avg[2]=sum(eff[i60_i:i60_f+1])/(i60_f+1-i60_i)
		avg[3]=sum(eff[i68_i:i68_f+1])/(i68_f+1-i68_i)
		avg[4]=sum(eff[i78_i:i78_f+1])/(i78_f+1-i78_i)
		avg[5]=sum(eff[i89_i:i89_f+1])/(i89_f+1-i89_i)
		avg[6]=sum(eff[i100_i:i100_f+1])/(i100_f+1-i100_i)
		avg[7]=sum(eff[i119_i:i119_f+1])/(i119_f+1-i119_i)
		avg[8]=sum(eff[i140_i:i140_f+1])/(i140_f+1-i140_i)
