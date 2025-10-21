import numpy as np
import sys

file = sys.argv[1]
data = []
with open(file, 'r') as f:
    for line in f:
        line = line.strip()
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            th = float(parts[0])
            I = float(parts[1])
            data.append((th, I))
        except:
            pass

data = np.array(data)
I_vals = data[:,1]
#in case theta is in rad, comment below
thetas_rad = data[:,0]
#uncomment here
#thetas_deg = data[:,0]
#thetas_rad = np.deg2rad(thetas_deg)

#get cross section from diff cross section
sin_theta = np.sin(thetas_rad)
integrand = I_vals * sin_theta
sigma_tot = 2 * np.pi * np.trapz(integrand, thetas_rad) #same unit as input
r_eff = np.sqrt(sigma_tot / np.pi)  # same unit as input

print(np.round(r_eff*1e6,5))#,np.round(sigma_tot*1e12,3))
