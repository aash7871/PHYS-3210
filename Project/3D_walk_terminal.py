import numpy as np
import matplotlib.pyplot as plt
import random as rand
from cycler import cycler
from mpl_toolkits import mplot3d
import glob as glob
import matplotlib.patches as mpatches

def density_delta(radius, r_star):
    p0 = 100
    frac_radius = radius/r_star
    p = p0*(1-frac_radius**2)
    return p

G = 6.67e-8
k = 1.38e-16
mu = 1.3
mp = 1.67e-24

def pressure_delta(radius, central_density, R_star): 
    
    P_c = (4/15)*np.pi*G*(central_density**2)*(R_star**2)
    term1 = (radius**2)/6
    term2 = (2*radius**4)/(15*R_star**2)
    term3 = (radius**6)/(30*R_star**4)
    
    P_r = P_c - (4*np.pi*G*(central_density**2)*(term1-term2+term3))
    
    return P_r

def temp_alpha(density, pressure):

    T_r = (mu*mp*pressure)/(k*density)
    
    return T_r

def rosseland_mean_opacity(density, temperature):
    X = 0.7381
    Y = 0.2485
    Z = 0.0134
    
    C_ff = 3.68e18
    gff = 1
    
    K_ff = C_ff*(X+Y)*(1+X)*(density/(temperature**3.5))
    
    C_bf = 4.34e21
    g_bf = 1
    t = 1
    
    K_bf = C_bf*(g_bf/t)*Z*(1+X)*(density/(temperature**3.5))
    
    C_es = 0.02
    K_es = C_es*(1+X)
    
    #we're able to ignore H- opacity
    """
    if temperature <= 6000 and density <= 1e-5:
        C_H_ = 7.9e-34
        K_H_ = C_H_*(Z/0.02)*(density**-0.5)*(temperature**9)
    elif temperature >= 3000 and density >= 1e-10:
        C_H_ = 7.9e-34
        K_H_ = C_H_*(Z/0.02)*(density**-0.5)*(temperature**9)
    else:
        K_H_ = 0
    
    #K = np.mean([K_ff, K_bf, K_es]) 
    #K = ((K_ff**-1)+(K_bf**-1)+(K_es**-1))**-1
    """
    K = K_ff+K_bf+K_es
    
    return K

def mean_free_path_beta(opacity, density):
    l = 1/(density*opacity)
    return l

def mean_free_path_gaussian(opacity, density, sigma):
    l = 1/(density*opacity)
    mfp = np.random.normal(l, sigma, 1)
    return float(mfp[0])

def photon_escape(N_walkers, N_layers, N_steps):
    c = 3e10
    sigma = 1
    
    x = 0
    y = 0
    z = 0

    x_val = [0]
    y_val = [0]
    z_val = [0]

    r_sample = 0.1
    r_sun = 6.957e10
    section = r_sun/N_layers

    epsilon= 1e-7

    radius_bisection = np.arange(epsilon, r_sun, section)

    escape_times = []
    for radius in radius_bisection:
        layer_escape = []
        for walker in range(N_walkers):
            x = radius
            y = 0
            z = 0
            t = 0
            r = radius
            x_walk = [x]
            y_walk = [y]
            z_walk = [z]

            for n in range(N_steps):
                density = density_delta(r, r_sun)
                pressure = pressure_delta(r, 100, r_sun)
                temperature = temp_alpha(density, pressure)
                opacity = rosseland_mean_opacity(density, temperature)

                theta = rand.uniform(0,2*np.pi)
                theta_walk.append(theta)
                phi = rand.uniform(0,2*np.pi)
                phi_walk.append(phi)
            
                if gaussian == True:
                    l = mean_free_path_gaussian(opacity, density, sigma)
                else:
                    l = mean_free_path_beta(opacity, density)
            
                l_all.append(l)
            
                dx = l * np.sin(phi) * np.cos(theta)
                dy = l * np.sin(phi) * np.sin(theta)
                dz = l * np.cos(phi)
           
                x = x + dx
                y = y + dy
                z = z + dz
            
                r = np.sqrt(x**2 + y**2 + z**2)
            
                x_walk.append(x)
                y_walk.append(y)
                z_walk.append(z)
            
                r = np.sqrt((x**2) + (y**2) + (z**2))

                d = r - radius

                dt = l/c
                dt_all.append(dt)
                t = t + dt
   
                if d>(r_sample):
                    escape_step = n
                    if escape_step == 0:
                        escape_time = np.mean(dt_all)
                        layer_escape.append(escape_time/((3.14e9)))
                        
        
                    else:
                        n_steps_escape = escape_step*(((section/r_sample))**2)
                        escape_time = n_steps_escape*np.mean(dt_all)
                        layer_escape.append(escape_time/((3.14e9)))
                        
                    break
        
        escape_times.append(np.mean(layer_escape))   
        
    return np.sum(escape_times)


file_gaussian = open('escape_times_gaussian.txt', '+w')
file_non = open('escape_times_non.txt', '+w')
N_tests = 100
N_walkers = 1000
N_layers = 10
N_steps = int(1e5)

escape_times_gaussian = []
escape_times_non = []
for test in range(N_tests):
    escape_g = photon_escape(N_walkers, N_layers, N_steps)
    escape_n = photon_escape(N_walkers, N_layers, N_steps, False)
    escape_times_gaussian.append(escape_g)
    escape_times_non.append(escape_n)
    
np.savetxt(file_gaussian, np.array(escape_times_gaussian))
np.savetxt(file_non, np.array(escape_times_non))
file.close()