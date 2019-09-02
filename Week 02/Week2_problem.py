#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 09:39:12 2019

@author: amandaash
"""
import numpy as np
import matplotlib.pyplot as plt
#equations:
#d_x = v_0,x*t + x_0
# d_y = v_0,y*t - 1/2gt^2 + y0
#vf^2 = v0^2 + 2ad

def projectile(initial_velocity, angle, y_initial, x_initial, accelaration_gravity, time, mass, drag_coefficient, delta_t, plot = False):
    """
    Input:
        initial_velocity [m/s]
        angle [degrees]
        initial y position [m]
        initial x position [m]
        acceleration due to gravity [m/s^2]
        time - array [s]
        mass [kg]
        drag coefficient [kg/m]
        time step [s]
        plot [True/False]
    Returns:
        ideal trajectory array [time, distance x, distance y, velocity x, velocity y]
        drag trajectory array [time, distance x, distance y, velocity x, velocity y]
        plot of distance y vs. distance in x (if plot condition = True)
        """
    #trajectory in the y:
    y0 = y_initial
    x0 = x_initial
    g = accelaration_gravity #[m/s^2]
    theta = (angle*np.pi)/180
    vy_0 = initial_velocity * np.sin(theta)#[m/s]
    vx_0 = initial_velocity*np.cos(theta)#[m/s]
    A = drag_coefficient/mass
    
    vy_ideal = vy_0 - g*time
    y_ideal = y0 + (vy_0*time) + ((1/2)*(-g)*time**2)
    
    vx_ideal = np.full(np.shape(time), vx_0)
    x_ideal = x0 + (vx_0*time)
    
    vx_id = []
    x_id = []
    vy_id = []
    y_id = []
    time_id = []
    
    ideal_motion = np.column_stack((x_ideal, y_ideal, vx_ideal, vy_ideal, time))
    for line in ideal_motion:
        if np.sign(line[1]) == 1:
            vx_id.append(line[2])
            x_id.append(line[0])
            vy_id.append(line[3])
            y_id.append(line[1])
            time_id.append(line[4])
    
    #beta = (np.arctan((np.sqrt(A)*vy_0)/np.sqrt(g)))/(np.sqrt(A)*np.sqrt(g))
    #vy_drag = (np.sqrt(g)*np.tan((time+beta)*np.sqrt(A)*np.sqrt(g)))/np.sqrt(A)
    #alpha = y0 + ((mass*np.log(np.abs(np.cos(np.sqrt(A)*np.sqrt(g)*beta))))/drag_coefficient)
    #alpha = y0 + ((np.log(np.cos(np.sqrt(A)*np.sqrt(g)*(time+beta))))/np.sqrt(A))
    #y_drag = ((-mass*np.log(np.abs(np.cos(np.sqrt(A)*np.sqrt(g)*(time+beta)))))/drag_coefficient)+alpha
    #y_drag = alpha = (-np.log(np.cos(np.sqrt(A)*np.sqrt(g)*(time+beta))))/np.sqrt(A) + alpha
    #a_drag_y = (-drag_coefficient*vy**2)/mass
    #vy_drag = (np.sqrt(g)*np.tan(np.sqrt(A)*np.sqrt(g)*time))/(np.sqrt(A))
    #y_drag = ((-1)/A)*np.log(np.cos(np.sqrt(A)*np.sqrt(g)*time))+y0
    
    t0 = 0
    vx = []
    vy = []
    dx = []
    dy = []
    time_drag = []
    for t in time:
        
        #y-direction:
        y0 = y0 + (vy_0*delta_t)
        if np.sign(y0)==-1:
            continue
        a_y0 = -g + ((-drag_coefficient*vy_0**3)/(mass*(np.abs(vy_0)+0.0000000000001)))
        vy_0 = vy_0 + (a_y0*delta_t)
        
        #x-direction:
        x0 = x0 + (vx_0*delta_t)
        a_x0 = -(drag_coefficient*vx_0**2)/mass
        vx_0 = vx_0 + (a_x0*delta_t)
     
        vx.append(vx_0)
        dx.append(x0)
        vy.append(vy_0)
        dy.append(y0)
        time_drag.append(t)
    
    #trajectory in the x:
    
    #a_drag_x = (-drag_coefficient*vx_0**2)/mass
    #x_drag = x0 + (vx_0*time) + ((1/2)*(a_drag_x)*(time**2))
    #vx_drag = -1/((-A*time)-(1/vx_0))
    #x_drag = ((1/A)*np.log((-A*time)-(1/vx_0))) + (A*(x0/(mass*np.log(-1/vx_0))))
    #x_drag = ((-np.log((A*time)+(1/vx_0)))/A) + ((mass*np.log(1/vx_0))/drag_coefficient) + x0
    #vx_drag = (-1)/((-A*time)+vx_0)
    #x_drag = (1/A)*(np.log((-A*time)+vx_0))+x0
    if plot == True:
        plt.plot(x_id,y_id, '.', label = 'Ideal', color = '#81BEF7')
        plt.plot(dx, dy, '.', label = 'Drag', color = '#5F04B4')
        #plt.suptitle(r'initial_velocity = {0}, $\theta$ = {1}, $x_0$ = {2}, \n$y_0$ = {3}, mass = {4}, $\Delta$t = {5}'.format(initial_velocity, angle, y0,x0,mass, delta_t))
        plt.title('trajectory given: initial velocity={0}[m/s], $\Theta ={1}^\circ$, $x_0$={2}[m], \n $y_0$={3}[m], $a_g$={4}$[m*s^-2]$, mass={5}[kg], \n $\Delta t={6}[s]$, $C_d$={7}$[kg*m^-1]$'.format(initial_velocity, angle, x_initial,y_initial,g, mass, delta_t,drag_coefficient))
        #plt.hlines(0, np.min(x_ideal)-3, np.max(x_ideal)+3)
        plt.ylabel('y(t) [m]')
        plt.xlabel('x(t) [m]')
        plt.legend(loc = 'upper right')
        plt.tight_layout()
        plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 02/trajectory_plot.pdf')
        plt.show()
       
    return (np.column_stack((time_id, x_id, y_id, vx_id, vy_id))), (np.column_stack((time_drag,dx,dy,vx,vy)))

ideal, drag = projectile(30,45,0,0,9.81,np.arange(0,100,0.001),2,0.02,0.002, True)
np.savetxt('/Users/amandaash/Desktop/PHYS_3210/Week 02/drag_trajectory.txt', drag, header = 'time, distance x, distance y, velocity x velocity y')