#
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def walk_2D(N):
    """ Function to compute an N-step random walk
    
        Input:
            N  ::  Total number of steps

        Output:
            x  ::  Array of all x positions
            y  ::  Array of all y positions

    """
    # seed random number generator
    rand.seed() #this initializes the first value in the rand.random func'n. 
    

    # initialize x, y
    x = [0.0] 
    y = [0.0] #creates x and y array where 0,0 is the first corrdinate in the walk

    # step in x-y space N times
    for n in range(N):
        x.append(x[-1] + (rand.random() - 0.5)*2.0) #appends the x and y array with the random value
        y.append(y[-1] + (rand.random() - 0.5)*2.0)
    d_origin = np.sqrt(x[-1]**2+y[-1]**2)
    
    return np.array(x), np.array(y), d_origin


# Example simulation
walker_1 = walk_2D(1000)   # compute path for 1000 steps

# Example plot of (x, y) pairs from example simulation 
plt.plot(walker_1[0], walker_1[1], '-')
plt.xlabel("X")
plt.ylabel("Y")
plt.text(np.min(walker_1[0])+1, np.min(walker_1[1])+1, '$d_t$ = {0}'.format(walker_1[2]))
plt.show()

for q in range(5):
    walker = walk_2D(1000)
    plt.plot(walker[0], walker[1], '-', alpha = 0.5)
    plt.xlabel("X")
    plt.ylabel("Y")
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 03/rand_walk.pdf')
plt.show()

d_walker = []
for q in range(1000):
    walker = walk_2D(1000)
    plt.plot(walker[0], walker[1], '-', alpha = 0.5)
    d_walker.append(walker[2])
    plt.xlabel("X")
    plt.ylabel("Y")
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 03/rand_walk_100.pdf')
plt.show()

plt.hist(d_walker)
#plt.text(3,0.5, '$d_a = {0}$'.format(np.mean(d_walker)))
#plt.text(3,1.5, '$d_s = {0}$'.format(np.std(d_walker)))
#plt.text(3,2.5, '$d_m = {0}$'.format(np.median(d_walker)))
plt.axvline(np.median(d_walker), color = '#8258FA', linestyle = 'dashed')
plt.axvline(np.median(d_walker) - np.std(d_walker), color = '#A9D0F5', linestyle = 'dashed')
plt.axvline(np.median(d_walker) + np.std(d_walker), color = '#A9D0F5', linestyle = 'dashed')
plt.title('distance walker travelled from origin : 2D')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 03/walker_hist.pdf')
plt.show()
    
def walk_3D(N):
    """ Function to compute an N-step random walk
    
        Input:
            N  ::  Total number of steps

        Output:
            x  ::  Array of all x positions
            y  ::  Array of all y positions
            z  ::  Array of all z positions
    """
    # seed random number generator
    rand.seed() #this initializes the first value in the rand.random func'n. 
    

    # initialize x, y
    x = [0.0] 
    y = [0.0] #creates x and y array where 0,0 is the first corrdinate in the walk
    z = [0.0]
    # step in x-y space N times
    for n in range(N):
        x.append(x[-1] + (rand.random() - 0.5)*2.0) #appends the x and y array with the random value
        y.append(y[-1] + (rand.random() - 0.5)*2.0)
        z.append(z[-1] + (rand.random() - 0.5)*2.0)
    d_origin = np.sqrt(x[-1]**2+y[-1]**2 + z[-1]**2)
    
    return np.array(x), np.array(y), np.array(z), d_origin


fig = plt.figure()
ax = plt.subplot(111, projection='3d')
d_travel = []
for w in range(1000):
    walker = walk_3D(1000)
    ax.plot(walker[0], walker[1], walker[2])
    d_travel.append(walker[3])
plt.title('distance walker travelled from origin : 3D')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 03/walk_3D.pdf')
plt.show()

plt.hist(d_travel)
#plt.text(25,0.5, '$d_a = {0}$'.format(np.mean(d_travel)))
#plt.text(25,1.5, '$d_s = {0}$'.format(np.std(d_travel)))
#plt.text(25,2.5, '$d_m = {0}$'.format(np.median(d_travel)))
plt.axvline(np.median(d_travel), color = '#8258FA', linestyle = 'dashed')
plt.axvline(np.median(d_travel) - np.std(d_travel), color = '#A9D0F5', linestyle = 'dashed')
plt.axvline(np.median(d_travel) + np.std(d_travel), color = '#A9D0F5', linestyle = 'dashed')
plt.title('distance walker travelled from origin : 3D')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 03/walker_hist3D.pdf')
plt.show()