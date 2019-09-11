#
def powerResidue(N, seed=None, a=273673163155, c=13, M=2**48, to_plot = False):
    """ Calculate a series of random numbers
    """
    import datetime
    import matplotlib.pyplot as plt
    import numpy as np
    if seed == None:
        print("Seed value set to NONE, defaulting to system time.")
        seed=int(datetime.datetime.now().strftime("%Y%m%d%H%M%s"))
    else:
        pass

    r = seed
    rand = []
    for i in range(N):
        rand.append(((a*r + c) % M)/M)

        r = (a*r + c) % M
        
    if to_plot == True:
        
        plt.plot(np.arange(0,N,1),rand[0:N])
        #plt.show()
    return rand[0:N]
