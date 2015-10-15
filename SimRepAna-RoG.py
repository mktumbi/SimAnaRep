import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline


def Rgyr(u):
    """This function will plot Radius of gyration for a give universe (trajectory)
    :input
        universe: A universe with both PSF and DCD Files.
    :return 
        1) matplotlib object for figure. 
        2) array with data. 
    """
    
    Rgyr = []
    protein = u.select_atoms("protein")
    for ts in u.trajectory:
       Rgyr.append((u.trajectory.time, protein.radius_of_gyration()))
    Rgyr = np.array(Rgyr)

    #print Rgyr
    #ax = plt.subplot(111)
    ax.plot(Rgyr[:,0], Rgyr[:,1], 'r--', lw=2, label=r"$R_G$")
    ax.set_title("Radius of Gyration")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
    #ax.figure.savefig("Rgyr.pdf")
    #plt.draw()
    return  ax, Rgyr