import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np


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
    ax = plt.subplot(111)
    ax.plot(Rgyr[:,0], Rgyr[:,1], 'r--', lw=2, label=r"$R_G$")
    ax.set_title("Radius of Gyration")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
    #ax.figure.savefig("Rgyr.pdf")
    #plt.draw()
    return  ax, Rgyr


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This function will plot Radius of gyration for a given universe (trajectory).')
    parser.add_argument('-j', '--jobname', help='Enter your job name and it will appear as first coloumn in the result file', default='Test')
    parser.add_argument('-trj', '--trajectory', help='Filename of Trajecotry file.', required=True)
    parser.add_argument('-top', '--topology', help='Filename of psf/topology file', required=True)
    args = parser.parse_args()
    
    u = MDAnalysis.Universe(args.topology, args.trajectory)
    fig, Rgyrdata = Rgyr(u)
    np.savetxt(args.jobname+"-RoG.data", Rgyrdata)
    #print Rgyrdata
    fig.figure.savefig(args.jobname+"-RoG.pdf")