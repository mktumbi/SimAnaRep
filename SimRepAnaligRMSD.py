import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis.rms import rmsd

def ligRMSD(u,ref):
    """
    This function produces RMSD data and plots for ligand. 
    :input 
        1) Universe of Trajectory
        2) reference universe
    :return
        1) matplot object
        2) array for RMSD data.
        
    """
    RMSD_lig = []
    ligand = u.select_atoms("(resid 142:146) and not name H*") ## include selection based on user description
    #current = u.select_atoms("segname BGLC and not name H*")
    reference = ref.select_atoms("(resid 142:146) and not name H*")
    for ts in u.trajectory:
        A = ligand.coordinates()
        B = reference.coordinates()
        C = rmsd(A,B)
        RMSD_lig.append((u.trajectory.frame, C))
    RMSD_lig = np.array(RMSD_lig)
    #print RMSD_lig
    import matplotlib.pyplot as plt
    ax = plt.subplot(111)
    ax.plot(RMSD_lig[:,0], RMSD_lig[:,1], 'r--', lw=2, label=r"$R_G$")
    ax.set_xlabel("Frame")
    ax.set_ylabel(r"RMSD of ligand ($\AA$)")
    #ax.figure.savefig("RMSD_ligand.pdf")
    #plt.draw()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc = 'lower left')
    return ax,  RMSD_lig
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This function will plot RMSD for a given universe (trajectory).')
    parser.add_argument('-j', '--jobname', help='Enter your job name and it will appear as first coloumn in the result file', default='Test')
    parser.add_argument('-trj', '--trajectory', help='Filename of Trajecotry file.', required=True)
    parser.add_argument('-top', '--topology', help='Filename of psf/topology file', required=True)
    args = parser.parse_args()
    
    u = MDAnalysis.Universe(args.topology, args.trajectory)
    ref = MDAnalysis.Universe(args.topology, args.trajectory)
    ligandRMSD = []
    fig,ligandRMSD = ligRMSD(u,ref)
    #print caRMSD
    np.savetxt(args.jobname+"_ligRMSD.data", ligandRMSD)
    fig.figure.savefig(args.jobname+"_ligRMSD.pdf")