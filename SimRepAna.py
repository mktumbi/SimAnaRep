import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis.rms import rmsd
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from  SimAnaRepproRMSD import *
from SimAnaRepRoG import *
from SimRepAnaligRMSD import *


parser = argparse.ArgumentParser(description='SimRepAna is a programme to generate report for your Ligand-Protein MD Simulations.')
parser.add_argument('-j', '--jobname', help='Enter your job name and it will appear as first coloumn in the result file', default='Test')
parser.add_argument('-trj', '--trajectory', help='Filename of Trajecotry file.', required=True)
parser.add_argument('-top', '--topology', help='Filename of psf/topology file', required=True)
args = parser.parse_args()

u = MDAnalysis.Universe(args.topology, args.trajectory)
ref = MDAnalysis.Universe(args.topology, args.trajectory)



ligandRMSD = []
ligfig,ligandRMSD = ligRMSD(u,ref)
np.savetxt(args.jobname+"_ligRMSD.data", ligandRMSD)
ligfig.figure.savefig(args.jobname+"_ligRMSD.pdf")
caRMSD =[]
allRMSD = []
profig,caRMSD,allRMSD = proRMSD(u,ref)
np.savetxt(args.jobname+"-caRMSD-pro.data", caRMSD)
np.savetxt(args.jobname+"-allRMSD-pro.data", allRMSD)
profig.figure.savefig(args.jobname+"-proRMSD.pdf")

Rogfig, Rgyrdata = Rgyr(u)
np.savetxt(args.jobname+"-RoG.data", Rgyrdata)
Rogfig.figure.savefig(args.jobname+"-RoG.pdf")


    

