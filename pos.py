import os 
import shutil
import sys
import subprocess

"""
From each location of Janus particle to equilibrate the system
and then measure the force between the Janus particle and solvent.
"""

def inputfile(i):
    content = "units lj \n" +\
              "atom_style bond \n" +\
              "boundary p p f \n" +\
              "processors * * * \n" +\
              "read_restart res.hollow.loc.{0}0000\n\n".format(i) +\
              "reset_timestep 0 \n" +\
              "timestep 0.005 \n" +\
              "neighbor 1 multi \n" +\
              "neigh_modify page 10000000 one 20000\n" +\
              "comm_modify mode multi \n " +\
              "\n\n" +\
              "pair_style lj/cut/opt 1.122462\n" +\
              "pair_coeff 1 1 1.0 1.0 2.5\n" +\
              "pair_coeff 1 2 0.2 1.0 4.0\n" +\
              "pair_coeff 1 3 0.4 1.0 4.0\n" +\
              "pair_coeff 2 2 1.0 1.0 1.122462\n" +\
              "pair_coeff 2 3 1.0 1.0 1.122462\n" +\
              "pair_coeff 3 3 1.0 1.0 1.122462\n" +\
              "pair_modify shift yes \n" +\
              "\n\n" +\
              "bond_style fene\n" +\
              "bond_coeff 1 30 1.5 1.0 1.0\n" +\
              "special_bonds fene\n" +\
              "\n\n" +\
              "balance 1.0 shift xyz 20 1.02\n" +\
              "\n\n" +\
              "group tetramer type 1 \n"+\
              "group hydrophobic type 2\n" +\
              "group hydrophilic type 3\n" +\
              "group janus type 2 3\n" +\
              "\n\n" +\
              "fix 1 all wall/lj93 zlo 0.0 2.0 1.0 3.0 zhi 100.0 2.0 1.0 0.858374 units box\n" +\
              "region 1 cylinder z 50.0 50.0 50.0 0.0 100.0\n" +\
              "fix 3 tetramer wall/region 1 lj93 2.1 1.0 2.5\n" +\
              "fix 5 tetramer langevin 0.7 0.7 10.0 724323\n" +\
              "fix 6 tetramer nve\n" +\
              "\n\n" +\
              "dump 1 all custom 40000 dump_janus_relax.1 id type xs ys zs \n" +\
              "thermo 1000\n" +\
              "thermo_style custom step atoms temp epair etotal cpu\n" +\
              "thermo_modify lost warn flush yes\n" +\
              "restart 500000 res.janus.relax\n" +\
              "run     1000000 \n\n" +\
              "undump 1\n" +\
              "compute ferg janus group/group tetramer \n" +\
              "fix 7 janus ave/time 50 20 1000 c_ferg c_ferg[1] c_ferg[2] c_ferg[3] file ferg.txt\n\n" +\
              "fix 8 hydrophobic ave/time 50 20 1000 c_ferg c_ferg[1] c_ferg[2] c_ferg[3] file ferg_phobic.txt\n\n" +\
              "fix 9 hydrophilic ave/time 50 20 1000 c_ferg c_ferg[1] c_ferg[2] c_ferg[3] file ferg_philic.txt\n\n" +\
              "dump 1 all custom 40000 dump_janus_measure.1 id type xs ys zs\n" +\
              "dump_modify 1 sort id\n" +\
              "thermo 1000\n" +\
              "thermo_style custom step atoms temp epair etotal cpu\n" +\
              "thermo_modify lost warn flush yes\n" +\
              "restart 400000 res.janus.measure \n"  +\
              "run 2000000\n"

    dst = "./pos{0}".format(i)
    ifile = open(dst + "/input.txt", 'w')
    ifile.write(content)
    ifile.close()
    
    content = "#!/bin/bash -l\n"+\
              "#SBATCH --job-name=lammps-janus-particle-job\n"+\
              "#SBATCH -N 1\n" +\
              "#SBATCH -n 32\n" +\
              "#SBATCH --time=3-00:00:00\n" +\
              "#SBATCH -p debug\n"+\
              "#SBATCH -q normal\n\n" +\
              "srun -n 32 /sharedata01/yanfeitang/bin/lmp_mpi_v01 -sf opt -in input.txt -log log.1 > outfile.1"
    
    
    slurmfile = open(dst +"/mj.janus{0}.sh".format(i), 'w')
    slurmfile.write(content)
    slurmfile.close()


for i in xrange(2, 62, 2):
    dst = "./pos{0}".format(i)
    if not os.path.isdir(dst):
        os.mkdir(dst)
    shutil.copy("./loc/res.hollow.loc.{0}0000".format(i), dst)
    inputfile(i)
    p = subprocess.Popen(['sbatch', 'mj.janus{0}.sh'.format(i)], cwd=dst)
