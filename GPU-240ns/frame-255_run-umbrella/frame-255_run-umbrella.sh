#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --gres=gpu:1
#SBATCH --job-name=PMF_Vesicle
#SBATCH --account=nawimem
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out

# Now list your module

ml use -a /nopt/nrel/apps/modules/test/modulefiles
ml cuda/11.2 openmpi/4.1.3/gcc-11.3.0-cuda-11.7 gcc/8.4.0 fftw/3.3.8/intel-impi intel-mpi/2020.1.217

# Short equilibration
/home/dtnguyen28/Software/GROMACS/gromacs-2022/bin/gmx_mpi grompp -f npt_umbrella.mdp -c conf255.gro -r conf255.gro -p system.top -o npt255.tpr  -n index.ndx -maxwarn 10
srun /home/dtnguyen28/Software/GROMACS/gromacs-2022/bin/gmx_mpi mdrun -v -deffnm npt255 

# Umbrella run
/home/dtnguyen28/Software/GROMACS/gromacs-2022/bin/gmx_mpi grompp -f md_umbrella.mdp -c npt255.gro -r npt255.gro -t npt255.cpt -p system.top -o umbrella255.tpr -n index.ndx -maxwarn 10
srun /home/dtnguyen28/Software/GROMACS/gromacs-2022/bin/gmx_mpi mdrun -v -deffnm umbrella255

echo "Job is done, thanks Dan!"
