#!/bin/bash
#SBATCH --time=RUN_TIME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NCPUS
#SBATCH --mem=MEM_ALLOC
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=JOB_NAME
#SBATCH --output=../output/job_logs/%x-%j.out

# INPUT=../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox000-000_m200b_halos.list
# M=300
# OUTPUT=../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox000-000/m200b_MASS_vp_n1_c2.dat

echo $INPUT $M $OUTPUT| ../code/01-calc_vp_halo_llist_pbc.exe
