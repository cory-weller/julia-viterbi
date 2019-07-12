#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time=02:00:00
#SBATCH -A biol4585
#SBATCH --output=viterbi_run_on_medium_6-28_output.txt
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=16GB

module purge
module load julia

julia Viterbi-First-Implementation.jl observations.dat medium.dat
