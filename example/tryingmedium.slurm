#!/usr/bin/env bash
#SBATCH --ntasks 2
#SBATCH --mem 16G
#SBATCH --time 0-00:45:00
#SBATCH --partition standard
#SBATCH --account biol4585

echo "working!"

julia Joseph-Viterbi-First-Implementation.jl observations.dat medium.dat
