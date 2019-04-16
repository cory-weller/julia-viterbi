#!/usr/bin/env bash
#SBATCH --ntasks 1
#SBATCH --mem 8G
#SBATCH --time 0-00:45:00
#SBATCH --partition standard
#SBATCH --account biol4585

julia Joseph-Viterbi-First-Implementation.jl observations.dat medium.dat
