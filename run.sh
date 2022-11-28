#!/bin/bash

#SBATCH -p gpu       # Partición (cola)
#SBATCH -N 1                # Numero de nodos
#SBATCH -n 47         # Numero de cores(CPUs)
#SBATCH --mem 16GB
#SBATCH -t 2-00:00     # Duración (D-HH:MM)
#SBATCH -o %j.out  #STDOUT
#SBATCH -e %j.err   #STDERR

python main.py reflectores