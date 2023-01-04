#!/bin/bash

#SBATCH -p gpu       # Partición (cola)
#SBATCH -N 1                # Numero de nodos
#SBATCH -n 55
#SBATCH -t 2-00:00     # Duración (D-HH:MM)
#SBATCH -o %j.out  #STDOUT
#SBATCH -e %j.err   #STDERR

python main_rev.py tejido