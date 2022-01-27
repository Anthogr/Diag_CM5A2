# Load modules and run script BasicDiag_CMA5A2.py 

# The two following lines unload all modules and reload specific ones
# compatible with the used libraries in BasicDiag_CMA5A2.py
ml purge
ml load intel/19.0.5.281 mpi/openmpi/4.0.5 python3/3.7.5

# Run script
python3 Diag_CMA5A2.py
