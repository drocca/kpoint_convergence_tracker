# Kpoint convergence tracker

Python package to automatically test k-point convergence in crystalline compounds.

The code can be used as 
./find_optimal_k  </path/pw_input_file>  \<cutoff>

WARNING: Cutoff denotes the wavefunction cutoff (ecutwfc). In ultrasoft or PAW calculations the density cutoff (ecutrho) will not be modified in pw_input_file; it is important to propose a value of cutoff compatible with the ecutrho provided in pw_input_file.

If find_optimal_k is not recognized as an executable you will need to run the command “chmod +x find_optimal_k”. In order to use this script on your machine it is necessary to edit the file run_qe.py to add the correct instructions to run pw.x on your computer. This script also requires the f90nml module, that can be installed with the command "pip install f90nml".  

At this stage the code uses only the total energy divided by the number of electrons as a convergence criterion. When this energy changes by less than 5e-5 Ry the code stops and provides the converged k-point grid.  
