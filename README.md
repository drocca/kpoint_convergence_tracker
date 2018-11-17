# kpoint_convergence_tracker
Python package to automatically test k-point convergence in crystalline compounds.

The code can be used as 
./find_optimal_k  "/path/pw_input_file"  "cutoff"

If find_optimal_k is not recognized as an executable you will need to run the command “chmod +x find_optimal_k”. In order to use this script on your machine you will need to edit the file run_qe.py to add the correct instructions to run pw.x on your computer. It is also necessary to install the f90nml module with the command "pip install f90nml".  

At this stage the code uses only the total energy divided by the number of electrons as a convergence criterion. When this energy changes by less than 0.00005 Ry the code stops and provides the converged k-point grid.  

For the moment the code uses only an isotropic mech.
