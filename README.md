# Code from ACSE publication "[Anderson Acceleration for Perturbed Newton Methods](https://doi.org/10.3934/acse.2025024)" 

## Matt Dallas 

### 28 November 2025

The code above allows the reader to reproduce the numerical 
experiments presented in Section 5.1 and Section 5.3. All 
scripts are written in MATLAB/GNU Octave (version 9.4.0). 

 - `chand_solver.m` solves the Chandrasekhar H-equation with 
		Anderson or $\gamma$-safeguarded Anderson applied to 
		Newton or Levenberg-Marquardt. The user may adjust the 
		problem parameters in the script. 
 - `Beh1_solver.m`, `Beh2_solver.m`, `Beh3_solver.m`, and 
		`Beh4_solver.m` solve the corresponding test problems 
		from Section 5.3. 

