function out = dBeh4(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 5.3 on page 1119 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	% J. Op. and App., 2019.  
	% --------------------------------------------------------
	J1 = [-2*x(1) 1]; 
	J2 = [2*x(1) 1]; 
	out = [J1;J2];
