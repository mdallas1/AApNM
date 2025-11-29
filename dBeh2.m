function out = dBeh2(x)
	% --------------------------------------------------------
	% Implmentation of Jacobian of Ex. 5.2 on page 1118 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	% J. Op. and App., 2019.  
	%
	% --------------------------------------------------------
	J1 = [3*x(1).^2-x(2), -x(1)];
	J2 = [3*x(1)^2 + x(2), x(1)];
	out = [J1;J2];
