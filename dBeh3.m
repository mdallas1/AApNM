function out = dBeh3(x)
	% --------------------------------------------------------
	% Implmentation of Jacobian Ex. 5.3 on page 1119 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	% J. Op. and App., 2019.  
	%
	% --------------------------------------------------------
	J1 = [-(1/9)*sin(x(1))-x(2).*cos(x(1)), -sin(x(1))];
	J2 = [(1/9)*cos(x(1))-x(2).*sin(x(1)), cos(x(1))];
	out = [J1;J2];
