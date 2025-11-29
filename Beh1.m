function out = Beh1(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 5.1 on page 1117 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	%	Problems Under an Error Bound Condition*, 
	%	J. Op. and App., 2019.  
	%
	% Has a set of non-isolated minimizers, 
	% the local error bound holds, and 
	% the Jacobian is rank 1 everywhere except the origin. 
	% 
	% x0 = [0;sqrt(5)+0.03];
	% --------------------------------------------------------
	y1 = x(1).^2 + x(2).^2-1; 
	y2 = x(1).^2 + x(2).^2-9;
	out = [y1;y2];
