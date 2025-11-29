function out = Beh3(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 5.3 on page 1119 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	% J. Op. and App., 2019.  
	%
	% Has a set of non-isolated minimizers x(2) = 0.
	% The local error bound holds, and 
	% the rank decreases as x approaches X^*.
	% LM converges linearly since mu_k bounded away from zero. 
	%
	% x0 = [pi;0.001]; 
	% mu = 0.2; 
	% --------------------------------------------------------
	y1 = (1/9)*cos(x(1))-x(2).*sin(x(1));
	y2 = (1/9)*sin(x(1))+x(2).*cos(x(1));
	out = [y1;y2];
