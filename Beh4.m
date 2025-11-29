function out = Beh4(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 5.3 on page 1119 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	% J. Op. and App., 2019.  
	%
	% Unique minimizer at origin. 
	% The local error bound holds, and the rank varies.
	% Theory presented in paper does not guarantee convergence, 
	% and reported experiments fail except when mu is taking to 
	% be sufficiently large. Linear convergence is observed. 
	%
	% x0 = [0.01;0]; 
	% mu = 5 gives linear convergence
	% mu = res fails to converge, but convergence 
	% is recovered with AA using depth m > 1.
	% --------------------------------------------------------
	y1 = x(2)-x(1).^2-1;
	y2 = x(2)+x(1).^2+1;
	out = [y1;y2];

