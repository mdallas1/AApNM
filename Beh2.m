function out = Beh2(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 5.2 on page 1118 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares 
	% Problems Under an Error Bound Condition*, 
	%	J. Op. and App., 2019.  
	%
	% Isolated Global minimizer at [-1;0]. 
	% Local error bound holds. 
	% Nonisolated set of local minimzers 
	% when x_1 = 0. Rank of Jacobian varies. 
	%
	% x0 = [0.008;2.0];
	% --------------------------------------------------------
	y1 = x(1).^3 - x(1).*x(2)+1; 
	y2 = x(1).^3 + x(1).*x(2)+1;
	out = [y1;y2];
