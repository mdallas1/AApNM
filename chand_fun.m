function out = chand_fun(x,p,i)
	% ===========================================
	% IMPLEMENTATION OF CHANDRASEKHAR H-EQUATION 
	% AND IT'S JACOBIAN. 
	%
	% - CHAN_FUN(U,P,I)
	% 	- X IS INPUT VARIABLE, P IS THE 
	%			PARAMETER TYPICALLY DENOTED \OMEGA. 
	% 	- WHEN P\IN(0,1), THE PROBLEM IS 
	%			NONSINGULAR, AND WHEN P=1 THE 
	%			PROBLEM IS SINGULAR. 
	%  	- I == 0 THE FUNCTION IS RETURNED, IF 
	%			I == 1 THEN THE JACOBIAN IS RETURNED. 
	%
	% NOTES:
	%
	% 	- THE NONTRIVIAL ROOT IS A FIRST ORDER 
	%			SINGULARITY AND THE JACOBIAN AT THE 
	% 	  SOLUTION HAS A ONE-DIMENSIONAL NULL 
	%			SPACE. FOR MORE INFORMATION CONCERNING 
	%			THE CHANDRASEKHAR H-EQUATION, SEE 
	%			P.229FF OF 
	%
	%			KELLEY, C. (2018). NUMERICAL METHODS 
	%			FOR NONLINEAR EQUATIONS. 
	%			ACTA NUMERICA, 27, 207-287. 
	%			DOI:10.1017/S0962492917000113
	%
	% 		and references therein. 
	% ===========================================
	N=length(x);
	I=eye(N);
	%---------------------
	% FORM DIAGONAL MATRIX
	%---------------------
	v = [1:N]';
	v = v-0.5*ones(N,1);
	D = (p/(2*N))*diag(v);
	%---------------------
	% FORM HANKEL MATRIX
	%---------------------
	w = [2:2*N]';
	w = w - ones(2*N-1,1);
	w = 1./w;
	a = w(1:N);
	b = w(N:2*N-1);
	H = hankel(a,b);
	L = D*H;
	u = 1./(1-L*x);
	if i==0
		%----------
		% DEFINE F
		%----------
		y = zeros(N,1);
		y = x-u;
		out = y;
	else
		%-----------
		% DEFINE DF
		%-----------
		z = zeros(N,1);
		D_hat = diag(u.^2);
		z = I-D_hat*L;
		out = z;
	end
end
	
	
	


