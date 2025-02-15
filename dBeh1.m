function out = Beh1(x)
	% --------------------------------------------------------
	% Implmentation Jacobian of Ex. 5.2 on page 1118 of 
	%
	%	Behling R., Goncalves, D.S., and Santos, S.A., 
	%	*Local Convergence Analysis of the Levenberg–Marquardt 
	%	Framework for Nonzero-Residue Nonlinear Least-Squares Problems 
	%	Under an Error Bound Condition*, J. Op. and App., 2019.  
	%
	% Has a set of non-isolated minimizers, and a global 
	% minimizer (-1,0)^T. The local error bound holds, and 
	% the rank varies at the axis $x_1=0$. 
	%
	% Does local error hold at origin? 2-regular?
	% --------------------------------------------------------
	% 	
	y1 = [3*x(1)^2 - x(2), -x(1)];
	y2 = [3*x(1)^2+x(2),x(1)];
	out = [y1;y2];
