function out = dBeh2(x)
	% --------------------------------------------------------
	% Implmentation of Jacobian of function given on p. 7 of 
	%
	%	Behling R., Fischer A., Herrich M., Iusem A., and Ye Y.,
	%	*A Levenberg-Marquardt method with approximate projections*, 
	%	Comput. Optim. Appl., 2014.  
	%
	% local error bound holds on solution set, but Jacobian 
	% is rank deficient. 
	% --------------------------------------------------------

	out1 = [0 1];
	out2 = [x(2)*exp(x(1)) exp(x(1))];
	out = [out1;out2];
