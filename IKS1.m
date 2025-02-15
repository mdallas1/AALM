function out = IKS1(x)
	% --------------------------------------------------------
	% Implmentation of Ex. 7 on page 487 of 
	%
	%	Izmailov A.F., Kurennoy A.S., Solodov M.V., Critical 
	%	solutions of nonlinear equations: stability issues, 
	%	Math. Program. Ser. B, 2018.
	%
	% Any nonzero solution to IKS1(x) = 0 is noncritical, but 
	% the zero solution is critical.
	% --------------------------------------------------------
	out1 = x(1)*x(2);
	out2 = x(1)*x(3);
	out3 = x(2)*x(3);
	out = [out1;out2;out3];
end
