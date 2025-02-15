function out = chand_fun(x,p,i)
	% ==============================================================================
	% Implementation of Chandrasekhar H-equation and it's Jacobian 
	%
	% ------------------------------------------------------------------------------
	% chan_fun(x,p,i)
	% 	- x is input variable, p is the parameter typically denoted \omega. 
	% 	- When p\in(0,1), the problem is nonsingular, and when p=1 the problem is 
	% 	  singular. 
	%  	- i == 0 the function is returned, if i == 1 then the jacobian is returned. 
	% ------------------------------------------------------------------------------
	%
	% Notes:
	% 	- The nontrivial root is a first order singularity and the Jacobian at the 
	% 	  solution has a one-dimensional null space. 
	% 	- For more information concerning the Chandrasekhar H-equation, see p.229ff of 
	%
	%		Kelley, C. (2018). Numerical methods for nonlinear equations. 
	%		Acta Numerica, 27, 207-287. doi:10.1017/S0962492917000113
	%
	% 	  and references therein. 
	% ==============================================================================

	N=length(x);
	I=eye(N);
	%---------------------
	% Form diagonal matrix
	%--------------------
	v = [1:N]';
	v = v-0.5*ones(N,1);
	D = (p/(2*N))*diag(v);
	%---------------------
	% Form Hankel matrix
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
		%--------------------
		% Define F
		%--------------------
		y = zeros(N,1);
		y = x-u;
		out = y;
	else
		%--------------------
		% Define dF
		%--------------------
		z = zeros(N,1);
		D_hat = diag(u.^2);
		z = I-D_hat*L;
		out = z;
	end
end
	
	
	


