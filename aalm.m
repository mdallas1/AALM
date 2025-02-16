function [out1,out2] = aalm(f,df,x,varargin)
	%{

	AALM(f,df,x)
	 implements Anderson accelerated Levenberg-Marquardt (AALM) 
	 with algorithmic depth 1. 
	 
	-- f and it's Jacobian df are function handles.
	-- x initial iterate.

	-- optional args {'solver' (char), 'tol' (float), and 'maxiters' (float)}
	 solvers = {'newton','newt','lm','LM'}.

	-- aalm(f,df,x,'aa') implements Anderson acceleration. 
	-- aalm(f,df,x,'aa',m) sets depth to m.
	-- aalm(f,df,x,'gsg',r) implements adaptive gamma safeguarding with tolerance r.
	-- aalm(f,df,x,'print') prints residual at each iteration. 

	%}

	solvers = {'newton','newt','lm','LM'};

	f = @(x) f(x); df = @(x) df(x);
	x = x ; fx = f(x); dfx = df(x); 
	iters = 0; res = norm(fx); I = eye(length(x));

	if nargin < 4	
		tol = 1e-8; maxiters = 100; 
		solver = 'lm'; aa = 1; m = 1; gsg = 0;
	end
	
	if nargin >= 4 
		solver_index = find(strcmp('solver',varargin));
		tol_index = find(strcmp('tol',varargin));
		maxiters_index = find(strcmp('maxiters',varargin));
		aa_index = find(strcmp('aa',varargin));

		if !isempty(solver_index)
			solver = varargin{solver_index(1)+1};
			if isempty(find(strcmp(solver,solvers)))
				error('Solver not in solver list. Type "help AALM."')
			end
		end

		if !isempty(tol_index)
			tol = varargin{tol_index(1)+1};
		else 
			tol = 1e-8;
		end
	
		if !isempty(maxiters_index) 
			maxiters = varargin{maxiters_index(1)+1};
		else 
			maxiters = 100;
		end
	
	if !isempty(find(strcmp('aa',varargin)))
		aa = 1;	
		aaindex = find(strcmp('aa',varargin))(1);		
		try 	
			if isfloat(varargin{aaindex+1})
				m = varargin{aaindex+1};
			end
		catch 
			disp('!! ERROR: "aa" option requires algorithmic depth.');	
			return
		end
	else 
		aa = 0;
	end

	if !isempty(find(strcmp('gsg',varargin)))
		gsg = 1; 
		gsgindex = find(strcmp('gsg',varargin))(1);		
		try 	
			if isfloat(varargin{gsgindex+1})
				gamma_tol = varargin{gsgindex+1};
			end
		catch 
			disp('!! ERROR: "gsg" option requires safeguarding tolerance.');	
			return
		end
	else 
		gsg = 0;
	end

	end

	res_arr = []; x_arr = []; d_arr = []; xaa_arr = [];

	if res < tol 
		disp("Initial iterate satisfies tolerance.");
	end

	% SOLVER
	while (iters < maxiters && res > tol)

		
		res_arr = [res res_arr];
		
		if iters < 1

			x_arr = [x x_arr]; 

			if strcmp(solver,'lm') || strcmp(solver,'LM')  
				mu = 0.5*1e-8*norm(fx)^2; % LM parameter (from KaYaFu03)
				d = - (dfx'*dfx + mu * I)\(dfx'*fx);
			elseif strcmp(solver,'newt') || strcmp(solver,'newton')
				d = -dfx\fx;
			end
			%d_arr = [d d_arr(1:mm-1)];
			x_arr = [x x_arr];
			d_arr = [d d_arr];
			x = x+d;

		else 
			if strcmp(solver,'lm') || strcmp(solver,'LM')
				mu = min(mu,nfx^2);						
				d = - (dfx'*dfx + mu * I)\(dfx'*fx);
			elseif strcmp(solver,'newton') || strcmp(solver,'newt')
				d = -dfx\fx;
			end

			% ANDERSON	
			if aa 
				mm = min(m,iters);
				x_arr = [x x_arr(:,1:mm)];
				d_arr = [d d_arr(:,1:mm)];
				%E = x_arr(:,1:end-1) - x_arr(:,2:end);
				%F = d_arr(:,1:end-1) - d_arr(:,2:end);
				E = x_arr(:,1:end-1) - x_arr(:,2:end);
				F = d_arr(:,1:end-1) - d_arr(:,2:end);
				gamma = F(:,1:mm)\d_arr(:,1);
			end

			if gsg 
				%[sg_gamma,lam,r] = gamma_safgrd(d_arr(:,1),d_arr(:,2),gamma,gamma_tol,1,eye(length(d)));
				[gamma,lam,r] = gamma_safeguard(d_arr(:,1),d_arr(:,2),gamma,gamma_tol,1,eye(length(d)));
			end

			if aa 
				 %x = x + d - (E(:,1)+F(:,1))*gamma;
			     x = x + d - (E(:,1:mm)+F(:,1:mm))*gamma;
			else 
				x = x+d;
			end
			
		end
		iters = iters+1;
		fx = f(x);
		nfx = norm(fx);
		dfx = df(x); 
		res = nfx; 
		%if iters > 1
		%	log(res_arr(:,1))/log(res_arr(:,2))
		%end
		if !isempty(find(strcmp('print',varargin)))
			fprintf("iter \t res \n-----------------------------------\n%g \t %.3e\n",iters,res);	
		end
	end
	res_arr = [res res_arr];
	semilogy([0:iters],flip(res_arr),'r-*','Linewidth',1.5);

end
			
