function aalm(f,df,x,varargin)
	%{

	AALM(f,df,x)
	 implements Anderson accelerated Levenberg-Marquardt (AALM) 
	 with algorithmic depth 1. 
	 
	-- f and it's Jacobian df are function handles.
	-- x initial iterate.

	-- aalm(f,df,x,'anderson',m) implements Anderson acceleration with depth m.
	-- aalm(f,df,x,'safeguard',r) implements adaptive gamma safeguarding with tolerance r.
	-- aalm(f,df,x,'print') prints residual at each iteration. 
	-- aalm(f,df,x,'tol',tol) sets tolerance. Default tol = 1e-8. 
	-- aalm(f,df,x,'maxiters',maxiters) sets maximum iterations. Default maxiters = 100.

	-- Plot parameters may be passed with aalm(...,'plot',plot_options), where 
	plot_options is a cell array. 

	%}

	assert(isa(f,'function_handle'));
	assert(isa(df,'function_handle'));
	assert(isa(x,'double'));

	[x_m,x_n] = size(x); 
	if x_n > x_m 
		x = x';
	end

	p = inputParser(); 
	p.FunctionName = 'aalm' ; 

	default_tol = 1e-8;
	p.addParameter('tol',default_tol,@isdouble);

	default_maxiters = 100;
	p.addParameter('maxiters',default_maxiters,@isdouble);

	vld_solver = @(s) any(strcmp(s,{'newton','newt','lm','LM'}));
	default_solver = 'lm';
	p.addParameter('solver',default_solver,vld_solver);

	default_depth= 1;
	p.addParameter('anderson',default_depth,@isfloat);

	default_gsg_tol = 0; 
	p.addParameter('safeguard',default_gsg_tol,@isfloat);

	default_plot = {}; 
	p.addParameter('plot',default_plot,@(a) isa(a,'cell'));

	p.addSwitch('print');

	p.parse(varargin{:}); 
	solver = p.Results.solver; m = p.Results.anderson; gsg_tol = p.Results.safeguard;
	tol = p.Results.tol; maxiters = p.Results.maxiters; 
	plot_options = p.Results.plot; 

	if any(strcmp('anderson',varargin));
		aa = 1;
	else 
		aa = 0;
	end 
	if any(strcmp('safeguard',varargin));
		gsg = 1;
	else 
		gsg = 0;
	end

	f = @(x) f(x); df = @(x) df(x);
	x = x ; fx = f(x); dfx = df(x); 
	iters = 0; res = norm(fx); I = eye(length(x));

	res_arr = []; x_arr = []; d_arr = []; xaa_arr = [];

	if res < tol 
		disp("Initial iterate satisfies tolerance.");
	end

	% SOLVER
	while (iters < maxiters && res > tol)

		
		res_arr = [res res_arr];
		
		if iters < 1

			x_arr = [x x_arr]; 

			if any(strcmp(solver,{'lm','LM'}))
				mu = 0.5*1e-8*norm(fx)^2; % LM parameter (from KaYaFu03)
				d = - (dfx'*dfx + mu * I)\(dfx'*fx);
			elseif any(strcmp(solver,{'newt','newton'}))
				d = -dfx\fx;
			end
			%d_arr = [d d_arr(1:mm-1)];
			x_arr = [x x_arr];
			d_arr = [d d_arr];
			x = x+d;

		else 
			if any(strcmp(solver,{'lm','LM'}))
				mu = min(mu,nfx^2);						
				d = - (dfx'*dfx + mu * I)\(dfx'*fx);
			elseif any(strcmp(solver,{'newt','newton'}))
				d = -dfx\fx;
			end

			% ANDERSON	
			if aa 
				mm = min(m,iters);
				x_arr = [x x_arr(:,1:mm)];
				d_arr = [d d_arr(:,1:mm)];
				E = x_arr(:,1:end-1) - x_arr(:,2:end);
				F = d_arr(:,1:end-1) - d_arr(:,2:end);
				gamma = F(:,1:mm)\d_arr(:,1);
			end

			if gsg 
				[gamma,lam,r] = gamma_safeguard(d_arr(:,1),d_arr(:,2),gamma,gsg_tol,1,I);
			end

			if aa 
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
		if p.Results.print
			fprintf("iter \t res \n-----------------------------------\n%g \t %.3e\n",iters,res);	
		end
	end
	res_arr = [res res_arr];
	%semilogy([0:iters],flip(res_arr),'r-*','Linewidth',1.5);
	semilogy([0:iters],flip(res_arr),plot_options{:});
	x

end
			
