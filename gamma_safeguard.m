function [out1,out2,out3] = gamma_safeguard(w1,w0,gamma,tol,adap,M)
	% =================================================================
	% Safeguarding scheme for Newton-Anderson with Anderson depth 1
	% Author: Matt Dallas, University of Dallas Mathematics Dept
	% email: mdallas@udallas.edu
	%
	% Introduced in  
	%
	%   Dallas M., Pollock S., Newton-Anderson at Singular Points, 
	%   International Journal of Numerical Analysis and Modeling, 
	%   doi: 10.4208/ijnam2023-1029
	% 
	% Present form as well as adaptive version analyzed in 
	%
	%	Dallas M., Pollock S., Rebholz L. G., ANALYSIS OF AN ADAPTIVE SAFEGUARDED 	 % 	 NEWTON-ANDERSON ALGORITHM OF DEPTH ONE WITH APPLICATIONS TO FLUID 
	%	PROBLEMS. ACSE. 2024. doi:10.3934/acse.2024013
	%
	% -----------------------------------------------------------------
	% gamma_safeguard(w1,w0,gamma,tol,adap,M)
	%	- w1 and w0 are the previous two Newton residuals
	%	- gamma is the gamma value computed by Newton-Anderson
	%	- tol is a user chosen parameter between 0 and 1. 			 
	%	  The closer tol is to 0, the stricter the safegurading, 
	%     in which case the iterates behave more like standard 
	%     Newton.
	% 	- adap is set to 1 to run adaptive gamma-safeguarding, and set 
	%	  to 0 to run standard gamma-safeguarding		
	% 	- M is an SPD matrix used to compute a weighted norm 
	%  	  of the residuals. For example, M may be the stiffness 
	% 	  matrix if the norm used is the H1 seminorm. 
	% -----------------------------------------------------------------
	% Notes:
	% 	- A weakened form of safeguarding is obtained by replacing 
	% 	  "gamma>=1" with "cos(w_k+1,w_k)>0.942", or more precisely, 
	% 	  |\beta|<2*sqrt(2)/3. This allows gamma to reach 2. 
	% =================================================================
	normw1 = sqrt(w1'*M*w1);
	normw0 = sqrt(w0'*M*w0);
	if adap == 0
		r = tol;
	elseif adap == 1
		r = min(normw1/normw0,tol);
	else 
		error("adap must be 0 or 1");
	end
	beta = r*(normw1/normw0);
	lambda=1;
    if abs(gamma)<1e-14 || gamma>=1 
        lambda=0;
    elseif abs(gamma)/(1-gamma)>beta
		lambda = beta/(gamma*(beta+sign(gamma)));
    end
	out3 = r;
	out2 = lambda;
    out1 = lambda*gamma;
end



