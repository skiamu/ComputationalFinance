function [price] = FDLogPrice1( S0,K,r,T,N,M,sigma,optionType,NumMethod,theta,...
	barrierType,barrier)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
dt = T / M;
t = 0:dt:T;
xmin = log(S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T)));
xmax = log(S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T)));
if nargin > 10
	switch barrierType
		case {'DO','UI'}
			xmin = log(barrier);
		case {'UO','DI'}
			xmax = log(barrier);
		otherwise
	end
end
dx = (xmax - xmin) / N;
x = (xmin:dx:xmax)';
%% matrix definition
e = ones([N-1 1]);
switch NumMethod
	case 'explicit'
		alpha_tilde = (-(r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		beta_tilde = 1 / dt + (-sigma^2/dx^2 - r);
		gamma_tilde = ((r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		M2 = spdiags([alpha_tilde*e beta_tilde*e gamma_tilde*e],-1:1,N-1,N-1);
		M1 = spdiags((1/dt)*e,0,N-1,N-1);
	case 'implicit'
		alpha = -(-(r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		beta = 1 / dt - (-sigma^2/dx^2 - r);
		gamma = -((r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		M1 = spdiags([alpha*e beta*e gamma*e],-1:1,N-1,N-1);
		M2 = spdiags((1/dt)*e,0,N-1,N-1);
	case 'CN'
		alpha = -theta * (-(r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		beta = 1 / dt - theta * (-sigma^2/dx^2 - r);
		gamma = -theta * ((r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		alpha_tilde = (1 - theta) * (-(r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		beta_tilde = 1 / dt + (1 - theta) * (-sigma^2/dx^2 - r);
		gamma_tilde = (1 - theta) * ((r - sigma^2/2) / (2*dx) + sigma^2/(2*dx^2));
		M1 = spdiags([alpha*e beta*e gamma*e],-1:1,N-1,N-1);
		M2 = spdiags([alpha_tilde*e beta_tilde*e gamma_tilde*e],-1:1,N-1,N-1);
	otherwise
		error('NumMethod invalid')
end
%% backward in time
BC = zeros([N-1 1]);
if nargin == 10 % european case
	switch optionType
		case 'Call'
			V = subplus(exp(x(2:end-1)) - K);
			for j = M : -1 : 1
				BC(end) = -M1(1,2) * (exp(xmax) - exp(-r*(T - t(j)))*K) + ...
					M2(1,2) * (exp(xmax) - exp(-r*(T - t(j+1)))*K);
				V = M1 \(M2 * V + BC);
			end
		case 'Put'
			V = subplus(K - exp(x(2:end-1)));
			for j = M : -1 : 1
				BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - exp(xmin)) + ...
					M2(2,1) * (exp(-r*(T - t(j+1)))*K - exp(xmin));
				V = M1 \(M2 * V + BC);
			end
		otherwise
	end
elseif nargin > 10 % barrier case
	switch optionType
		case 'Call'
			switch barrierType
				case 'DO'
					V = subplus(exp(x(2:end-1)) - K);
					for j = M : -1 : 1
						BC(end) = -M1(1,2) * (exp(xmax) - exp(-r*(T - t(j)))*K) + ...
							M2(1,2) * (exp(xmax) - exp(-r*(T - t(j+1)))*K);
						V = M1 \(M2 * V + BC);
					end
				case 'UO'
					V = subplus(exp(x(2:end-1)) - K);
					for j = M : -1 : 1
						V = M1 \(M2 * V + BC);
					end
			end
		case 'Put'
			switch barrierType
				case 'DO'
					V = subplus(K - exp(x(2:end-1)));
					for j = M : -1 : 1
						V = M1 \(M2 * V + BC);
					end
				case 'UO'
					V = subplus(K - exp(x(2:end-1)));
					for j = M : -1 : 1
						BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - exp(xmin)) + ...
							M2(2,1) * (exp(-r*(T - t(j+1)))*K - exp(xmin));
						V = M1 \(M2 * V + BC);
					end
			end
		otherwise
	end
end
price = interp1(exp(x(2:end-1)),V,S0,'spline');
end

