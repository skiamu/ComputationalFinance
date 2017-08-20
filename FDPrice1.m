function [price] = FDPrice1( S0,K,r,T,N,M,sigma,optionType,NumMethod,theta,...
	barrierType,barrier)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
dt = T / M;
t = 0:dt:T;
Smin = S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T));
Smax = S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T));
if nargin > 10
	switch barrierType
		case {'DO','UI'}
			Smin = barrier;
		case {'UO','DI'}
			Smax = barrier;
		otherwise
	end
end
dS = (Smax - Smin) / N;
S = (Smin:dS:Smax)';
%% matrix definition
e = ones([N-1 1]);
switch NumMethod
	case 'explicit'
		alpha_tilde = (1 - theta) * (-(r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		beta_tilde = 1 / dt + (1 - theta) * (-sigma^2 * S(2:end-1).^2/dS^2 - r);
		gamma_tilde = (1 - theta) * ((r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		M2 = spdiags([alpha_tilde beta_tilde gamma_tilde],[1 0 -1],N-1,N-1)';
		M1 = spdiags((1/dt)*e,0,N-1,N-1);
	case 'implicit'
		alpha = -(-(r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		beta = 1 / dt -(-sigma^2 * S(2:end-1).^2/dS^2 - r);
		gamma = -((r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		M1 = spdiags([alpha beta gamma],[1 0 -1],N-1,N-1)';
		M2 = spdiags((1/dt)*e,0,N-1,N-1);
	case 'CN'
		alpha = -theta * (-(r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		beta = 1 / dt - theta * (-sigma^2 * S(2:end-1).^2/dS^2 - r);
		gamma = -theta * ((r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		alpha_tilde = (1 - theta) * (-(r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		beta_tilde = 1 / dt + (1 - theta) * (-sigma^2 * S(2:end-1).^2/dS^2 - r);
		gamma_tilde = (1 - theta) * ((r * S(2:end-1)) / (2*dS) + sigma^2 * S(2:end-1).^2/(2*dS^2));
		M1 = spdiags([alpha beta gamma],[1 0 -1],N-1,N-1)';
		M2 = spdiags([alpha_tilde beta_tilde gamma_tilde],[1 0 -1],N-1,N-1)';
	otherwise
		error('NumMethod invalid')
end
%% backward in time
BC = zeros([N-1 1]);
if nargin == 10 % european case
	switch optionType
		case 'Call'
			V = subplus(S(2:end-1) - K);
			for j = M : -1 : 1
				BC(end) = -M1(1,2) * (Smax - exp(-r*(T - t(j)))*K) + ...
					M2(1,2) * (Smax - exp(-r*(T - t(j+1)))*K);
				V = M1 \(M2 * V + BC);
			end
		case 'Put'
			V = subplus(K - S(2:end-1));
			for j = M : -1 : 1
				BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - Smin) + ...
					M2(2,1) * (exp(-r*(T - t(j+1)))*K - Smin);
				V = M1 \(M2 * V + BC);
			end
		otherwise
	end
elseif nargin > 10 % barrier case
	switch optionType
		case 'Call'
			switch barrierType
				case 'DO'
					V = subplus(S(2:end-1) - K);
					for j = M : -1 : 1
						BC(end) = -M1(1,2) * (Smax - exp(-r*(T - t(j)))*K) + ...
							M2(1,2) * (Smax - exp(-r*(T - t(j+1)))*K);
						V = M1 \(M2 * V + BC);
					end
				case 'UO'
					V = subplus(S(2:end-1) - K);
					for j = M : -1 : 1
						V = M1 \(M2 * V + BC);
					end
			end
		case 'Put'
			switch barrierType
				case 'DO'
					V = subplus(K - S(2:end-1));
					for j = M : -1 : 1
						V = M1 \(M2 * V + BC);
					end
				case 'UO'
					V = subplus(K - S(2:end-1));
					for j = M : -1 : 1
						BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - Smin) + ...
							M2(2,1) * (exp(-r*(T - t(j+1)))*K - Smin);
						V = M1 \(M2 * V + BC);
					end
			end
		otherwise
	end
end
price = interp1(S(2:end-1),V,S0,'spline');
end

