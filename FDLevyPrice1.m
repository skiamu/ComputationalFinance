function [price] = FDLevyPrice1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
	barrierType,barrier)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
sigma = param(1); % diffusion volatility
lambda = param(2); % jump intensity

dt = T / M;
t = 0:dt:T;
Smin = S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T));
Smax = S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T));
if nargin > 11
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
S = S(2:end-1);
nodes = S;
%% compute alpha
k = LevyDensity(param, model);
[alpha,lambda_num,Bl,Br] = LevyIntegral1(k,N)
%% matrix definition
e = ones([N-1 1]);
switch NumMethod
	case 'explicit'
		A_tilde = (0.5*sigma^2*S.^2/dS^2 - (r-alpha)*S/(2*dS));
		B_tilde = 1/dt + (-sigma^2*S.^2/dS^2 - (r + lambda));
		C_tilde = (0.5*sigma^2*S.^2/dS^2 + (r-alpha)*S/(2*dS));
		M2 = spdiags([A_tilde B_tilde C_tilde],[1 0 -1], N-1, N-1)';
		M1 = spdiags(1/dt*e,0,N-1,N-1);
	case 'implicit'
		A = -(0.5*sigma^2*S.^2/dS^2 - (r-alpha)*S/(2*dS));
		B = 1/dt -(-sigma^2*S.^2/dS^2 - (r + lambda));
		C = -(0.5*sigma^2*S.^2/dS^2 + (r-alpha)*S/(2*dS));
		M1 = spdiags([A B C],[1 0 -1], N-1, N-1)';
		M2 = spdiags(1/dt*e,0,N-1,N-1);
	case 'CN'
		A = -theta * (0.5*sigma^2*S.^2/dS^2 - (r-alpha)*S/(2*dS));
		B = 1/dt - theta * (-sigma^2*S.^2/dS^2 - (r + lambda));
		C = -theta * (0.5*sigma^2*S.^2/dS^2 + (r-alpha)*S/(2*dS));
		A_tilde = (1-alpha) * (0.5*sigma^2*S.^2/dS^2 - (r-alpha)*S/(2*dS));
		B_tilde = 1/dt + (1-alpha) * (-sigma^2*S.^2/dS^2 - (r + lambda));
		C_tilde = (1-alpha) * (0.5*sigma^2*S.^2/dS^2 + (r-alpha)*S/(2*dS));
		M1 = spdiags([A B C],[1 0 -1], N-1, N-1)';
		M2 = spdiags([A_tilde B_tilde C_tilde],[1 0 -1], N-1, N-1)';
	otherwise
		error('NumMethod invalid')
end
%% backward in time
transf = 'Price';
BC = zeros([N-1 1]);
if nargin == 11 % european case
	switch optionType
		case 'Call'
			u = subplus(nodes - K);
			h = waitbar(0,'wait please...');
			for j = M : -1 : 1
				BC(end) = -M1(1,2) * (Smax - exp(-r*(T - t(j)))*K) + ...
					M2(1,2) * (Smax - exp(-r*(T - t(j+1)))*K);
				if lambda == 0
					I = 0;
				else
					I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf);
				end
				u = M1 \(M2 * u + BC + I);
				waitbar((M-j)/M, h);
			end
			close(h)
		case 'Put'
			u = subplus(K - nodes);
			h = waitbar(0,'wait please...');
			for j = M : -1 : 1
				BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - Smin) + ...
					M2(2,1) * (exp(-r*(T - t(j+1)))*K - Smin);
				if lambda == 0
					I = 0;
				else
					I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf);
				end
				u = M1 \(M2 * u + BC + I);
				waitbar((M-j)/M, h);
			end
			close(h)
		otherwise
	end
elseif nargin > 11 % barrier case
	switch optionType
		case 'Call'
			switch barrierType
				case 'DO'
					u = subplus(nodes - K);
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						BC(end) = -M1(1,2) * (Smax - exp(-r*(T - t(j)))*K) + ...
							M2(1,2) * (Smax - exp(-r*(T - t(j+1)))*K);
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
				case 'UO'
					u = subplus(nodes - K);
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
			end
		case 'Put'
			switch barrierType
				case 'DO'
					u = subplus(K - nodes);
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
				case 'UO'
					u = subplus(K - nodes);
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - Smin) + ...
							M2(2,1) * (exp(-r*(T - t(j+1)))*K - Smin);
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
			end
		otherwise
	end
end
S = nodes;
C = u;
price = interp1(S,C,S0,'spline');
end


function I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
qnodes = linspace(Bl,Br,N+1);
dq = qnodes(2) - qnodes(1);
I = zeros([N-1 1]);
if nargin == 10 % european option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) .* exp(qnodes),nodes,u,K,S0,optionType,transf) .* k(qnodes));
	end
elseif nargin > 10 % barrier option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) .* exp(qnodes),nodes,u,K,S0,optionType,transf,barrierType) .* k(qnodes));
	end
end

end % f
