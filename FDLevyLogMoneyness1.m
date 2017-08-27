function [price] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
sigma = param(1); % diffusion volatility
lambda = param(2); % jump intensity
dt = T / M;
Smin = S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T));
Smax = S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T));
xmin = log(Smin / K); % logmoneyness transformation
xmax = log(Smax / K);
dx = (xmax - xmin) / N;
x = (xmin:dx:xmax)';
nodes = x(2:end-1);
%% compute alpha
k = LevyDensity(param, model);
[alpha,lambda_num,Bl,Br] = LevyIntegral1(k,N)
%% matrix definition
e = ones([N-1 1]);
switch NumMethod
	case 'explicit'
		A_tilde = -(-sigma^2/(2*dx^2) - (sigma^2/2 + alpha) / (2*dx));
		B_tilde = 1 / dt -(sigma^2 / dx^2 + lambda);
		C_tilde = -(-sigma^2/(2*dx^2) + (sigma^2/2 + alpha) / (2*dx));
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e], -1:1, N-1, N-1);
		M1 = spdiags(1/dt*e, 0, N-1, N-1);
	case 'implicit'
		A = -sigma^2/(2*dx^2) - (sigma^2/2 + alpha) / (2*dx);
		B = 1 / dt + sigma^2 / dx^2 + lambda;
		C = -sigma^2/(2*dx^2) + (sigma^2/2 + alpha) / (2*dx);
		M1 = spdiags([A*e B*e C*e], -1:1, N-1, N-1);
		M2 = spdiags(1/dt*e, 0, N-1, N-1);
	case 'CN'
		A = theta * (-sigma^2/(2*dx^2) - (sigma^2/2 + alpha) / (2*dx));
		B = 1 / dt + theta * (sigma^2 / dx^2 + lambda);
		C = theta * (-sigma^2/(2*dx^2) + (sigma^2/2 + alpha) / (2*dx));
		A_tilde = -(1 - theta) * (-sigma^2/(2*dx^2) - (sigma^2/2 + alpha) / (2*dx));
		B_tilde = 1 / dt -(1 - theta) * (sigma^2 / dx^2 + lambda);
		C_tilde = -(1 - theta) * (-sigma^2/(2*dx^2) + (sigma^2/2 + alpha) / (2*dx));
		M1 = spdiags([A*e B*e C*e], -1:1, N-1, N-1);
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e], -1:1, N-1, N-1);
	otherwise
		error('NumMethod invalid')
end
%% forward in time
transf = 'LogMoneyness';
BC = zeros([N-1 1]);
switch optionType
	case 'Call'
		u = subplus(exp(nodes) - 1);
		h = waitbar(0,'wait please...');
		for j = 1 : M
			BC(end) = -M1(1,2) * (exp(xmax) - 1) + ...
				M2(1,2) * (exp(xmax) - 1);
			if lambda == 0
				I = 0;
			else
				I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf);
			end
			u = M1 \(M2 * u + BC + I);
			waitbar(j/M, h);
		end
		close(h)
	case 'Put'
		u = subplus(1 - exp(nodes));
		h = waitbar(0,'wait please...');
		for j = 1 : M
			BC(1) = -M1(2,1) * (1 - exp(xmin)) + ...
				M2(2,1) * (1 - exp(xmin));
			if lambda == 0
				I = 0;
			else
				I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf);
			end
			u = M1 \(M2 * u + BC + I);
			waitbar(j/M, h);
		end
		close(h)
	otherwise
		error('optionType')
end
S = K * exp(nodes - r * T);
C = K * u * exp(-r * T);
price = interp1(S,C,S0,'spline');
end

