function [price] = FDLevyLogMoneyness1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
	epsilon)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
if length(param) == 4
	sigmaDiff = param(4);
else
	sigmaDiff = 0;
end
dt = T / M;
% Smin = S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T));
% Smax = S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T));
Smin = 10;
Smax = 300;
xmin = log(Smin / K); % logmoneyness transformation
xmax = log(Smax / K);
dx = (xmax - xmin) / N;
x = (xmin:dx:xmax)';
nodes = x(2:end-1);
%% compute alpha
k = LevyDensity(param, model);
[alpha,lambda,sigma2_eps,Bl,Br] = LevyIntegral1AI(k,N,epsilon)
if lambda * dt < 1
	disp('sufficient condition for convergence satisfied')
else
	disp('sufficient condition for convergence NOT satisfied')
end
sigma = sqrt(sigma2_eps + sigmaDiff^2);
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
				I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf);
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
				I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf);
			end
			u = M1 \(M2 * u + BC + I);
			waitbar(j/M, h);
		end
		close(h)
	otherwise
end
S = K * exp(nodes - r * T);
C = K * u * exp(-r * T);
price = interp1(S,C,S0,'spline');

end % function
