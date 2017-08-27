function [price] = FDLevyLogpriceFwdprice1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
	epsilon,barrierType,barrier)
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
xmin = log(Smin); % logmoneyness transformation
xmax = log(Smax);
if nargin > 12
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
		A_tilde = (1 - theta) * (0.5 * sigma^2 / dx^2 - (r - 0.5 * sigma^2 - alpha) / (2*dx));
		B_tilde = 1 / dt - (1 - theta) * (-sigma^2/dx^2 - lambda);
		C_tilde = (1 - theta) * (0.5 * sigma^2 / dx^2 + (r - 0.5 * sigma^2 - alpha) / (2*dx));
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e],-1:1,N-1,N-1);
		M1 = spdiags(1/dt*e,0,N-1,N-1);
	case 'implicit'
		A = -(0.5 * sigma^2 / dx^2 - (r - 0.5 * sigma^2 - alpha) / (2*dx));
		B = 1 / dt -(-sigma^2/dx^2 - lambda);
		C = -(0.5 * sigma^2 / dx^2 + (r - 0.5 * sigma^2 - alpha) / (2*dx));
		M1 = spdiags([A*e B*e C*e],-1:1,N-1,N-1);
		M2 = spdiags(1/dt*e,0,N-1,N-1);
	case 'CN'
		A = -theta * (0.5 * sigma^2 / dx^2 - (r - 0.5 * sigma^2 - alpha) / (2*dx));
		B = 1 / dt - theta * (-sigma^2/dx^2 - lambda);
		C = -theta * (0.5 * sigma^2 / dx^2 + (r - 0.5 * sigma^2 - alpha) / (2*dx));
		A_tilde = (1 - theta) * (0.5 * sigma^2 / dx^2 - (r - 0.5 * sigma^2 - alpha) / (2*dx));
		B_tilde = 1 / dt - (1 - theta) * (-sigma^2/dx^2 - lambda);
		C_tilde = (1 - theta) * (0.5 * sigma^2 / dx^2 + (r - 0.5 * sigma^2 - alpha) / (2*dx));
		M1 = spdiags([A*e B*e C*e],-1:1,N-1,N-1);
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e],-1:1,N-1,N-1);
	otherwise
		error('NumMethod invalid')
end
%% forward in time
transf = 'LogPriceFwdPrice';
BC = zeros([N-1 1]);
if nargin == 12 % european case
	switch optionType
		case 'Call'
			u = subplus(exp(nodes) - K);
			h = waitbar(0,'wait please...');
			for j = M : -1: 1
				BC(end) = -M1(1,2) * (Smax*exp(r*(T-j*dt))-K) + ...
					M2(1,2) * (Smax*exp(r*(T-j*dt))-K);
				if lambda == 0
					I = 0;
				else
					I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf);
				end
				u = M1 \(M2 * u + BC + I);
				waitbar((M-j)/M, h);
			end
			close(h)
		case 'Put'
			u = subplus(K - exp(nodes));
			h = waitbar(0,'wait please...');
			for j = M : -1: 1
				BC(1) = -M1(2,1) * (K-Smin*exp(r*(T-j*dt))) + ...
					M2(2,1) * (K-Smin*exp(r*(T-j*dt)));
				if lambda == 0
					I = 0;
				else
					I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf);
				end
				u = M1 \(M2 * u + BC + I);
				waitbar((M-j)/M, h);
			end
			close(h)
		otherwise
	end
elseif nargin > 12 % barrier case
	switch optionType
		case 'Call'
			switch barrierType
				case 'DO'
					u = subplus(exp(nodes) - K);
					h = waitbar(0,'wait please...');
					for j = M : -1: 1
						BC(end) = -M1(1,2) * (Smax*exp(r*(T-j*dt))-K) + ...
							M2(1,2) * (Smax*exp(r*(T-j*dt))-K);
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
				case 'UO'
					u = subplus(exp(nodes) - K);
					h = waitbar(0,'wait please...');
					for j = M : -1: 1
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
			end
		case 'Put'
			switch barrierType
				case 'DO'
					u = subplus(K - exp(nodes));
					h = waitbar(0,'wait please...');
					for j = M : -1: 1
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
				case 'UO'
					u = subplus(K - exp(nodes));
					h = waitbar(0,'wait please...');
					for j = M : -1: 1
						BC(1) = -M1(2,1) * (K-Smin*exp(r*(T-j*dt))) + ...
							M2(2,1) * (K-Smin*exp(r*(T-j*dt)));
						if lambda == 0
							I = 0;
						else
							I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf,barrierType);
						end
						u = M1 \(M2 * u + BC + I);
						waitbar((M-j)/M, h);
					end
					close(h)
			end
		otherwise
	end
end
S = exp(nodes);
C = u * exp(-r * T);
price = interp1(S,C,S0,'spline');
end

