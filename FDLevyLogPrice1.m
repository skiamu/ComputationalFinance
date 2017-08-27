function [price] = FDLevyLogPrice1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
	barrierType,barrier)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
sigma = param(1); % diffusion volatility
lambda = param(2); % jump intensity
dt = T / M;
Smin = S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T));
Smax = S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T));
xmin = log(Smin / S0); % logmoneyness transformation
xmax = log(Smax / S0);
if nargin > 11
	switch barrierType
		case {'DO','UI'}
			xmin = log(barrier / S0);
		case {'UO','DI'}
			xmax = log(barrier / S0);
		otherwise
	end
end
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
		A_tilde = (0.5*sigma^2/dx^2 - (r - 0.5*sigma^2 - alpha)/(2*dx));
		B_tilde = 1/dt + (-sigma^2/dx^2 - (lambda + r));
		C_tilde = (0.5*sigma^2/dx^2 + (r - 0.5*sigma^2 - alpha)/(2*dx));
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e], -1:1, N-1, N-1);
		M1 = spdiags(1/dt*e, 0, N-1, N-1);
	case 'implicit'
		A = -(0.5*sigma^2/dx^2 - (r - 0.5*sigma^2 - alpha)/(2*dx));
		B = 1/dt - (-sigma^2/dx^2 - (lambda + r));
		C = -(0.5*sigma^2/dx^2 + (r - 0.5*sigma^2 - alpha)/(2*dx));
		M1 = spdiags([A*e B*e C*e], -1:1, N-1, N-1);
		M2 = spdiags(1/dt*e, 0, N-1, N-1);
	case 'CN'
		A = -theta * (0.5*sigma^2/dx^2 - (r - 0.5*sigma^2 - alpha)/(2*dx));
		B = 1/dt - theta * (-sigma^2/dx^2 - (lambda + r));
		C = -theta * (0.5*sigma^2/dx^2 + (r - 0.5*sigma^2 - alpha)/(2*dx));
		A_tilde = (1 -theta) * (0.5*sigma^2/dx^2 - (r - 0.5*sigma^2 - alpha)/(2*dx));
		B_tilde = 1/dt + (1 -theta) * (-sigma^2/dx^2 - (lambda + r));
		C_tilde = (1 -theta) * (0.5*sigma^2/dx^2 + (r - 0.5*sigma^2 - alpha)/(2*dx));
		M1 = spdiags([A*e B*e C*e], -1:1, N-1, N-1);
		M2 = spdiags([A_tilde*e B_tilde*e C_tilde*e], -1:1, N-1, N-1);
	otherwise
		error('NumMethod invalid')
end
%% backward in time
transf = 'LogPrice';
BC = zeros([N-1 1]);
if nargin == 11 % european case
	switch optionType
		case 'Call'
			u = subplus(S0 * exp(nodes) - K);
			h = waitbar(0,'wait please...');
			for j = M : -1 : 1
				BC(end) = -M1(1,2) * (S0*exp(xmax) - K) + ...
					M2(1,2) * (S0*exp(xmax) - K);
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
			u = subplus(K - S0 * exp(nodes));
			h = waitbar(0,'wait please...');
			for j = M : -1 : 1
				BC(1) = -M1(2,1) * (K - S0 * exp(xmin)) + ...
					M2(2,1) * (K - S0 * exp(xmin));
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
					u = subplus(S0 * exp(nodes) - K);
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						BC(end) = -M1(1,2) * (S0*exp(xmax) - K) + ...
							M2(1,2) * (S0*exp(xmax) - K);
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
					u = subplus(S0 * exp(nodes) - K);
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
					u = subplus(K - S0 * exp(nodes));
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
					u = subplus(K - S0 * exp(nodes));
					h = waitbar(0,'wait please...');
					for j = M : -1 : 1
						BC(1) = -M1(2,1) * (K - S0 * exp(xmin)) + ...
					M2(2,1) * (K - S0 * exp(xmin));
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
S = S0 * exp(nodes);
C = u;
price = interp1(S,C,S0,'spline');
end

