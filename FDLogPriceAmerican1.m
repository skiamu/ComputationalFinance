function [price] = FDLogPriceAmerican1( S0,K,r,T,N,M,sigma,NumMethod,theta,...
	maxiter,tol,omega)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
%% grid
dt = T / M;
t = 0:dt:T;
xmin = log(S0 * exp((r - sigma^2/2) * T -6 * sigma * sqrt(T)));
xmax = log(S0 * exp((r - sigma^2/2) * T +6 * sigma * sqrt(T)));
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
V = subplus(K - exp(x(2:end-1)));
h = waitbar(0,'Please wait...');
for j = M : -1 : 1
	BC(1) = -M1(2,1) * (exp(-r*(T - t(j)))*K - exp(xmin)) + ...
		M2(2,1) * (exp(-r*(T - t(j+1)))*K - exp(xmin));
	% 	V = M1 \(M2 * V + BC);
	b = M2 * V + BC;
	V = PSORalgorithm1(M1,b,V,maxiter,tol,omega,K,N,x);
	waitbar((M - j)/M, h);
end
price = interp1(exp(x(2:end-1)),V,S0,'spline');
close(h)
end

