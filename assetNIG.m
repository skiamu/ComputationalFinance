function [ S, check ] = assetNIG( S0,r,T,param,Nsim,Nstep )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
sigma = param(1);
kappa = param(2);
Psi_X = @(v) 1/kappa - 1/kappa*sqrt(1 + sigma^2*kappa*v.^2);
theta = -Psi_X(-1i);
dt = T / Nstep;
dX = zeros(Nsim, Nstep+1);
Z = randn([Nsim Nstep]);
mu = dt;  lambda = dt^2/kappa; % gamma parameters
dS = simulateIG(mu,lambda,Nsim,Nstep);
dX(:,2:end) = r * dt + theta * dS + sigma * sqrt(dS) .* Z;
X = cumsum(dX, 2);
S = S0 .* exp(X);
check = mean(S(:,end) / exp(r * T));

	function IG = simulateIG(mu, lambda, N, M)
		method = 1;
		if method == 1
		Y = randn([N M]).^2;
		IG = mu + mu^2 * Y / (2*lambda) - mu / (2*lambda) .* ...
			(sqrt(4*mu*lambda*Y + mu^2*Y.^2));
		U = rand([N M]);
		IG(U > mu./(IG + mu)) = mu^2 ./ IG(U > mu./(IG + mu));
		else
			IG = icdf('InverseGaussian',rand([N M]),mu, lambda);
		end
	end
end

