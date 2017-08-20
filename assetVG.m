function [ S, check ] = assetVG( S0,r,T,param,Nsim,Nstep )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
sigma = param(1); 
kappa = param(2);
Psi_X = @(v) -1/kappa*log(1 + kappa*sigma^2/2*v.^2);
theta = -Psi_X(-1i); % drift risk neutral
dt = T / Nstep;
dX = zeros([Nsim Nstep+1]);
dS = kappa * simulateGamma(dt / kappa, Nsim, Nstep);
Z = randn([Nsim Nstep]);
dX(:,2:end) = r * dt + theta * dS + sigma * sqrt(dS) .* Z;
X = cumsum(dX,2);
S = S0 .* exp(X);
check = mean(S(:,end) / exp(r * T));

	function G = simulateGamma(a,N,M)
		X_loop = rand([N*M 1]).^(1/a);
		Y_loop = rand([N*M 1]).^(1/(1-a));
		while any(X_loop + Y_loop > 1)
			ii = find(X_loop + Y_loop > 1);
			X_loop(ii) = rand([length(ii) 1]).^(1/a);
			Y_loop(ii) = rand([length(ii) 1]).^(1/(1-a));
		end
		X_loop = reshape(X_loop,[N M]);
		Y_loop = reshape(Y_loop,[N M]);
		E = icdf('Exponential', rand([N M]), 1);
		G = X_loop .* E ./ (X_loop + Y_loop);
	end

end

