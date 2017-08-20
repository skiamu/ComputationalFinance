function [ S, check ] = assetMerton(S0,r,T,param,Nsim,Nstep)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
method = 1;
sigma = param(1); % diffusion volatility
lambda = param(2); % jump intensity
mu = param(3); % mean jump size
delta = param(4); % std jump size
Psi_X = @(v) -sigma^2/2*v.^2 + lambda*(exp(-delta^2/2*v.^2 + 1i*mu*v) - 1);
driftRN = -Psi_X(-1i); % drift risk neutral
dt = T / Nstep;
t = 0:dt:T;
dX = zeros([Nsim Nstep+1]);
Z = randn([Nsim Nstep]);
dX(:,2:end) = r * dt + driftRN * dt + sigma * sqrt(dt) * Z;
X = cumsum(dX,2);
CP = simulateCP(method);
S = S0 .* exp(X + CP);
check = mean(S(:,end) / exp(r * T));

	function CP = simulateCP(method)
		% method = 1 if incremental 2 if conditional
		CP = zeros([Nsim Nstep+1]); % same dimension of X
		if method == 1 % incremental method
			for i = 1 : Nsim
				T_i = 0;
				while T_i < T
					T_i = T_i + icdf('Exponential',rand, 1 / lambda);
					Y_i = mu + delta * randn;
					idxAfterTi = t >= T_i;
					CP(i,idxAfterTi) = CP(i,idxAfterTi) + Y_i;
				end
			end
		else % conditional method
			N = icdf('Poisson', rand([Nsim 1]), lambda * T); % number of jumps
			for i = 1 : Nsim
				T_i = rand([1 N(i)]) * T; % jump instant
				Y = mu + delta * randn([1 N(i)]); % jump size
				for j = 1 : N(i)
					idxAfterTi = t >= T_i(j);
					CP(i,idxAfterTi) = CP(i,idxAfterTi) + Y(j);
				end
			end
		end	
	end
end % function

