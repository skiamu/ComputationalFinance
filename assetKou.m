function [ S, check ] = assetKou(S0,r,T,param,Nsim,Nstep)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
sigma = param(1);
lambda = param(2);
lambda_p = param(3);
lambda_m = param(4);
p = param(5);
Psi_X = @(v) -sigma^2/2*v.^2 + 1i*v*lambda.*(p./(lambda_p - 1i*v) -...
	(1-p)./(lambda_m + 1i*v)); % characteristic exponent
driftRN = -Psi_X(-1i);
method = 2;

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
		CP = zeros([Nsim Nstep+1]);
		if method == 1 % incremental
			for i = 1 : Nsim
				T_i = 0;
				while T_i < T
					T_i = T_i + icdf('Exponential', rand, 1 / lambda);
					if rand < p
						Y_i = icdf('Exponential', rand, 1 / lambda_p);
					else
						Y_i = -icdf('Exponential', rand, 1 / lambda_m);
					end
					idxAfterTi = t >= T_i;
					CP(i,idxAfterTi) = CP(i,idxAfterTi) + Y_i;
				end
			end
		else % conditional
			N = icdf('Poisson', rand([Nsim 1]), lambda * T);  % number of jumps
			for i = 1 : Nsim
				T_i = rand([1 N(i)]) * T; % jump instant
				Y = zeros([1 N(i)]);
				for ii = 1 : N(i)
					if rand < p
						Y(ii) = icdf('Exponential', rand, 1 / lambda_p);
					else
						Y(ii) = -icdf('Exponential', rand, 1 / lambda_m);
					end
					idxAfterTi = t >= T_i(ii);
					CP(i,idxAfterTi) = CP(i,idxAfterTi) + Y(ii);
				end
			end	
		end
	end % end simulateCP

end % end function

