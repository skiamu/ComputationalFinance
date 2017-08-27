function [ Psi ] = CharacteristicExp( model, param )
%Characteristic_Exp restituisce l'esponente caratteristico del modello
%indicato in input
%
%   INPUT:
%         model = striga che indica il modello di cui si vuole l'espoinente
%                 caratteristico
%
%         param =  vettore parametri del modello
%         --> 'Heston'
%                     param.k;
%                     param.theta;
%                     param.epsilon;
%                     param.rho;
%                     param.V0;
%
%         --> 'BS'
%                     param(1) = sigma
%         --> 'Kou'
%                     param.sigma;
%                     param.lambda;
%                     param.lambda_p;
%                     param.lambda_m;
%                     param.p;
%
%         --> 'Merton'
%                     param(1) = lambda;
%                     param(2) = sigma;
%                     param(3) = mu;
%                     param(4) = delta;
%
%         --> 'VG'
%                     param(1) = sigma
%                     param(2) = kappa
%   OUTPUT:
%         Psi = esponente caratteristico (function handle)


switch model
	
	case 'Heston'
		
		% parametri
		k = param.k;
		theta = param.theta;
		rho = param.rho;
		epsilon = param.epsilon;
		V0 = param.V0;
		
		gamma = @(u) k - rho * epsilon .* u * 1i;
		
		alpha = @(u) - 0.5 .* (u.^2 + 1i .* u);
		
		d = @(u) sqrt( gamma(u).^2 - 2 * (epsilon^2) *  alpha(u));
		
		C = @(t,u) (-2 * k * theta / (epsilon^2)) .* ...
			(2 * log( (2 * d(u) - (d(u) - gamma(u)) .* ...
			(1 - exp(-d(u) * t)))./ (2*d(u)) ) + (d(u) - gamma(u)) * t  );
		
		B = @(t,u) ((2 * alpha(u) .* (1 - exp(- d(u) * t))) * V0) ./ ...
			((2 * d(u) - (d(u) - gamma(u)) .* (1 - exp(-d(u) * t))));
		
		Psi = @(t,u) B(t,u) + C(t,u);
		
	case 'Merton'
		% estrazione parametri
		% parametri Merton: 1) sigma = diffusion volatility
		%                   2) lambda = jump intensity
		%                   3) mu = mean jump
		%                   4) delta = jump std
		sigma = param(1);
		lambda = param(2);
		mu = param(3);
		delta = param(4);
		
		% è l'esponente caratteristico con drift nullo
		Psi = @(u) - 0.5 * u.^2 * sigma^2 + lambda * (exp(- 0.5 * delta^2 *...
			u.^2 + 1i * mu * u) - 1);
		
	case 'Kou'
		% estrazione parametri
		% parametri Kou: 1)   sigma = diffusion volatility;
		%                     lambda = jump intensity
		%                     lambda_p = positive jump parameter
		%                     param_m = negative jump parameter
		%                     p = probability positive jump
		sigma = param(1);
		lambda = param(2);
		lambda_p = param(3);
		lambda_m = param(4);
		p = param(5);
		
		% esponente caratteristico con drift nullo
		Psi = @(u) (-sigma^2 .* u.^2) / 2 + 1i .* u .* lambda .* (p ./ (lambda_p - ...
			1i .* u) - (1 - p) ./ (lambda_m + 1i .* u));
		
	case 'NIG'
		sigma = param(1);
		kappa = param(2);
		
		% esponente caratteristico con drift nullo
		Psi = @(u) 1 / kappa - (1 / kappa) * sqrt( 1 + u.^2 * sigma^2 * kappa);
		
	case 'VG'
		% parametri VG: 1) theta = Brownian drift to be subordinated
		%               2) sigma = Brownian volatility to be subordinated
		%               3) kappa = subordinator volatility
		theta = param(1);
		sigma = param(2);
		kappa = param(3);
		if length(param) == 4
			sigmaDiff = param(4);
		else
			sigmaDiff = 0;
		end
		% esponente caratteristico con drift nullo
		Psi = @(u) - (1 / kappa) * log(1 + 0.5 * u.^2 * sigma^2 * kappa - ...
			1i * theta * kappa * u) - 0.5 * sigmaDiff^2 * u.^2;
		
	case 'BS'
		sigma = param(1);
		
		% esponente caratteristico con drift nullo
		Psi = @(v)  - (sigma^2) / 2 .* v.^2;
	case 'alpha'
		alpha = param(1);
		sigma = param(2);
		
		% esponente caratteristico con drift nullo
		Psi = @(u) -sigma^alpha * abs(u).^alpha;
	otherwise
		error('model invalid : %s', model);
		
end


end

