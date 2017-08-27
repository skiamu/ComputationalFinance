function [ Payoff ] = optionPayoff(K,optionType,barrierType, barrier)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

switch optionType
	case 'Call'
		switch barrierType
			case 'DO'
				Payoff = @(x) subplus(x(:,end) - K) .* (min(x, [], 2) > barrier);
			case 'UO'
				Payoff = @(x) subplus(x(:,end) - K) .* (max(x, [], 2) < barrier);
			case 'DI'
				Payoff = @(x) subplus(x(:,end) - K) .* (min(x, [], 2) < barrier);
			case 'UI'
				Payoff = @(x) subplus(x(:,end) - K) .* (max(x, [], 2) > barrier);
			otherwise
				error('invalid barrierType)');
		end
	case 'Put'
		switch barrierType
			case 'DO'
				Payoff = @(x) subplus(K - x(:,end)) .* (min(x, [], 2) > barrier);
			case 'UO'
				Payoff = @(x) subplus(K - x(:,end)) .* (max(x, [], 2) < barrier);
			otherwise
				error('invalid barrierType)');	
		end
	otherwise
		error('optioType')
end

end

