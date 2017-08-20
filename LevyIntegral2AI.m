function I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,epsilon,optionType,barrierType)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
qnodes1 = linspace(Bl,-epsilon,round(N/2));
qnodes2 = linspace(epsilon,Br,round(N/2));
dq1 = qnodes1(2) - qnodes1(1);
dq2 = qnodes2(2) - qnodes2(1);
I = zeros([N-1 1]);
if nargin == 8 % european option
	for i = 1 : N-1
		I(i) = dq1 * trapz(f_u(nodes(i) + qnodes1,nodes,u,optionType) .* k(qnodes1)) + ...
			dq2 * trapz(f_u(nodes(i) + qnodes2,nodes,u,optionType) .* k(qnodes2));
	end
elseif nargin > 8 % barrier option
	for i = 1 : N-1
		I(i) = dq1 * trapz(f_u(nodes(i) + qnodes1,nodes,u,optionType,barrierType) .* k(qnodes1)) + ...
			dq2 * trapz(f_u(nodes(i) + qnodes2,nodes,u,optionType,barrierType) .* k(qnodes2));
	end
end

end % function

function f = f_u(z,nodes,u,optionType,barrierType)
f = zeros(size(z));
if nargin == 4
	switch optionType
		case 'Call'
			index = find(z > nodes(end));
			f(index) = exp(z(index)) - 1;
			index = find(z > nodes(1) & z < nodes(end));
			f(index) = interp1(nodes,u,z(index));
		case 'Put'
			index = find(z < nodes(1));
			f(index) = 1 - exp(z(index));
			index = find(z > nodes(1) & z < nodes(end));
			f(index) = interp1(nodes,u,z(index));
		otherwise
			error('invalid optionType');
	end
elseif nargin > 4
	switch optionType
		case 'Call'
			switch barrierType
				case 'DO'
					index = find(z > nodes(end));
					f(index) = exp(z(index)) - 1;
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'UO'
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid barrierType')
			end
		case 'Put'
			switch barrierType
				case 'DO'
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'UO'
					index = find(z < nodes(1));
					f(index) = 1 - exp(z(index));
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid barrierType')
			end
		otherwise
			error('invalid optionType');
	end
	
end

end
