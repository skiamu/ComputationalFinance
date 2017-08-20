function I = LevyIntegral2(k,nodes,u,Br,Bl,N,optionType,barrierType)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
qnodes = linspace(Bl,Br,N);
dq = qnodes(2) - qnodes(1);
I = zeros([N-1 1]);
if nargin == 7 % european option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) + qnodes,nodes,u,optionType) .* k(qnodes));
	end
elseif nargin > 7 % barrier option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) + qnodes,nodes,u,optionType,barrierType) .* k(qnodes));
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
