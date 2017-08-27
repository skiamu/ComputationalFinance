function f = f_u(z,nodes,u,K,S0,optionType,transf,barrierType)
f = zeros(size(z));
if nargin == 7
	switch transf
		case 'LogMoneyness'
			switch optionType
				case 'Call'
					index = find(z > nodes(end));
					f(index) = exp(z(index)) - 1;
					index = find(z >= nodes(1) & z <= nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'Put'
					index = find(z < nodes(1));
					f(index) = 1 - exp(z(index));
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid optionType');
			end
		case 'LogPrice'
			switch optionType
				case 'Call'
					index = find(z > nodes(end));
					f(index) = S0 * exp(z(index)) - K;
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'Put'
					index = find(z < nodes(1));
					f(index) = K - S0 * exp(z(index));
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid optionType');
			end
		case 'LogPriceFwdPrice'
			switch optionType
				case 'Call'
					index = find(z > nodes(end));
					f(index) = exp(z(index)) - K;
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'Put'
					index = find(z < nodes(1));
					f(index) = K - exp(z(index));
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid optionType');
			end
		case 'Price'
			switch optionType
				case 'Call'
					index = find(z > nodes(end));
					f(index) = z(index) - K;
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				case 'Put'
					index = find(z < nodes(1));
					f(index) = K - z(index);
					index = find(z > nodes(1) & z < nodes(end));
					f(index) = interp1(nodes,u,z(index));
				otherwise
					error('invalid optionType');
			end
		otherwise
	end
	
elseif nargin > 7
	switch transf
		case 'LogMoneyness'
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
		case 'LogPrice'
			switch optionType
				case 'Call'
					switch barrierType
						case 'DO'
							index = find(z > nodes(end));
							f(index) = S0 * exp(z(index)) - K;
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
							f(index) = K - S0 * exp(z(index));
							index = find(z > nodes(1) & z < nodes(end));
							f(index) = interp1(nodes,u,z(index));
						otherwise
							error('invalid barrierType')
					end
				otherwise
					error('invalid optionType');
			end
		case 'LogPriceFwdPrice'
			switch optionType
				case 'Call'
					switch barrierType
						case 'DO'
							index = find(z > nodes(end));
							f(index) = exp(z(index)) - K;
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
							f(index) = K - exp(z(index));
							index = find(z > nodes(1) & z < nodes(end));
							f(index) = interp1(nodes,u,z(index));
						otherwise
							error('invalid barrierType')
					end
				otherwise
					error('invalid optionType');
			end
		case 'Price'
			switch optionType
				case 'Call'
					switch barrierType
						case 'DO'
							index = find(z > nodes(end));
							f(index) = z(index) - K;
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
							f(index) = K - z(index);
							index = find(z > nodes(1) & z < nodes(end));
							f(index) = interp1(nodes,u,z(index));
						otherwise
							error('invalid barrierType')
					end
				otherwise
					error('invalid optionType');
			end
			
		otherwise
			error('inavalid transf')
	end
	
	
end

end
