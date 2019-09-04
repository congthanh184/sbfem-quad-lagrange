function [output] = convertSingleArrayToEyeMatrix(a, p)
	output = [];
	for i = 1:numel(a)
		output = [output; eye(p) *  a(i)];
	end
end