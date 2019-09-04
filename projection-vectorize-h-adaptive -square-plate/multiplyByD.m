function [out] = MultiplyByD(D, inp)
	sizeInp = size(inp);
	lengthReshape = (sizeInp(1) * sizeInp(2))/3;
	out = D * reshape(inp, [3, lengthReshape]);
	out = reshape(out, sizeInp);
end