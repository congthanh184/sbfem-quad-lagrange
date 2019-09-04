function [childLarge, childSmall] = separateMatrix(parentMat, threshold)
% separate parentMat to 2 matrix depend on threshold
% parentMat > threshold -> childLarge
% parentMat < threshold -> childSmall

idxMatLarger1 = abs(real(parentMat)) > threshold;
childLarge = zeros(size(parentMat));
childSmall = zeros(size(parentMat));
childLarge(idxMatLarger1) = parentMat(idxMatLarger1);
childSmall(~idxMatLarger1) = parentMat(~idxMatLarger1);
end