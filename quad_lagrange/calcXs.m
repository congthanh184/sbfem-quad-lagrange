%% calcXs: calculating Xs and Ys from X, Y and center coordinators
function [output] = calcXs(x, y, x0, y0, xi)
	xs = (x - x0) / xi;
	ys = (y - y0) / xi;
    output = [xs; ys];
end