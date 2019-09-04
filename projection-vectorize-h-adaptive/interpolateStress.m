function [interStress] = interpolateStress(sampleX, sampleStress, interX, orderP)
	% sampleX = [2;7;13];
	% sampleStress = [5;10;20];
	% interX = [0; 4; 10; 15];	

	sampleX = sampleX(:);
	interX = interX(:);
    
	xmax = max([sampleX; interX]);
	xmin = min([sampleX; interX]);

	interX = -1 + 2 * (interX - xmin)/(xmax - xmin);
	sampleX = -1 + 2 * (sampleX - xmin)/(xmax - xmin);

	% assembly P depend on order P
	if orderP == 1
		P = [ones(numel(sampleX), 1), sampleX];
		interP = [ones(numel(interX), 1), interX];
	elseif orderP == 2
		P = [ones(numel(sampleX), 1), sampleX, sampleX.^2];		
		interP = [ones(numel(interX), 1), interX, interX.^2];
	end

	A = P'*P;
	b = P'*sampleStress;

	a = A\b;

	interStress = interP * a;
%     figure;
% 	plot(sampleX, sampleStress(:,1), 'x');
% 	hold on;
% 	plot(interX, interStress(:,1), 'o-')
end