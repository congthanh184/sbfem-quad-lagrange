function [recover_str] = recoverFromModalStress_p2_vectorize(sigma_star, supel, local, p)

	nnel = supel.nnel;
	nnode = supel.nnode;
	nel = supel.nel;
	sdoff = nnode*2 + 2*(sum(p)-nel);
	supelvalues = supel.values;
	supelc = supel.c;

% 	nitanew=[-1:0.05:1];
% 	sinew=[1:-0.1:0];%later should have bounded and unbounded as input, and if need should have different interval for dif ele to increase no of points for certain ele
    nitanew = local(1).nitanew;
    sinew = local(1).sinew;
    rec_stry = 0;
    rec_strx = 0;
    
    global_p = 2;
	for iel=1:nel		
		% assump all element has p = 2 -> sigma has 3 entries, dof of each node = 2 => (p+1)*(nnel)
		% each sigma is a column 3x1 [sigmaX; sigmaY; sigmaT], therefore element matrix size = [3x3, number of phi]
		sig_star_idx = 1 + (iel-1)*(global_p+1)*(nnel);
		element_sigma_star = sigma_star(sig_star_idx:sig_star_idx + (global_p+1)^2 - 1, :);

		% set-up shape function matrix 
	    N_mat = [];

	    for i = 1:1:numel(nitanew)
           	shapeLagr = lagrangian2(nitanew(i));
			N_for_sigma = [diag(ones(1, 3)*shapeLagr(1)), diag(ones(1, 3)*shapeLagr(2)), diag(ones(1, 3)*shapeLagr(end))];
			N_mat = [N_mat; N_for_sigma];
	    end

	    % init matrix 
		c_mat = repmat(supelc(:), [1 numel(sinew)]);
    	si_mat = repmat(sinew(:)', [numel(supelc) 1]);
    	lambda_mat = repmat(supelvalues(:), [1 numel(sinew)]);

    	c_sii_supelvalues_minus1_mat = c_mat .* (si_mat .^ (-lambda_mat - 1));

    	% calculate sigma_star_nita[i, j] = (N * sigma_star), with i is nita, j is number of phi
    	% stressXYT[i, j], with i is nita, j is sii
    	stressXYT = real((N_mat * element_sigma_star) * c_sii_supelvalues_minus1_mat);
    	% transpose XYT to have (i, j) with i is sii, and j is nita
    	stressXYT = stressXYT';
    	recover_str(iel).strx = stressXYT(:, 1:3:end);
    	recover_str(iel).stry = stressXYT(:, 2:3:end);

% 		for i=1:1:length(sinew)
%         	for j=1:1:length(nitanew)
%             		sii= sinew(i);
%             		nitaa= nitanew(j);
% 					% [N,shapelq,dhdnitalq,Nder]= lobatto(nitaa,p(iel));%row-vectors
% % 					[modij, invmodij]=sbfemodijcirlipa(nnel,supel(1).centrecoord(iel).coord,shapelq,dhdnitalq,local(iel).xcrd,local(iel).ycrd);   
% % 					[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
					
% 					% shape function for sigma need to have size 3x9 as sigma has 3 entries - p = 2 => 6 for 2 end node + 3 for middle node
% 					% N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(end))];
%                     shapeLagr = lagrangian2(nitaa);
% 					N_for_sigma = [diag(ones(1, 3)*shapeLagr(1)), diag(ones(1, 3)*shapeLagr(2)), diag(ones(1, 3)*shapeLagr(end))];

% 					% N_for_sigma = [];
% 					% for iN = 1:numel(shapelq)
% 					% 	N_for_sigma = [N_for_sigma, diag(ones(1,3))*shapelq(iN)];
% 					% end
% 					% si_star = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, i);
% 					si_star = 0;
% 					for jj= 1:sdoff
% 						si_star = real( si_star + supel.c(jj) * sii^(-eigvalf(jj) -1) * N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+8, jj) );
% 					end
% 					recover_str(iel).strx(i, j) = si_star(1);
% 					recover_str(iel).stry(i, j) = si_star(2);
% 			end    
% 		end
		
	end
end