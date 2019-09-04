function [recover_str] = recoverFromModalStress(sigma_star, supel, local, p)

	nnel = supel.nnel;
	nnode = supel.nnode;
	nel = supel.nel;
	sdoff = nnode*2 + 2*(sum(p)-nel);
	eigvalf = supel.values;

	nitanew=[-1:0.05:1];
	sinew=[1:-0.1:0];%later should have bounded and unbounded as input, and if need should have different interval for dif ele to increase no of points for certain ele
    rec_stry = 0;
    rec_strx = 0;
    
	for iel=1:nel		
		% sigma has 3 entries
		sig_star_idx = 1 + (iel-1)*6;

		for i=1:1:length(sinew)
        	for j=1:1:length(nitanew)
            		sii= sinew(i);
            		nitaa= nitanew(j);
					% [N,shapelq,dhdnitalq,Nder]= lobatto(nitaa,p(iel));%row-vectors
% 					[modij, invmodij]=sbfemodijcirlipa(nnel,supel(1).centrecoord(iel).coord,shapelq,dhdnitalq,local(iel).xcrd,local(iel).ycrd);   
% 					[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
					
					% shape function for sigma need to have size 3x6 as sigma has 3 entries - assump p only 1
					% N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(end))];
                    shapeLagr = lagrangian2(nitaa);
					N_for_sigma = [diag(ones(1, 3)*shapeLagr(1)), diag(ones(1, 3)*shapeLagr(2)), diag(ones(1, 3)*shapeLagr(end))];

					% N_for_sigma = [];
					% for iN = 1:numel(shapelq)
					% 	N_for_sigma = [N_for_sigma, diag(ones(1,3))*shapelq(iN)];
					% end
					% si_star = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, i);
					si_star = 0;
					for jj= 1:sdoff
						si_star = real( si_star + supel.c(jj) * sii^(-eigvalf(jj) -1) * N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+8, jj) );
					end
					recover_str(iel).strx(i, j) = si_star(1);
					recover_str(iel).stry(i, j) = si_star(2);
			end    
		end
		
	end
end