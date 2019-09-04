function [recover_str]= recover_modal_stress(pf, localfine, supelf, uf)

	%prog for ref sol using the pseudo code in book
	%start with 
	nel     = supelf.nel;
	nnel    = supelf.nnel;
	nnode   = supelf.nnode;
	emodule = supelf.emodule;
	poisson = supelf.poisson; 

	% [localfine,supelf,uf]= today3(p+ones(nel,1));
	cf = supelf(1).c;
	eigvalf = supelf(1).values;
	sdoff = nnode*2 + 2*(sum(pf)-nel);

	%uniform refined mesh all upgrade to p+1, using this as the correct solution to derived energyStress for coarse mesh and define which to refine
	ngl=16; %can go up to 10 N'*N*det(J) ~ 3*p
	[p64,w64] = sbfeglqd1(ngl);
	[matmtrx,G] = sbfematisolq(1,emodule,poisson);
	D = G*matmtrx;
	invD = inv(D);

	energyStress = zeros(nel,1);% contribution of each ele in energyStress eq 6.6
	deltaerror = zeros(nel,1); % energyStress rate for each ele eq 6.5

	% calculate sigma_star at each node for each model
	% sigma_star has number_nodes rows and n mode column
	sigma_star = [];
	ssNode_1_save = [];
	for iel = 1:nel
		ssNode_1 = [];
		ssNode_m1 = [];
		for i=1:(sdoff) %compute e1t    FINE mesh
			[N,shapelq,dhdnitalq,Nder]= lobatto(-1,pf(iel));%row-vectors
			[modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
			[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);					
			ssNode_m1 = [ssNode_m1,  D * (-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i)] ;

			[N,shapelq,dhdnitalq,Nder]= lobatto(1,pf(iel));%row-vectors

			[modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
			[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);					
			ssNode_1 = [ssNode_1,  D * (-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i)] ;
		end
		if ~isempty(ssNode_1_save)
			ssNode_m1 = (ssNode_m1 + ssNode_1_save)./2;
		end
		sigma_star = [sigma_star; ssNode_m1];
		ssNode_1_save = ssNode_1;
	end
	sigma_star = [sigma_star; ssNode_1_save];

	nitanew=[-1:0.05:1];
	sinew=[1:-0.1:0];%later should have bounded and unbounded as input, and if need should have different interval for dif ele to increase no of points for certain ele
    rec_stry = 0;
    rec_strx = 0;
    
	for iel=1:nel		
		% sigma has 3 entries
		sig_star_idx = 1 + (iel-1)*3;

		for i=1:1:length(sinew)
        	for j=1:1:length(nitanew)
            		sii= sinew(i);
            		nitaa= nitanew(j);
					[N,shapelq,dhdnitalq,Nder]= lobatto(nitaa,pf(iel));%row-vectors
					[modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
					[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
					
					% shape function for sigma need to have size 3x6 as sigma has 3 entries - assump p only 1
					N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(2))];
					% si_star = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, i);
					si_star = 0;
					for jj= 1:sdoff
						si_star = real( si_star + supelf.c(jj) * sii^(-eigvalf(jj) -1) * N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, jj) );
					end
					recover_str(iel).rec_strx(i, j) = si_star(1);
					recover_str(iel).rec_stry(i, j) = si_star(2);
			end    
		end
		
	end

end