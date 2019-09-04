function [deltaerror, energyStress, nitaError]= generate_phi(sigma_star, pf, localfine, supelf)

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
    
    % c_save = [];
    % sig_star_save = [];
    % err_star_save = [];
    % c_plus = [];
    % eig_plus = [];
    
	% for iel=1:nel		
	% 	sig_star_el = 0; 
	% 	err_star_el = 0;
		
	% 	% sigma has 3 entries
	% 	sig_star_idx = 1 + (iel-1)*3;

	% 	for i=1:(sdoff-2) %compute sig_star_el    FINE mesh
	% 		for j=1:(sdoff-2)
	% 			sig_star_save(i, j) = 0;
	% 			err_star_save(i, j) = 0;
	% 			for int=1:ngl
	% 				nita = p64(int); 
	% 				wt = w64(int);
	% 				[N,shapelq,dhdnitalq,Nder]= lobatto(nita,pf(iel));%row-vectors
	% 				[modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
	% 				[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
					
	% 				si = D * (-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i) ;
	% 				sj = D * (-eigvalf(j) *B1 +B2) * localfine(iel).phi(:,j) ;
					
	% 				% shape function for sigma need to have size 3x6 as sigma has 3 entries - assump p only 1
	% 				N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(end))];
	% 				si_star = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, i);
	% 				sj_star = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+5, j);
	% 				ei = si_star - si;
	% 				ej = sj_star - sj;
					
	% 				sig_star_el = sig_star_el + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * transpose(si_star) * invD * sj_star * det(modij);
	% 				sig_star_save(i, j) = sig_star_save(i, j) + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * transpose(si_star) * invD * sj_star * det(modij);
	% 				err_star_el = err_star_el + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * transpose(ei) * invD * ej * det(modij);
	% 				err_star_save(i, j) = err_star_save(i, j) + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * transpose(ei) * invD * ej * det(modij);
	% 			end
	% 			c_save(i, j) = (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j)));
 %                c_plus(i, j) = cf(i)*cf(j);
 %                eig_plus(i, j) = eigvalf(i)+eigvalf(j);
	% 		end    
	% 	end
		
	% 	energyStress(iel)= real(-(sig_star_el)); %LHS of eq 6.6
		
	% 	deltaerror(iel) = (real(-err_star_el)); 
		
	% end

	wt_array = ones(3, 1) * w64';	% duplicate so it has 3 entries for each sigma
	wt_array = wt_array(:) * ones(1, sdoff-2); % duplicate the column for sdof phi (sigma)

	c_matrix = supelf.c(1:sdoff-2) * ones(1, sdoff-2);
	eig_matrix = transpose(eigvalf(1:sdoff-2)) * ones(1, sdoff-2);

	c_combine = transpose(c_matrix) .* c_matrix;
	eig_combine = transpose(eig_matrix) + eig_matrix;

	c_over_eig = c_combine ./ (eig_combine);

	for iel=1:nel		
		if max(pf) == 1
			sig_star_idx = 1 + (iel-1)*3;
		else
			sig_star_idx = 1 + (iel-1)*6;		% for quadratic, better change when use different pf for different segments
		end

		si_D = []; si_noD = [];
		si_star_D = []; si_star_noD = [];
		detmodij_col = [];
		for int=1:ngl
			nita = p64(int); 
			[N,shapelq,dhdnitalq,Nder]= lobatto(nita,pf(iel));%row-vectors
			[modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
			[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);


			modalStressTemp_D = [];
			modalStressTemp_noD = [];
			for idoff = 1:sdoff-2
				stressTemp = (-eigvalf(idoff) *B1 +B2) * localfine(iel).phi(:,idoff);
				modalStressTemp_D = [modalStressTemp_D,  D * stressTemp] ;
				modalStressTemp_noD = [modalStressTemp_noD,  stressTemp] ;		
			end
			si_D = [si_D; modalStressTemp_D];
			si_noD = [si_noD; modalStressTemp_noD];

			modalStressTemp_D = [];
			modalStressTemp_noD = [];

			shapeLagr = lagrangian2(nita);
			if max(pf) == 1
				N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(end))];
				shift_constant = 5;
			else
				N_for_sigma = [diag(ones(1, 3)*shapeLagr(1)), diag(ones(1, 3)*shapeLagr(2)), diag(ones(1, 3)*shapeLagr(end))];
				shift_constant = 8;
			end
			
			for idoff = 1:sdoff-2
				si_star_temp = N_for_sigma * sigma_star(sig_star_idx:sig_star_idx+shift_constant, idoff);
				modalStressTemp_D = [modalStressTemp_D,  si_star_temp] ;
				modalStressTemp_noD = [modalStressTemp_noD,  invD*si_star_temp] ;
			end
			si_star_D = [si_star_D; modalStressTemp_D];
			si_star_noD = [si_star_noD; modalStressTemp_noD];
			detmodij_col = [detmodij_col; det(modij)];
		end
		error_si_D = si_star_D - si_D;
		error_si_noD = si_star_noD - si_noD;

		detmodij_mat = ones(3, 1) * detmodij_col';	% duplicate so it has 3 entries for each sigma
		detmodij_mat = detmodij_mat(:) * ones(1, sdoff-2); % duplicate the column for sdof phi (sigma)

        temp = transpose(wt_array .* detmodij_mat .* si_star_D) * si_star_noD;
		energyStress(iel) = real(-(sum(sum(c_over_eig .* temp))));
        temp = transpose(wt_array .* detmodij_mat .* error_si_D) * error_si_noD;
		deltaerror(iel) = real(-(sum(sum(c_over_eig .* temp))));
    end
    totalEnergy = sqrt( mean(energyStress) );
	nitaError = sqrt(deltaerror) ./ totalEnergy;
end