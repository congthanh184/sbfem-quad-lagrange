%prog for using energy norm in stress as in paper

%start with 
%p=[1;1;1;1];

function [localcoarse,array,deltaerror,rhs65,error_array,rhs66,errest, ET, c_c_over_eig]= h_project_error_estimator_vectorize(p, supel1, pf, supel2)

%prog for ref sol using the pseudo code in book
%start with 
% nel=4;
% nnel=2;
% nnode=5;
% emodule=250000;
% poisson=0.3; 

nel     = supel1.nel;
nnel    = supel1.nnel;
nnode   = supel1.nnode;
emodule = supel1.emodule;
poisson = supel1.poisson; 


% [localcoarse,supelc,uc]= today3(p);%current mesh
[localcoarse,supelc,uc] = today3_modify(p, supel1);
% [localfine,supelf,uf]= today3(p+ones(nel,1));
[localfine,supelf,uf] = today3_modify(pf, supel2);

pc = p; % pf = p+ones(nel,1);
cc = supelc(1).c; cf = supelf(1).c;
eigvalc = supelc(1).values; eigvalf = supelf(1).values;
sdofc = nnode*2 + 2*(sum(pc)-nel);
sdoff = nnode*2 + 2*(sum(pf)-nel);

% project localf on localc, c and lambda of finecoarse actually the same with the fine
[localfc] = project_local(pf, localfine, supelf, pc, localcoarse, supelc);


% Need to look
% ==================================================================
for iel=1:nel
    pfsum(iel)= sum(pf(1:iel));
end

ndf2c = supelf(1).nodes; % =ndf need to extract to obtain ndf2c
for i=length(pfsum):-1:1
    ndf2c(pfsum(i))= [];
end

indexndf2c = sbfeeldof(ndf2c,2);
%temp1 = supelf(1).phi1(:,indexndf2c);
%temp2 = temp1(indexndf2c,:);
%cf2c = inv(temp2)*supelf(1).globnodesdisp(indexndf2c);
% cf2c = inv(supelc(1).phi1)*uf(indexndf2c);
% cc = cf2c;
% ==================================================================

%uniform refined mesh all upgrade to p+1, using this as the correct solution to derived error for coarse mesh and define which to refine
ngl=16; %can go up to 10 N'*N*det(J) ~ 3*p
[p64,w64] = sbfeglqd1(ngl);
[matmtrx,G] = sbfematisolq(1,emodule,poisson);
D = G*matmtrx;
tol= 1;%tolerate error =1%

invD = inv(D);

error_array = zeros(nel,1);% contribution of each ele in error eq 6.6
deltaerror = zeros(nel,1); % error rate for each ele eq 6.5
errestop = 0;
errestbot = 0;
errest = 0; %this the relative error norm to count for stop criteria 

ET = [];
% for iel=1:nel

%     e1t=0; e2t=0; e3t=0; e4t=0; %error1,2,3,4 to correspond to 4 sub integration
%     E =zeros(sdoff-2,sdoff-2);
%     delta=0;
    
%     for i=1:(sdoff-2) %compute e1t    FINE mesh
%         for j=1:(sdoff-2)
%             for int=1:ngl
%                 nita = p64(int); 
%                 wt = w64(int);
%                 [N,shapelq,dhdnitalq,Nder]= lobatto(nita,pf(iel));%row-vectors
%                 [modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
%                 [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
                
%                 e1i = cf(i) * D * (-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i) ;
%                 e1j = cf(j) * D * (-eigvalf(j) *B1 +B2) * localfine(iel).phi(:,j) ;
                
%                 %e1t = e1t + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * ((-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i))' * (G*matmtrx)' * ((-eigvalf(j) *B1 +B2) * localfine(iel).phi(:,j)) * det(modij);
%                 e1t = e1t + wt * (1/(eigvalf(i)+eigvalf(j))) * transpose(e1i) * inv(D) * e1j * det(modij);
                
%                 e5i = cf(i) * ((-eigvalf(i) * B1 + B2) * localfine(iel).coeffstr(:,i) - (-eigvalf(i) * B1 + B2) * localfine(iel).coefcstr(:,i));
%                 e5j = cf(j) * ((-eigvalf(j) * B1 + B2) * localfine(iel).coeffstr(:,j) - (-eigvalf(j) * B1 + B2) * localfine(iel).coefcstr(:,j));
%                 delta = delta + wt * (1/(eigvalf(i)+eigvalf(j))) * transpose(e5i) * (G*matmtrx)' * e5j * det(modij);
%             end
%         end    
%     end
    
%     for i=1:(sdofc-2) %compute e4t    COARSE mesh
%         for j=1:(sdofc-2)
%             for int=1:ngl
%                 nita = p64(int); 
%                 wt = w64(int);
%                 [N,shapelq,dhdnitalq,Nder]= lobatto(nita,pc(iel));%row-vectors
%                 [modij, invmodij]=sbfemodijcirlipa(nnel,supelc(1).centrecoord(iel).coord,shapelq,dhdnitalq,localcoarse(iel).xcrd,localcoarse(iel).ycrd);   
%                 [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
                
%                 e4i = cc(i) * D * (-eigvalc(i) *B1 +B2) * localcoarse(iel).phi(:,i) ;
%                 e4j = cc(j) * D * (-eigvalc(j) *B1 +B2) * localcoarse(iel).phi(:,j) ;
%                 e4t = e4t + wt * (1/(eigvalc(i)+eigvalc(j))) * transpose(e4i) * inv(D) * e4j * det(modij);
                
%             end
%         end    
%     end
    
%     for i=1:(sdoff-2) %compute e2t & e3t    FINE&COARSE meshes
%         for j=1:(sdofc-2)
%             for int=1:ngl
%                 nita = p64(int); 
%                 wt = w64(int);
%                 [Nf,shapelqf,dhdnitalqf,Nderf]= lobatto(nita,pf(iel));%row-vectors
%                 [modijf, invmodijf]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelqf,dhdnitalqf,localfine(iel).xcrd,localfine(iel).ycrd);   
%                 [b1f, b2f, B1f, B2f]=sbfederiv2(nnel,shapelqf,dhdnitalqf,invmodijf);
                
%                 [Nc,shapelqc,dhdnitalqc,Nderc]= lobatto(nita,pc(iel));%row-vectors
%                 [modijc, invmodijc]=sbfemodijcirlipa(nnel,supelc(1).centrecoord(iel).coord,shapelqc,dhdnitalqc,localcoarse(iel).xcrd,localcoarse(iel).ycrd);   
%                 [b1c, b2c, B1c, B2c]=sbfederiv2(nnel,shapelqc,dhdnitalqc,invmodijc);
                
%                 e2t = e2t + wt * (cf(i)*cc(j)/(eigvalf(i)+eigvalc(j))) * transpose((-eigvalf(i) *B1f +B2f) * localfine(iel).phi(:,i)) * (G*matmtrx)' * ((-eigvalc(j) *B1c +B2c) * localcoarse(iel).phi(:,j)) * det(modijf);
%                 e3t = e3t + wt * (cf(i)*cc(j)/(eigvalf(i)+eigvalc(j))) * transpose((-eigvalc(j) *B1c +B2c) * localcoarse(iel).phi(:,j)) * (G*matmtrx)' * ((-eigvalf(i) *B1f +B2f) * localfine(iel).phi(:,i)) * det(modijf);
%             end
%         end    
%     end
    
%     error_array(iel)= sqrt(real(-(e1t-e2t-e3t+e4t))); %LHS of eq 6.6
%     errestop = errestop + real(-(e1t-e2t-e3t+e4t));
%     errestbot = errestbot + real(-e1t);
    
%     deltaerror(iel) = sqrt(real(-delta)); 
    
% end

wt_array = ones(3, 1) * w64';   % duplicate so it has 3 entries for each sigma
fwt_array = wt_array(:) * ones(1, sdoff); % duplicate the column for sdof phi (sigma)
cwt_array = wt_array(:) * ones(1, sdofc); % duplicate the column for sdof phi (sigma)
%====================
f_c_matrix = supelf.c(1:sdoff) * ones(1, sdoff);
c_c_matrix = supelc.c(1:sdofc) * ones(1, sdofc);

f_eig_matrix = transpose(eigvalf(1:sdoff)) * ones(1, sdoff);
c_eig_matrix = transpose(eigvalc(1:sdofc)) * ones(1, sdofc);

f_c_combine = transpose(f_c_matrix) .* f_c_matrix;
f_eig_combine = transpose(f_eig_matrix) + f_eig_matrix;
f_c_over_eig = f_c_combine ./ (f_eig_combine);

[zeroRow, zeroCol] = find(isnan(f_c_over_eig));
f_c_over_eig(zeroRow, zeroCol) = 0;
%====================
c_c_combine = transpose(c_c_matrix) .* c_c_matrix;
c_eig_combine = transpose(c_eig_matrix) + c_eig_matrix;
c_c_over_eig = c_c_combine ./ (c_eig_combine);

[zeroRow, zeroCol] = find(isnan(c_c_over_eig));
c_c_over_eig(zeroRow, zeroCol) = 0;
%====================
fc_c_matrix = supelf.c(1:sdoff) * ones(1, sdofc);
cf_c_matrix = supelc.c(1:sdofc) * ones(1, sdoff);
fc_eig_matrix = transpose(eigvalf(1:sdoff)) * ones(1, sdofc);
cf_eig_matrix = transpose(eigvalc(1:sdofc)) * ones(1, sdoff);

fc_c_combine = transpose(cf_c_matrix) .* fc_c_matrix;
fc_eig_combine = transpose(cf_eig_matrix) + fc_eig_matrix;
fc_c_over_eig = fc_c_combine ./ (fc_eig_combine);

[zeroRow, zeroCol] = find(isnan(fc_c_over_eig));
fc_c_over_eig(zeroRow, zeroCol) = 0;

cf_c_over_eig = transpose(fc_c_over_eig);

for iel = 1:nel
    % fsigma_D = []; fsigma_noD = [];
    % sigma_project_D = []; sigma_project_noD = [];
    % csigma_D = []; csigma_noD = [];
    % fdetmodij_col = []; cdetmodij_col = [];
    one_fsigma_noD = []; two_fsigma_noD = [];
    one_sigma_project_noD = []; two_sigma_project_noD = [];
    half1_csigma_noD = []; half2_csigma_noD = [];
    full_csigma_noD = [];

    f1detmodij_col = []; f2detmodij_col = [];
    half1_cdetmodij_col = []; half2_cdetmodij_col = [];
    full_cdetmodij_col = [];

    for idx = 1:ngl
        nita = p64(idx);

        % nita is the index of local shape funtion [-1 1]
        % want to map nita to nita_global of coarse mesh, support by 2 fine element
        % nita_global = 1/2*(nita - 1) (nita = [-1 1] => nita_global = [-1 0])
        % nita_global = 1/2*(nita + 1) (nita = (-1 1] => nita_global = [0 1])
        nita_glob1 = 0.5 * (nita - 1);
        nita_glob2 = 0.5 * (nita + 1);
        iel_fine1 = 2 * iel - 1;
        iel_fine2 = 2 * iel;

        % calculate shape information for 1st half, 2nd half, full coarse element
        [F1_N, F1_modij, F1_B1, F1_B2] = lobatto_calcShapeFunctions(nita, pf, iel_fine1, supelf, localfine);
        [F2_N, F2_modij, F2_B1, F2_B2] = lobatto_calcShapeFunctions(nita, pf, iel_fine2, supelf, localfine);
        [C1_N, C1_modij, C1_B1, C1_B2] = lobatto_calcShapeFunctions(nita_glob1, pc, iel, supelc, localcoarse);
        [C2_N, C2_modij, C2_B1, C2_B2] = lobatto_calcShapeFunctions(nita_glob2, pc, iel, supelc, localcoarse);
        [C_N, C_modij, C_B1, C_B2] = lobatto_calcShapeFunctions(nita, pc, iel, supelc, localcoarse);

        one_modalStressTemp_noD = [];
        two_modalStressTemp_noD = [];

        one_projectionTemp_noD = [];
        two_projectionTemp_noD = [];

        for idoff = 1:sdoff
            % check nita to find the correct local element, here nita lead iel_local, and B1f B2f already map to local
            % element 1 of fine mesh 
            stressTemp = (-eigvalf(idoff) * F1_B1 + F1_B2) * localfine(iel_fine1).phi(:,idoff);
            one_modalStressTemp_noD = [one_modalStressTemp_noD,  stressTemp] ;      
            % element 2 of fine mesh
            stressTemp = (-eigvalf(idoff) * F2_B1 + F2_B2) * localfine(iel_fine2).phi(:,idoff);
            two_modalStressTemp_noD = [two_modalStressTemp_noD,  stressTemp] ;      

            % this line below in p-projection method, B1f & B2f should be B1c and B2c (according to the papers)
            % however, as we are using lobatto and the localfine.coefcstr has the same structure with the localfine.phi
            % but has zero rows at the (p+1) nodes, therefore we can use B1f and B2f
            % stressProjectionTemp = (-eigvalf(idoff) * B1f + B2f) * localfine(iel).coefcstr(:,idoff);

            % use the same c and lambda of finemesh, but B1 & B2 is from coarsemesh
            % 1st half of main element (-1, 0)
            stressProjectionTemp = (-eigvalf(idoff) * C1_B1 + C1_B2) * localfc(iel).phi(:, idoff);            
            one_projectionTemp_noD  = [one_projectionTemp_noD, stressProjectionTemp];
            % 2nd half of main element (0, +1)
            stressProjectionTemp = (-eigvalf(idoff) * C2_B1 + C2_B2) * localfc(iel).phi(:, idoff);            
            two_projectionTemp_noD  = [two_projectionTemp_noD, stressProjectionTemp];
        end

        one_fsigma_noD = [one_fsigma_noD; one_modalStressTemp_noD];
        two_fsigma_noD = [two_fsigma_noD; two_modalStressTemp_noD];

        one_sigma_project_noD   = [one_sigma_project_noD;   one_projectionTemp_noD];
        two_sigma_project_noD   = [two_sigma_project_noD;   two_projectionTemp_noD];

        half1_modalStressTemp_noD = [];
        half2_modalStressTemp_noD = [];
        full_modalStressTemp_noD = [];

        for idofc = 1:sdofc
            % half of the element
            stressTemp = (-eigvalc(idofc) * C1_B1 + C1_B2) * localcoarse(iel).phi(:,idofc);
            half1_modalStressTemp_noD = [half1_modalStressTemp_noD,  stressTemp] ;      

            % 2nd half
            stressTemp = (-eigvalc(idofc) * C2_B1 + C2_B2) * localcoarse(iel).phi(:,idofc);
            half2_modalStressTemp_noD = [half2_modalStressTemp_noD,  stressTemp] ;      

            % full 
            stressTemp = (-eigvalc(idofc) * C_B1 + C_B2) * localcoarse(iel).phi(:,idofc);
            full_modalStressTemp_noD = [full_modalStressTemp_noD,  stressTemp] ;      
        end

        half1_csigma_noD    = [half1_csigma_noD;    half1_modalStressTemp_noD];
        half2_csigma_noD    = [half2_csigma_noD;    half2_modalStressTemp_noD];
        full_csigma_noD     = [full_csigma_noD;     full_modalStressTemp_noD];

        f1detmodij_col = [f1detmodij_col; det(F1_modij)];
        f2detmodij_col = [f2detmodij_col; det(F2_modij)];
        half1_cdetmodij_col = [half1_cdetmodij_col; det(C1_modij)];
        half2_cdetmodij_col = [half2_cdetmodij_col; det(C2_modij)];
        full_cdetmodij_col = [full_cdetmodij_col; det(C_modij)];
    end


    %==================================================
%     if (numel(unique(f1detmodij_col)) > 1) 
%         disp 2
%     end
    if iel==14
        disp 2
    end

    f1detmodij_mat = ones(3, 1) * f1detmodij_col';  % duplicate so it has 3 entries for each sigma
    f1detmodij_mat = f1detmodij_mat(:) * ones(1, sdoff); % duplicate the column for sdof phi (sigma)
    f2detmodij_mat = ones(3, 1) * f2detmodij_col';  % duplicate so it has 3 entries for each sigma
    f2detmodij_mat = f2detmodij_mat(:) * ones(1, sdoff); % duplicate the column for sdof phi (sigma)

    % temp = transpose(fwt_array .* fdetmodij_mat .* fsigma_D) * fsigma_noD;

    % f12det_K = 10 * max([f1detmodij_col; f2detmodij_col]);
    % f1detTemp = f1detmodij_mat / f12det_K;
    % f2detTemp = f2detmodij_mat / f12det_K;

    temp1 = transpose(fwt_array .* f1detmodij_mat .* multiplyByD( D, one_fsigma_noD)) * one_fsigma_noD;
    temp2 = transpose(fwt_array .* f2detmodij_mat .* multiplyByD( D, two_fsigma_noD)) * two_fsigma_noD;
    % temp1 = transpose(fwt_array .* f1detTemp .* multiplyByD( matmtrx, one_fsigma_noD)) * one_fsigma_noD;
    % temp2 = transpose(fwt_array .* f2detTemp .* multiplyByD( matmtrx, two_fsigma_noD)) * two_fsigma_noD;
    e1t = sum(sum(f_c_over_eig .* temp1)) + sum(sum(f_c_over_eig .* temp2));
    % e1t = G * f12det_K * sum(sum(f_c_over_eig .* (temp1 + temp2)));

    % temp1 = transpose(fwt_array .* multiplyByD( matmtrx, one_fsigma_noD)) * one_fsigma_noD;
    % temp2 = transpose(fwt_array .* multiplyByD( matmtrx, two_fsigma_noD)) * two_fsigma_noD;
    % e1t = G * f1detmodij_col(1) * sum(sum(f_c_over_eig .* temp1)) + G * f2detmodij_col(1) * sum(sum(f_c_over_eig .* temp2));

    one_error_si_noD = one_fsigma_noD - one_sigma_project_noD;
    two_error_si_noD = two_fsigma_noD - two_sigma_project_noD;
    temp1 = transpose(fwt_array .* f1detmodij_mat .* multiplyByD( D, one_error_si_noD)) * one_error_si_noD;
    temp2 = transpose(fwt_array .* f2detmodij_mat .* multiplyByD( D, two_error_si_noD)) * two_error_si_noD;
    % temp1 = transpose(fwt_array .* f1detTemp .* multiplyByD( matmtrx, one_error_si_noD)) * one_error_si_noD;
    % temp2 = transpose(fwt_array .* f2detTemp .* multiplyByD( matmtrx, two_error_si_noD)) * two_error_si_noD;
    delta = sum(sum(f_c_over_eig .* temp1)) + sum(sum(f_c_over_eig .* temp2));
    
    % idxMatLarger1 = abs(real(f_c_over_eig)) > 1.0;
    % fcLarge = zeros(size(f_c_over_eig));
    % fcSmall = zeros(size(f_c_over_eig));

    % fcLarge(idxMatLarger1) = f_c_over_eig(idxMatLarger1);
    % fcSmall(~idxMatLarger1) = f_c_over_eig(~idxMatLarger1);

    % delta = G * f12det_K * ( sum(sum(fcLarge .* (temp1 + temp2))) + sum(sum(fcSmall .* (temp1 + temp2))) );

    % temp1 = transpose(fwt_array .* multiplyByD( matmtrx, one_error_si_noD)) * one_error_si_noD;
    % temp2 = transpose(fwt_array .* multiplyByD( matmtrx, two_error_si_noD)) * two_error_si_noD;
    % delta = G * f1detmodij_col(1) * sum(sum( f_c_over_eig .* temp1)) + G * f2detmodij_col(1) * sum(sum(f_c_over_eig .* temp2));
    %==================================================
    cdetmodij_mat = ones(3, 1) * full_cdetmodij_col';  % duplicate so it has 3 entries for each sigma
    cdetmodij_mat = cdetmodij_mat(:) * ones(1, sdofc); % duplicate the column for sdof phi (sigma)
    temp = transpose(cwt_array .* cdetmodij_mat .* multiplyByD( D, full_csigma_noD)) * full_csigma_noD;
    e4t = sum(sum(c_c_over_eig .* temp));
    %==================================================
    temp1 = transpose(fwt_array .* f1detmodij_mat .* multiplyByD( D, one_fsigma_noD)) * half1_csigma_noD;
    temp2 = transpose(fwt_array .* f2detmodij_mat .* multiplyByD( D, two_fsigma_noD)) * half2_csigma_noD;
    e2t = sum(sum(fc_c_over_eig .* temp1)) + sum(sum(fc_c_over_eig .* temp2));

    temp1 = transpose(cwt_array .* multiplyByD( D, half1_csigma_noD)) * (one_fsigma_noD .* f1detmodij_mat);
    temp2 = transpose(cwt_array .* multiplyByD( D, half2_csigma_noD)) * (two_fsigma_noD .* f2detmodij_mat);
    e3t = sum(sum(cf_c_over_eig .* temp1)) + sum(sum(cf_c_over_eig .* temp2));
    %==================================================
    error_array(iel)= sqrt(real(-(e1t-e2t-e3t+e4t))); %LHS of eq 6.6
    errestop = errestop + real(-(e1t-e2t-e3t+e4t));
    errestbot = errestbot + real(-e1t);
    ET = [ET; e1t, e2t, e3t, e4t, real(-(e1t-e2t-e3t+e4t)), errestop, errestbot, real(-delta)];
    deltaerror(iel) = sqrt(real(-delta)); 
end

errest = 100*(errestop^0.5)/(errestbot^0.5);

% RHS of equation 6.5 in book
rhs65= (1/3)*max(deltaerror);
% RHS of equation 6.6 in book
% rhs66= 10*mean(deltaerror)
rhs66 = 10*mean(error_array);

%check to see which element need to be refined
array= zeros(nel,1);
for iel=1:nel
    if or(deltaerror(iel) > rhs65, error_array(iel) > rhs66)
        array(iel)=1;%if =1 need to be refined
    end
end

end