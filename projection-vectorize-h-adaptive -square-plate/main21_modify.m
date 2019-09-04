%prog for using energy norm in stress as in paper

%start with 
%p=[1;1;1;1];

function [localcoarse,array,deltaerror,rhs65,error_array,rhs66,errest, ET, ccSave]= main21_modify(p, supel1, supel2)

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
[localfine,supelf,uf] = today3_modify(p + ones(size(p)), supel2);

pc = p; pf = p+ones(nel,1);
cc = supelc(1).c; cf = supelf(1).c;
eigvalc = supelc(1).values; eigvalf = supelf(1).values;
sdofc = nnode*2 + 2*(sum(pc)-nel);
sdoff = nnode*2 + 2*(sum(pf)-nel);

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
cf2c = inv(supelc(1).phi1)*uf(indexndf2c);
% cc = cf2c;

%uniform refined mesh all upgrade to p+1, using this as the correct solution to derived error for coarse mesh and define which to refine
ngl=16; %can go up to 10 N'*N*det(J) ~ 3*p
[p64,w64] = sbfeglqd1(ngl);
[matmtrx,G] = sbfematisolq(1,emodule,poisson);
D = G*matmtrx;
tol= 1;%tolerate error =1%

error_array = zeros(nel,1);% contribution of each ele in error eq 6.6
deltaerror = zeros(nel,1); % error rate for each ele eq 6.5
errestop = 0;
errestbot = 0;
errest = 0; %this the relative error norm to count for stop criteria 
ET = [];
for iel=1:nel

    e1t=0; e2t=0; e3t=0; e4t=0; %error1,2,3,4 to correspond to 4 sub integration
    E =zeros(sdoff-2,sdoff-2);
    delta=0;
    
    for i=1:(sdoff-2) %compute e1t    FINE mesh
        for j=1:(sdoff-2)
            for int=1:ngl
                nita = p64(int); 
                wt = w64(int);
                [N,shapelq,dhdnitalq,Nder]= lobatto(nita,pf(iel));%row-vectors
                [modij, invmodij]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelq,dhdnitalq,localfine(iel).xcrd,localfine(iel).ycrd);   
                [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
                
                e1i = cf(i) * D * (-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i) ;
                e1j = cf(j) * D * (-eigvalf(j) *B1 +B2) * localfine(iel).phi(:,j) ;
                
                %e1t = e1t + wt * (cf(i)*cf(j)/(eigvalf(i)+eigvalf(j))) * ((-eigvalf(i) *B1 +B2) * localfine(iel).phi(:,i))' * (G*matmtrx)' * ((-eigvalf(j) *B1 +B2) * localfine(iel).phi(:,j)) * det(modij);
                e1t = e1t + wt * (1/(eigvalf(i)+eigvalf(j))) * transpose(e1i) * inv(D) * e1j * det(modij);
                
                e5i = cf(i) * ((-eigvalf(i) * B1 + B2) * localfine(iel).coeffstr(:,i) - (-eigvalf(i) * B1 + B2) * localfine(iel).coefcstr(:,i));
                e5j = cf(j) * ((-eigvalf(j) * B1 + B2) * localfine(iel).coeffstr(:,j) - (-eigvalf(j) * B1 + B2) * localfine(iel).coefcstr(:,j));
                delta = delta + wt * (1/(eigvalf(i)+eigvalf(j))) * transpose(e5i) * (G*matmtrx)' * e5j * det(modij);
            end
        end    
    end
    ccSave = [];
    for i=1:(sdofc-2) %compute e4t    COARSE mesh
        for j=1:(sdofc-2)
            for int=1:ngl
                nita = p64(int); 
                wt = w64(int);
                [N,shapelq,dhdnitalq,Nder]= lobatto(nita,pc(iel));%row-vectors
                [modij, invmodij]=sbfemodijcirlipa(nnel,supelc(1).centrecoord(iel).coord,shapelq,dhdnitalq,localcoarse(iel).xcrd,localcoarse(iel).ycrd);   
                [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
                
                e4i = cc(i) * D * (-eigvalc(i) *B1 +B2) * localcoarse(iel).phi(:,i) ;
                e4j = cc(j) * D * (-eigvalc(j) *B1 +B2) * localcoarse(iel).phi(:,j) ;
                e4t = e4t + wt * (1/(eigvalc(i)+eigvalc(j))) * transpose(e4i) * inv(D) * e4j * det(modij);
                
            end
            ccSave(i, j) = cc(i)*cc(j) / (eigvalc(i)+eigvalc(j));
        end    
    end
    
    for i=1:(sdoff-2) %compute e2t & e3t    FINE&COARSE meshes
        for j=1:(sdofc-2)
            for int=1:ngl
                nita = p64(int); 
                wt = w64(int);
                [Nf,shapelqf,dhdnitalqf,Nderf]= lobatto(nita,pf(iel));%row-vectors
                [modijf, invmodijf]=sbfemodijcirlipa(nnel,supelf(1).centrecoord(iel).coord,shapelqf,dhdnitalqf,localfine(iel).xcrd,localfine(iel).ycrd);   
                [b1f, b2f, B1f, B2f]=sbfederiv2(nnel,shapelqf,dhdnitalqf,invmodijf);
                
                [Nc,shapelqc,dhdnitalqc,Nderc]= lobatto(nita,pc(iel));%row-vectors
                [modijc, invmodijc]=sbfemodijcirlipa(nnel,supelc(1).centrecoord(iel).coord,shapelqc,dhdnitalqc,localcoarse(iel).xcrd,localcoarse(iel).ycrd);   
                [b1c, b2c, B1c, B2c]=sbfederiv2(nnel,shapelqc,dhdnitalqc,invmodijc);
                
                e2t = e2t + wt * (cf(i)*cc(j)/(eigvalf(i)+eigvalc(j))) * transpose((-eigvalf(i) *B1f +B2f) * localfine(iel).phi(:,i)) * (G*matmtrx)' * ((-eigvalc(j) *B1c +B2c) * localcoarse(iel).phi(:,j)) * det(modijf);
                e3t = e3t + wt * (cf(i)*cc(j)/(eigvalf(i)+eigvalc(j))) * transpose((-eigvalc(j) *B1c +B2c) * localcoarse(iel).phi(:,j)) * (G*matmtrx)' * ((-eigvalf(i) *B1f +B2f) * localfine(iel).phi(:,i)) * det(modijf);
            end
        end    
    end
    ET = [ET; e1t e2t e3t e4t];
    error_array(iel)= sqrt(real(-(e1t-e2t-e3t+e4t))); %LHS of eq 6.6
    errestop = errestop + real(-(e1t-e2t-e3t+e4t));
    errestbot = errestbot + real(-e1t);
    
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