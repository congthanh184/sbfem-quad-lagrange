function [local,supel,globnodesdisp]= refsolpa_platesqrplobat(p, supel)
% need to have p= input value of this func
%18/9
%-----------------------------------------------------
%using sbfe to solve an infinite plate with 1unit udl load
%Plane stress problem. Quadratic line elements to account for curve
%units: MPa- N- mm
%bounded size: 3000*3000 with hole 1000*1000
%load: P= 1000(N)/1000mm
%---------------------------------------------------
nsupel = 1;
t=0.0001;
ff = supel.ff;
bcdof = supel.bcdof;
bcval = supel.bcval;
nnel = supel.nnel;
ndof = supel.ndof;
%-----------------------------------------------------
% Compute Ks for each superelement and assemble the sub matrices Ks into the whole domain stif matrx K
ff_length = length(supel.ft);
K = zeros(ff_length, ff_length);
% K=zeros(length(ff),length(ff));% sdof,sdof

for iel= 1:nsupel                         % loop for total no. of element
    
    [supel(iel).E0,supel(iel).E1,supel(iel).E2,supel(iel).values,supel(iel).phi1,supel(iel).q1,supel(iel).phit,supel(iel).qt,supel(iel).Ks]= today32(supel(iel).sidl,p,supel(iel).opt,supel(iel).nel,supel(iel).nnode,supel(iel).emodule,supel(iel).poisson,supel(iel).gcoord,supel(iel).centrecoord,supel(iel).localnodes,supel(iel).sidfares,t,supel(iel).ft);
    %E0,E1,E2 in reduced size, phi1,q1,Ks... in normal size
    
    nd= supel(iel).nodes; %now consider the superelement as single ele and assign nodes for that ele
       
    [supel(iel).index]= sbfeeldof(nd,2);%the element now not the basics 3nodes-6dofs but it is the supel(iel)
    
    [K]= sbfeasmbl(K,supel(iel).Ks,supel(iel).index);
    
end

%%
%--------------------------------------------------------------------------
% Apply boundary conditions. Solve the conventional finite element equation Ku= P for u. Define the
% displacements and stresses for the whole domain

[K,ff]=sbfeaplyc2(K,ff,bcdof,bcval);

globnodesdisp= K\ff;% depending on p but without derivatives
% globnodesdisp= pinv(K)*ff;
% globnodesdisp(bcdof) = 0;
%-------------------------------------------------------------------------

for iel= 1:nsupel   %extract globnodesdisp and c for each superelement
    supel(iel).globnodesdisp= zeros(length(supel(iel).index), 1);
    for i= 1:length(supel(iel).index)
    supel(iel).globnodesdisp(i)= globnodesdisp(supel(iel).index(i));
    end
    supel(iel).c = (supel(iel).phi1)\ supel(iel).globnodesdisp;
%     supel(iel).c = pinv(supel(iel).phi1) * supel(iel).globnodesdisp;
end
%cirhole1plot1== cirhole1plot to plot for the bounded 
%iel=1;
% [local]= today33(p,nnel,ndof,supel(iel).nnode,supel(iel).nel,supel(iel).localnodes,supel(iel).gcoord,supel(iel).phi1,supel(iel).centrecoord,supel(iel).c,supel(iel).values,supel(iel).emodule,supel(iel).poisson);
[local]= today33_vectorize(p,nnel,ndof,supel(iel).nnode,supel(iel).nel,supel(iel).localnodes,supel(iel).gcoord,supel(iel).phi1,supel(iel).centrecoord,supel(iel).c,supel(iel).values,supel(iel).emodule,supel(iel).poisson);
%cirhole1plot2 to plot the unbounded domain, domain restrain: x>3000, y>3000 

%-------------------------------------------------------------------------
%calculate projection based interpolation solution (displacement)
