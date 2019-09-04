function [local,supel,globnodesdisp]= refsolpa_platesqrplobat(p)
% need to have p= input value of this func
%18/9
%-----------------------------------------------------
%using sbfe to solve an infinite plate with 1unit udl load
%Plane stress problem. Quadratic line elements to account for curve
%units: MPa- N- mm
%bounded size: 3000*3000 with hole 1000*1000
%load: P= 1000(N)/1000mm
%---------------------------------------------------
%input data for control parameters
nsupel=1; 
nnel=2;
nel= [4];
nnode= [5];
opt1=1;%bounded
emodule= [250000]; % unit= MPa
poisson= [0.3];

orimesh= [1 1000 2000;...
          2 0 2000;...
          3 0 0;...
          4 2000 0;...
          5 2000 1000];
centrecoord=[1000 1000];
% p = [1;1;1;1];
%-----------------------------------------------------
%p=5;
for iel=1:nel
    ndof(iel,:)=[2 2*(p(iel)-1) 2];
end
s= [1 p(1) 1 2;...
    2 p(2) 2 3;...
    3 p(3) 3 4;...
    4 p(4) 4 5];

%input data for nodal coordinate values
%gcoord(i,j) where i-node no. and j-x or y coordinate values
%----------------------------------------------------- 
[nodes1,gcoord1,elecenter,ffval]= today31(p,orimesh,s,centrecoord,nel(1));
centrecoord1=elecenter;
%localcenter=center for all ele in the supel, later if p different need to
%adjust, ffval in local order 1:1:p+1,need adjust to put in ff
%-----------------------------------------------------
%define pltx,plty
%-----------------------------------------------------
% input data for nodal connectivity for each element
% localnodes(i,j) where i-element no. and j- superelement global nodes in
% that element
for iel=1:nel
    localnodes1=zeros(nel(1),max(p)+1);
end
start=1;
for i=1:nel
    temp= start:1:(start+p(i));
    localnodes1(i,:)= [start:1:(start+p(i)),zeros(1,max(p)+1-length(temp))];
    start= start+p(i);
end
%-----------------------------------------------------
%input data for constructing vector index for each superelement 
%-----------------------------------------------------
%input data for boundary conditions
for i=1:length(nodes1)
    if nodes1(i)==1
        posi1(1)=i;
    elseif nodes1(i)==2
        posi1(2)=i;
    end
end
for i=1:length(nodes1)
    if nodes1(i)==4
        posi4(1)=i;
    elseif nodes1(i)==5
        posi4(2)=i;
    end
end
len1=length(posi1(1):1:posi1(2));
len4=length(posi4(1):1:posi4(2));
nd1=zeros(1,len1);nd4=zeros(1,len4);

    nd1= nodes1(posi1(1):1:posi1(2));
    bcdof1=2*nd1;
    
    nd4= nodes1(posi4(1):1:posi4(2));
    bcdof4=2*nd4-1;

bcdof=[bcdof1,bcdof4];
bcval= zeros(1,length(bcdof));

%force vector, need to assemble ffval in 
total= sum(p)-nel;
ff= zeros(nnode(1)*2+2*total,1);
for i=1:length(nodes1)
    if nodes1(i)==2
        posi2(1)=i;
    elseif nodes1(i)==3
        posi2(2)=i;
    end
end

nd2= nodes1(posi2(1):1:posi2(2));
for i=1:1:length(nd2)
    ff(2*nd2(i)-1)=ffval(i);
end

t=0.0001;
ft1= zeros(length(ff),1);
sidfares1= [];
l1=0;
%%
%------------------------------------------------------------------

supel=struct('nel',{nel(1)},'nnode',{nnode(1)},'gcoord',{gcoord1},...
    'centrecoord',{centrecoord1},'localnodes',{localnodes1},'nodes',...
    {nodes1},'opt',{opt1},'emodule',{emodule(1)},'poisson',{poisson(1)},...
    'ft',{ft1},'sidfares',{sidfares1},'sidl',{l1});
%-----------------------------------------------------
% Compute Ks for each superelement and assemble the sub matrices Ks into the whole domain stif matrx K
K=zeros(length(ff),length(ff));% sdof,sdof

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

%-------------------------------------------------------------------------

for iel= 1:nsupel   %extract globnodesdisp and c for each superelement
    supel(iel).globnodesdisp= zeros(length(supel(iel).index), 1);
    for i= 1:length(supel(iel).index)
    supel(iel).globnodesdisp(i)= globnodesdisp(supel(iel).index(i));
    end
    supel(iel).c= inv(supel(iel).phi1)* supel(iel).globnodesdisp;
end
%cirhole1plot1== cirhole1plot to plot for the bounded 
%iel=1;
[local]= today33(p,nnel,ndof,supel(iel).nnode,supel(iel).nel,supel(iel).localnodes,supel(iel).gcoord,supel(iel).phi1,supel(iel).centrecoord,supel(iel).c,supel(iel).values,supel(iel).emodule,supel(iel).poisson);
%cirhole1plot2 to plot the unbounded domain, domain restrain: x>3000, y>3000 

%-------------------------------------------------------------------------
%calculate projection based interpolation solution (displacement)
