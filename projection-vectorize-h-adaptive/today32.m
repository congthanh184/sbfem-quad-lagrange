
function [E0,E1,E2,b, phi1, q1,phit,qt,globstifmatmtx]= sbfesidfaK(sidl,p,opt,nel,nnode,emodule,poisson,gcoord,centrecoord,nodes,sidfares,t,supelft)
%use for linear shape functions cirhole single unbounded domain
%-----------------------------------------------------
%Solving for global stiffness matrix of superelement or the whole domain (in case only one scaling 
%centre presented)
%units: MPa- N- mm
%-----------------------------------------------------
%input data for control parameters
nnel= 2;
for i=1:nel
    ndof(i,:)= [2 2*(p(i)-1) 2];
    edof(i)= ndof(i,1) + ndof(i,2) + ndof(i,3);
end
sdof= nnode*2 +2*(sum(p)-nel);
ngl= 64;
%-----------------------------------------------------
%input data for nodal coordinate values
%gcoord(i,j) where i-node no. and j-x or y coordinate values
%-----------------------------------------------------
%input data for nodal connectivity for each element
% nodes(i,j) where i-element no. and j-global nodes in that element
%-----------------------------------------------------
%initialization of matrices and vectors
globstifmatmtx= zeros(sdof,sdof);
b1= zeros(3,2);
b2= zeros(3,2);
%B1= zeros(3,edof);
%B2= zeros(3,edof);
E0= zeros(sdof,sdof);
E1= zeros(sdof,sdof);
E2= zeros(sdof,sdof);
matmtrx= zeros(3,3);
%-----------------------------------------------------
%compute element matrices and vector and assemble
[point1, weight1]= sbfeglqd1(ngl);
[matmtrx,G]=sbfematisolq(1,emodule,poisson); 

for iel= 1:nel                          % loop for total no. of element
    count=0;
    for i= 1:length(nodes(1,:))
        if nodes(iel,i)~=0
            count=count+1;
        end
    end
    local(iel).nd= zeros(1,count); 
    for i= 1:count
        local(iel).nd(i)= nodes(iel,i); 
        xcoord(i)= gcoord(local(iel).nd(i), 1);
        ycoord(i)= gcoord(local(iel).nd(i), 2);
    end
    
    localE0= zeros(edof(iel), edof(iel));           % for each element loop, local E0/E1/E2 startvalue= zeros(..,..)
    localE1= zeros(edof(iel), edof(iel)); 
    localE2= zeros(edof(iel), edof(iel)); 
    
    for int= 1:ngl
        nita= point1(int,1);%row vector
        wt= weight1(int,1);%row vector
        
        si=1;
        
        [N,shapelq,dhdnitalq]=lobatto(nita,p(iel));%row-vectors
        
        [modij, invmodij]=sbfemodijcirlipa(nnel,centrecoord(iel).coord,shapelq,dhdnitalq,xcoord,ycoord);%ok for this case cos' center=[0 0 0...0]   
        
        [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
        
        detj= det(modij);
        
        %compute the element matrices
        localE0= localE0+ B1'*matmtrx*B1*wt*detj;
        localE1= localE1+ B2'*matmtrx*B1*wt*detj;
        localE2= localE2+ B2'*matmtrx*B2*wt*detj;
      
    end
    
    [index]= sbfeeldof(local(iel).nd,2); %index of small element 3 nodes - 6 dofs
        
    [E0]= sbfeasmbl(E0,localE0,index);
    [E1]= sbfeasmbl(E1,localE1,index);
    [E2]= sbfeasmbl(E2,localE2,index);
    
end %end of 'iel=1:nel' loop, obtain E0,E1,E2

n=length(sidfares);% delete the corresponding 'sidfares' rows and cols from E0,E1,E2. Starting from the last 'sidfares' element
if n>0
   for i=n:-1:1
    c=sidfares(i);
    E0(:,c)=[]; 
    E0(c,:)=[];
    E1(:,c)=[]; 
    E1(c,:)=[];
    E2(:,c)=[]; 
    E2(c,:)=[];
    supelft(c)=[];
   end
else
end
% k=[inv(E0)*E1', -inv(E0); E1*inv(E0)*E1'-E2, -E1*inv(E0)];%reduced matrix or not depend on line90
k=[E0\E1', -inv(E0); E1*(E0\E1')-E2, -E1/E0];%reduced matrix or not depend on line90
% phitt=inv((t+1)^2*E0 + (t+1)*(E1'-E1)-E2)*(-supelft)*sidl/G;% Enow= Eactual/G--) inv(Eactual)= (1/G)*inv(Enow)
phitt=((t+1)^2*E0 + (t+1)*(E1'-E1)-E2)\(-supelft)*sidl/G;% Enow= Eactual/G--) inv(Eactual)= (1/G)*inv(Enow)
qtt=[(t+1)*E0 + E1']* G * phitt;%reduced matrix or not depend on line90
%-------------------------------------------------------------------------
% Compute V,D
for i=1:length(k) %trick1: shift theorem- add a const to diagonal of k to avoid giving 0s which leads to singular matrix
    k(i,i)= k(i,i)+2.7453;
end

[V,D] = eig(k);% V is the same; D(i,i)=D(i,i)+const

for i=1:length(D)
    D(i,i)= D(i,i)-2.7453;% actual D(i,i)
end

for i=((length(V))/2+1):length(V)
    for j=1:length(V)
        V(i,j)=V(i,j)*G;
    end
end
%-------------------------------------------------------------------------
% sort out phi1, q1 and estimate K by sub-function 
phi1= zeros(sdof,sdof);
q1= zeros(sdof,sdof);

a=zeros(1,length(D));
for i=1:length(D)
a(i)= D(i,i);
end
if opt==1
    [indices]=find(real(a)< -0.00001);%bounded
else
    [indices]=find(real(a)> 0.00001);%unbounded
end
%length(indices)= sdof - n - no. of rigid modes(2)
b= zeros(1, sdof);
bb= zeros(1, sdof-n);
for i=1:length(indices)
    bb(i)= a(indices(i));
end
for i=(length(indices)+1):(sdof-n)
     bb(i)= 0;
end

phi= zeros(sdof-n, sdof-n);
q= zeros(sdof-n, sdof-n);
for i= 1:(sdof-n)               %extract from V values for phi, q
    for j= 1:length(indices)
        phi(i,j)= V(i, indices(j));
        q(i,j)  = V((i+sdof-n), indices(j));
    end
end

for i= 1:(sdof-n)               %add rigid body motion in to fulfil the size 
    for j= (1+ length(indices)):(sdof-n)
        
            phi(i,j)=0;
            q(i,j)  = 0;
    end
end

temp=zeros(1,2*length(nodes(:,1)));
for k=1:length(nodes(:,1))
    temp(2*k-1)= local(k).nd(1);%nodes(k,1);
    temp(2*k)= local(k).nd(length(local(k).nd));
end
for i=length(temp):-1:2
    if temp(i)==temp(i-1)
        temp(i)=[];
    end
end
%index=sbfeeldof(temp,2);
if (1+ length(indices))<(sdof-n)
for k=1:length(temp)  
            phi(2*temp(k)-1,(1+ length(indices)):(sdof-n))=[1 0];
            phi(2*temp(k),(1+ length(indices)):(sdof-n))=[0 1];
end
else 
end
          

if n==0
    phi1=phi;
    q1=q;
    b=bb;
    phit=phitt;
    qt=qtt;
elseif n>0
    phit= zeros(sdof,1);
    qt= zeros(sdof,1); 
    redindex= [1:sdof];
    for i=n:-1:1
        c=sidfares(i);  %add 0s 1s to reach the size sdof and to account for the sidfa dof
        phi1(c,c)=1;
        q1(c,c)=0;
        redindex(:,c)=[];%need to define the reduced index for super-element
    end
    %assemble the sub matrices phi and q into phi1 and q1
    [phi1]= sbfeasmbl(phi1,phi,redindex);
    [q1]= sbfeasmbl(q1,q,redindex);
    %assemble the sub matrices redphit and redqt into phit and qt
    [phit]= sbfeasmblcol(phit,phitt,redindex);
    [qt]= sbfeasmblcol(qt,qtt,redindex);
    btp= b';
    [btp]= sbfeasmblcol(btp,bb',redindex);%btp = transpose
    b=btp';
end
if opt==1
        globstifmatmtx= real(q1/(phi1));
%         globstifmatmtx= real(q1 * pinv(phi1));
else
        globstifmatmtx= -real(q1/(phi1));
end


