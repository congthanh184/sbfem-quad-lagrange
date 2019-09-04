function [matmtrx,G]=sbfematisolq(iopt,elastic,poisson)

%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation/ material matrix for isotropic material
%
%  Synopsis:
%     [matmtrx]=sbfematisolq(iopt,elastic,poisson) 
%
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis
%------------------------------------------------------------------------

G=elastic/2/(1+poisson);

%in order to balance the orders of elements in matrix k to yield V displacements parts 
%different from 0s, have to divide the original matmtrx by G, see file notdivG for 
%matrices V,D when do not take matmtrx/G. The D matrix is still the same but the V's 
%elements are mainly of 0

%after obtain V,D matrices, have to adjust V to obtain the actual nodal displacements and forces 

 if iopt==1        % plane stress
   matmtrx= elastic/G/(1-poisson*poisson)* ...
   [1  poisson 0; ...
   poisson  1  0; ...
   0  0  (1-poisson)/2];

 elseif iopt==2     % plane strain
   matmtrx= elastic/G/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson 0; 
   poisson  (1-poisson)  0;
   0  0  (1-2*poisson)/2];

 elseif iopt==3     % axisymmetry
   matmtrx= elastic/G/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
 
 else     % three-dimension
   matmtrx= elastic/G/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];

 end
