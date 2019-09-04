
function [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,modij)

%------------------------------------------------------------------------
%  Purpose:
%     determine components B1 and B2 for the SBFE kinematic matrix
%
%  Synopsis:
%     [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,modij)  
%
%  Variable Description:
%     b1, b2 - b1=L1R11+L2R21; b2=L1R12+L2R22; 
%     B1, B2 - B1= b1(s)*N(s); B2= b2(s)* (N(s)',s) 
%     nnel - number of nodes per element   
%     shapelq - shape functions for 3 nodes line element 
%     dhdnita - derivative of shape functions w.r.t. natural coordinate
%     nita
%     invjacob - inverse of 2-D Jacobian at boundary 
%------------------------------------------------------------------------

L1=[1 0; 0 0; 0 1];
L2=[0 0; 0 1; 1 0];

      %invjacob at boundary, hence dhdsi& dhdnita depending only on nita
 b1 =L1*modij(1,1)+ L2*modij(2,1);
 b2 =L1*modij(1,2)+ L2*modij(2,2);   
 
 for i=1:length(shapelq)
     
     i1= (i-1)*2+1;
     i2= i1+1;
     
     B1(1,i1)= modij(1,1)*shapelq(i);
     B1(2,i2)= modij(2,1)*shapelq(i);
     B1(3,i1)= modij(2,1)*shapelq(i);
     B1(3,i2)= modij(1,1)*shapelq(i);
     
     B2(1,i1)= modij(1,2)*dhdnitalq(i);
     B2(2,i2)= modij(2,2)*dhdnitalq(i);
     B2(3,i1)= modij(2,2)*dhdnitalq(i);
     B2(3,i2)= modij(1,2)*dhdnitalq(i);
     
 end
