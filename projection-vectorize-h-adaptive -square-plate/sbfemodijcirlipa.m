 function [modij, invmodij]=sbfemodij(nnel,centrecoord,shape,dhdnita,xcoord,ycoord)

%------------------------------------------------------------------------
%  Purpose:
%     determine the similarJacobian for two-dimensional mapping
%
%  Synopsis:
%     [modij, invmodij]=sbfemodij(nnel,centrecoord,shape,dhdnita,xcoord,ycoord) 
%
%  Variable Description:
%     similarjacob2 - modified Jacobian for two-dimensional mapping but for boundary only
%     nnel - number of nodes per element   
%     dhdnita - derivative of shape functions w.r.t. natural coordinate
%     nita
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%     centrecoord - cartesian coordinates for centre
%------------------------------------------------------------------------

  
 modij=zeros(2,2);
 
 for i=1:length(shape)
 modij(1,1)=modij(1,1)+shape(i)*(xcoord(i)-centrecoord(i,1)); % si=1
 modij(1,2)=modij(1,2)+shape(i)*(ycoord(i)-centrecoord(i,2)); % si=1
 modij(2,1)=modij(2,1)+dhdnita(i)*(xcoord(i)-centrecoord(i,1)); % si=1
 modij(2,2)=modij(2,2)+dhdnita(i)*(ycoord(i)-centrecoord(i,2)); % si=1
 end
 
invmodij= inv(modij);