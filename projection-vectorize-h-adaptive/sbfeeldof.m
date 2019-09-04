function [index]= sbfeeldof(nd, ndof)
%-------------------------------------------------------------
%Purpose:
%    Compute system dofs associated for each element 
%
%Synopis:
%    [index]= sbfeeldof(nd, nnel, ndof)
%
%Variable description:
%    index - system dof vector assocciated with element iel
%    iel - element number whose system dofs are to be determined
%    nnel - number of nodes per element
%    ndof - number of dofs per node
%-------------------------------------------------------------

%edof= nnel*ndof;
k=0;

  for i=1:length(nd)
    start = (nd(i)-1)*ndof;
    for j=1:ndof
        k= k+1;
        index(k)= start+j;
    end
  end
