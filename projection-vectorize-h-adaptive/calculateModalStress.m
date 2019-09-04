function [x, y, modalStress] = calculateModalStress(D, iel, supel, local, p, posLob)
	modalStress = [];
	nnel = supel.nnel;
	nnode = supel.nnode;
	nel = supel.nel;

	sdoff = nnode*2 + 2*(sum(p)-nel);
	eigvalf = supel.values;

	[N,shapelq,dhdnitalq,Nder]= lobatto(posLob,p(iel));%row-vectors
	[modij, invmodij]=sbfemodijcirlipa(nnel,supel(1).centrecoord(iel).coord,shapelq,dhdnitalq,local(iel).xcrd,local(iel).ycrd);   
	[b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);					
	x = shapelq * local(iel).xcrd';
	y = shapelq * local(iel).ycrd';
	
	for idoff = 1:sdoff
		modalStress = [modalStress,  D * (-eigvalf(idoff) *B1 +B2) * local(iel).phi(:,idoff)] ;
	end

end