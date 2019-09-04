function [supel] = initialize_problem(orimesh, centrecoord, p)
	%input data for control parameters
	nsupel=1; 
	nnel=2;
	opt1=1;%bounded
	emodule= [250000]; % unit= MPa
	poisson= [0.3];

	nnode = [size(orimesh, 1)];
	nel = [nnode-1];

	% tempCoord = linspace(1000, 0, 5);
	% horiElem = [tempCoord; 1000*ones(size(tempCoord))];
	% vertElem = [1000*ones(size(tempCoord)); tempCoord];
	% orimesh = [horiElem'; vertElem'];
	% centrecoord=[0 0];
 %    p = ones( size(orimesh, 1)-1, 1);

	% calculate ndof = [2 sum_of_addtional_p 2]
	temp_ones = ones(nel, 1);
	ndof = [2*temp_ones, 2*(p(:)-1),  2*temp_ones];

	% s= [1 p(1) 1 2;...
	%     2 p(2) 2 3;...
	%     3 p(3) 3 4;...
	%     4 p(4) 4 5];

	temp_1_nel = orimesh(:,1);
	s = [temp_1_nel(1:end-1), p, temp_1_nel(1:end-1), temp_1_nel(2:end)];

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
	% for iel=1:nel
	%     localnodes1=zeros(nel(1),max(p)+1);
	% end
	localnodes1 = zeros(nel, max(p)+1);
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

	% for i=1:length(nodes1)
	%     if nodes1(i)==1
	%         posi1(1)=i;
	%     elseif nodes1(i)==2
	%         posi1(2)=i;
	%     end
	% end

% constrain x on node1, constrain y on node 9
	[~, posi1] = ismember([1 2], nodes1);

	% for i=1:length(nodes1)
	%     if nodes1(i)==4
	%         posi4(1)=i;
	%     elseif nodes1(i)==5
	%         posi4(2)=i;
	%     end
	% end

	[~, posi4] = ismember([4 5], nodes1);

	% len1=length(posi1(1):1:posi1(2));
	% len4=length(posi4(1):1:posi4(2));
	% nd1=zeros(1,len1);nd4=zeros(1,len4);

	% constrain x 
	nd1= nodes1(posi1(1):1:posi1(2));
	bcdof1=2*nd1;
	    
	% constrain y
	nd4= nodes1(posi4(1):1:posi4(2));
	bcdof4=2*nd4-1;

	bcdof=[bcdof1,bcdof4];
	bcval= zeros(size(bcdof));

	%force vector, need to assemble ffval in 
	total= sum(p)-nel;
	ff= zeros(nnode(1)*2+2*total,1);
	% for i=1:length(nodes1)
	%     if nodes1(i)==2
	%         posi2(1)=i;
	%     elseif nodes1(i)==3
	%         posi2(2)=i;
	%     end
	% end

	[~, posi2] = ismember([2 3], nodes1);

	nd2= nodes1(posi2(1):1:posi2(2));
	for i=1:1:length(nd2)
	    ff(2*nd2(i)-1)=ffval(i);
	end

	ft1= zeros(length(ff),1);
	sidfares1= [];
	l1=0;

	%%
	%------------------------------------------------------------------
	supel=struct('nel',{nel(1)},'nnode',{nnode(1)},'gcoord',{gcoord1},...
	    'centrecoord',{centrecoord1},'localnodes',{localnodes1},'nodes',...
	    {nodes1},'opt',{opt1},'emodule',{emodule(1)},'poisson',{poisson(1)},...
	    'ft',{ft1},'sidfares',{sidfares1},'sidl',{l1}, ...
	    'ff', {ff}, 'bcdof', bcdof, 'bcval', bcval, 'nnel', nnel, 'ndof', ndof);	
end