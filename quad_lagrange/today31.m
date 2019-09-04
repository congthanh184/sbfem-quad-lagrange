    function [nodes1,gcoord1,elecenter,ffval]= trans(p,orimesh,s,supelcentre,nel)
% transform for cirhole, linear shape, single unbounded 
%calculate the constants for only 1 unbounded cirhole domain using the
%spectral order shape functions from paper
%generate temps to create/ nodes1
countemps=length(orimesh(:,1));
for iel=1:nel
    local(iel).temps=zeros(1,s(iel,2)+1);
    local(iel).temps(1)=s(iel,3);
    local(iel).temps(s(iel,2)+1)=s(iel,4);
    for j=2:s(iel,2)
        countemps=countemps+1;
        local(iel).temps(j)=countemps;
    end
end
le = 0;
for i=1:nel
    le=le+length(local(i).temps);
end
nodes1=zeros(1,le);
count=0;
for i=1:nel
    temp=length(local(i).temps);
    for j=1:temp
        count=count+1;
        nodes1(count)=local(i).temps(j);
    end
end
for i=le:-1:2
    if nodes1(i)==nodes1(i-1)
        nodes1(i)=[];
    end
end
gcoord1=zeros(length(nodes1),2);
for i=1:length(orimesh(:,1))
    for j=1:length(gcoord1(:,1))
        if nodes1(j)==orimesh(i,1)
            gcoord1(j,:)=orimesh(i,2:3);
        end
    end
end
for iel=1:nel            
ndof(iel,:)=[2 2*(p(iel)-1) 2];
end
%1st col is ith element, 2nd col is p value for that ele
%centre
qy=-1; %1N/mm pointing in y negative direction

for i=1:nel
    elecenter(i).coord=zeros(p(i)+1,2);% later if p different need to take p-1 for each ele then sum
    elecenter(i).coord(1,:)= supelcentre;
    elecenter(i).coord(p(i)+1,:)= supelcentre;
    if p(i)>1
        for j=2:p(i)
            elecenter(i).coord(j,:)= [0 0];
        end
    end
end
      
for i=1:nel
    [~, nd] = ismember(s(i,3:4), orimesh(:, 1));
    local(i).gcoord = zeros(1+s(i,2),2);%coord = ele actual real xycoord
    local(i).gcoord(1,:)=orimesh(nd(1),2:3);
    local(i).gcoord(length(local(i).gcoord(:,1)),:)=orimesh(nd(2),2:3);
end

ffval= 0;
ngl= 64;
[point1, weight1]= sbfeglqd1(ngl);

% distributed force between node 1 and 2
idx_node1 = find(s(:,3)==2);
idx_node2 = find(s(:,4)==3);

for iel=idx_node1:idx_node2
    ffval_temp = 0;
    for int=1:ngl
        nita= point1(int,1);
        wt= weight1(int,1);
        % [NN,local(iel).N,local(iel).dhdnita]=lobatto(nita,p(iel));
        [NN,local(iel).N,dNN, local(iel).dhdnita]=lagrangian2(nita,p(iel));
        dyds= 0;
%         for j=1:length(local(iel).N)
%             dyds= dyds+local(iel).dhdnita(j)*(local(iel).gcoord(j,2)-elecenter(iel).coord(j,2)); % si=1
%         end
        dyds = local(iel).dhdnita*(local(iel).gcoord(:, 2) - elecenter(iel).coord(:, 2));
        ffval_temp = ffval_temp + (local(iel).N)'*qy*abs(dyds)*wt;%abs cos' the distance have to be positive
    end
    if ffval == 0
        ffval = ffval_temp;
    else
        ffval = [ffval(1:end-1); ffval(end) + ffval_temp(1); ffval_temp(2:end)];
    end
end


 
