% clear;
clc;

orimesh = [ 1 1000 2000; ...    
    2 0 2000; ...
    3 0 0; ...
    4 2000 0; ...
    5 2000 1000];

centrecoord=[1000 1000];

[matmtrx,G] = sbfematisolq( 1, 250000, 0.3);
% D = G*matmtrx;
D = matmtrx;

[pt, wt] = sbfeglqd1(64);
sdof = (size(orimesh, 1) + 4 + 2*4) * 2;

E0 = zeros(sdof, sdof);
E1 = zeros(sdof, sdof);
E2 = zeros(sdof, sdof);


p = ones( size(orimesh, 1)-1, 1) * 2;
[supel] = initialize_problem(orimesh, centrecoord, p);

for idRow = 2:size(orimesh, 1) 
	x_m1 	= orimesh(idRow-1, 2); 	y_m1 	= orimesh(idRow-1, 3);
	x_1 	= orimesh(idRow, 2); 	y_1 	= orimesh(idRow, 3);
	x_mid	= mean([x_m1, x_1]);	y_mid 	= mean([y_m1, y_1]);

	p_m1 	= calcXs(x_m1, y_m1, centrecoord(1), centrecoord(2), 1);
	p_1 	= calcXs(x_1, y_1, centrecoord(1), centrecoord(2), 1);
	p_mid 	= calcXs(x_mid, y_mid, centrecoord(1), centrecoord(2), 1);

	virtualPoint = [0;0];

	coordParent = [p_m1, virtualPoint, p_mid, virtualPoint, p_1];
	L1 = [1 0; 0 0; 0 1];
	L2 = [0 0; 0 1; 1 0];
	L = [L1 L2];

	localE0 = 0;
	localE1 = 0;
	localE2 = 0;

	for idPt = 1:numel(pt)
		[N, NN, dN, dNN] = quadLagrangHierarchial(pt(idPt), 1, 1);
		coordChildS = coordParent * NN';
		dCoordChildS = coordParent * dNN';
		Jacob = [coordChildS'; dCoordChildS'];
		detJ = det(Jacob);
		invJacob = Jacob \ eye(size(Jacob));
		b1 = L * convertSingleArrayToEyeMatrix(invJacob(:, 1), 2);
		b2 = L * convertSingleArrayToEyeMatrix(invJacob(:, 2), 2);
		B1 = b1 * N;
        B2 = b2 * dN;
        localE0 = localE0 + B1'*D*B1*detJ*wt(idPt);
        localE1 = localE1 + B2'*D*B1*detJ*wt(idPt);
        localE2 = localE2 + B2'*D*B2*detJ*wt(idPt);
	end

	idEBlock = idRow*8-6; idSBlock = idEBlock - size(localE0, 1) + 1;		% idEBlock = 2*3-2, 2*4-2, 2*5-2
	E0(idSBlock:idEBlock, idSBlock:idEBlock) = E0(idSBlock:idEBlock, idSBlock:idEBlock) + localE0;
	E1(idSBlock:idEBlock, idSBlock:idEBlock) = E1(idSBlock:idEBlock, idSBlock:idEBlock) + localE1;
	E2(idSBlock:idEBlock, idSBlock:idEBlock) = E2(idSBlock:idEBlock, idSBlock:idEBlock) + localE2;
end

invE0 = E0 \ eye(size(E0));

% solve 
k = [ invE0*E1', -invE0; E1*invE0*E1' - E2, -E1*invE0 ];
[eigVec, eigVal] = eig(k);

displacementModes = diag(eigVal);
indexBounded = real(displacementModes) < -0.1;
phi = zeros(sdof, sdof);
q = zeros(sdof, sdof);

phi(:, 1:sum(indexBounded==1)) = eigVec(1:sdof, indexBounded);
temp = repmat(eye(2), [sdof/2 1]);
phi(:, end-1:end) = temp;
phi(3:4:end, end-1:end) = 0;
phi(4:4:end, end-1:end) = 0;

q(:, 1:sum(indexBounded==1)) = eigVec( sdof+1:end, indexBounded) * G;	

globalStiffness = real(q/phi);
%=======================================================================================

ffval= 0;
ngl= 64;
[point1, weight1]= sbfeglqd1(ngl);
qy = -1;

for idRow = 3
    ffval_temp = 0;
    x_m1 	= orimesh(idRow-1, 2); 	y_m1 	= orimesh(idRow-1, 3);
	x_1 	= orimesh(idRow, 2); 	y_1 	= orimesh(idRow, 3);
	x_mid	= mean([x_m1, x_1]);	y_mid 	= mean([y_m1, y_1]);

	p_m1 	= calcXs(x_m1, y_m1, centrecoord(1), centrecoord(2), 1);
	p_1 	= calcXs(x_1, y_1, centrecoord(1), centrecoord(2), 1);
	p_mid 	= calcXs(x_mid, y_mid, centrecoord(1), centrecoord(2), 1);

	virtualPoint = [0;0];

	coordParent = [p_m1, virtualPoint, p_mid, virtualPoint, p_1];
    for int=1:ngl
        nita= point1(int,1);
        wt= weight1(int,1);
        % [NN,local(iel).N,local(iel).dhdnita]=lobatto(nita,p(iel));
        % [NN,local(iel).N,dNN, local(iel).dhdnita]=lagrangian2(nita,p(iel));
        [N, NN, dN, dNN] = quadLagrangHierarchial(nita, 1, 1);
        dyds= 0;
%         for j=1:length(local(iel).N)
%             dyds= dyds+local(iel).dhdnita(j)*(local(iel).gcoord(j,2)-elecenter(iel).coord(j,2)); % si=1
%         end
        % dyds = local(iel).dhdnita*(local(iel).gcoord(:, 2) - elecenter(iel).coord(:, 2));
        dyds = dNN * coordParent(2,:)';
        ffval_temp = ffval_temp + (NN)'*qy*abs(dyds)*wt;%abs cos' the distance have to be positive
    end
    if ffval == 0
        ffval = ffval_temp;
    else
        ffval = [ffval(1:end-1); ffval(end) + ffval_temp(1); ffval_temp(2:end)];
    end
end

elem1bcdof = 2:2:10;
elem4bcdof = 25:2:34;
bcdof = [elem1bcdof elem4bcdof];
ff = zeros(size(globalStiffness, 1), 1);
ff([9:2:18]) = ffval;

[globalStiffness,ff]=sbfeaplyc2(globalStiffness,ff,bcdof,zeros(size(bcdof)));
globalDisplacement = globalStiffness \ ff;
c = phi \ globalDisplacement;

%=======================================================================================
nitanew = [-1:0.05:1];
sinew = [1:-0.05:0];

if min(sinew)~= 0 
    sinew = [sinew 0];
end

local = [];
lambdaVal = zeros(numel(c), 1);
lambdaVal(1:sum(indexBounded)) = displacementModes(indexBounded);

% comparedata = load('today3_1111.mat');

for idRow = 2:size(orimesh, 1)
	x_m1 	= orimesh(idRow-1, 2); 	y_m1 	= orimesh(idRow-1, 3);
	x_1 	= orimesh(idRow, 2); 	y_1 	= orimesh(idRow, 3);
	x_mid	= mean([x_m1, x_1]);	y_mid 	= mean([y_m1, y_1]);

	p_m1 	= calcXs(x_m1, y_m1, centrecoord(1), centrecoord(2), 1);
	p_1 	= calcXs(x_1, y_1, centrecoord(1), centrecoord(2), 1);
	p_mid 	= calcXs(x_mid, y_mid, centrecoord(1), centrecoord(2), 1);

	virtualPoint = [0;0];

	coordParent = [p_m1, virtualPoint, p_mid, virtualPoint, p_1];
	% coordParent = [p_m1, p_1];
	L1 = [1 0; 0 0; 0 1];
	L2 = [0 0; 0 1; 1 0];
	L = [L1 L2];

	local(idRow-1).x = [];
	local(idRow-1).y = [];

	local(idRow-1).dispx = [];
	local(idRow-1).dispy = [];
	
	local(idRow-1).strx = [];
	local(idRow-1).stry = [];

	% idEBlock = idRow*8-6; idSBlock = idEBlock - size(localE0, 1) + 1;		
	idER = idRow*8-6; idSR = idER - size(localE0, 1) + 1; 
	localphi = phi(idSR:idER, :);
	% expectphi = comparedata.local(idRow-1).phi;
	% find(abs(localphi-expectphi) > 10e-5);

	for i = 1:numel(sinew)
		for j = 1:numel(nitanew)
			sii = sinew(i);
			nitaj = nitanew(j);
			[N, NN, dN, dNN] = quadLagrangHierarchial(nitaj, 1, 1);
			% [N, NN, dN, dNN] = lagrangian2(nitaj, 2);
			coordChildS = coordParent * NN';
			dCoordChildS = coordParent * dNN';
			Jacob = [coordChildS'; dCoordChildS'];
			detJ = det(Jacob);
			invJacob = Jacob \ eye(size(Jacob));
			b1 = L * convertSingleArrayToEyeMatrix(invJacob(:, 1), 2);
			b2 = L * convertSingleArrayToEyeMatrix(invJacob(:, 2), 2);
			B1 = b1 * N;
	        B2 = b2 * dN;

	        coordxy = centrecoord' + sii * coordParent * NN';
	        local(idRow-1).x(i, j) = coordxy(1);
	        local(idRow-1).y(i, j) = coordxy(2);

	        cmat = repmat(c, [1, 10])';
	        lmat = repmat(lambdaVal, [1, 10])';
	        pointdisp = real(N * sum(cmat .* sii.^(-lmat) .* localphi, 2));
	        local(idRow-1).dispx(i, j) = pointdisp(1);
	        local(idRow-1).dispy(i, j) = pointdisp(2);
	       
	       	pstr = 0;
	       	for jj = 1:sdof
	       		pstr = real(pstr + G* matmtrx * c(jj) * sii^(-lambdaVal(jj)-1) * (-lambdaVal(jj) * B1 + B2) * localphi(:, jj));
	       	end 
	        local(idRow-1).strx(i, j) = pstr(1);
	        local(idRow-1).stry(i, j) = pstr(2);

		end
	end
end

% plot_disp(local)
plot_stress(local, local)