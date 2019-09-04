function [local]= today33_vectorize(p,nnel,ndof,nnode,nel,nodes,gcoord,supelphi1,centrecoord,supelc,supelvalues,emodule,poisson)

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
         local(iel).xcrd(i)= gcoord(local(iel).nd(i), 1);
         local(iel).ycrd(i)= gcoord(local(iel).nd(i), 2); 
    end
end

nitanew=[-1:0.05:1];
sinew=[1:-0.01:0];%later should have bounded and unbounded as input, and if need should have different interval for dif ele to increase no of points for certain ele

for i=1:nel
    ndof(i,:)= [2 2*(p(i)-1) 2];
    edof(i)= ndof(i,1) + ndof(i,2) + ndof(i,3);
end
sdof= nnode*2 +2*(sum(p)-nel);

%index= zeros(edof,1);
b1= zeros(3,2);
b2= zeros(3,2);
%B1= zeros(3,edof);
%B2= zeros(3,edof);
%phi= zeros(edof,sdof);
%--------------------------------------------------------------------------

[matmtrx,G]=sbfematisolq(1,emodule,poisson); 
for iel=1:nel
    local(iel).phi= zeros(edof(iel),sdof); 
    local(iel).x=zeros(length(sinew),length(nitanew));
    local(iel).y=zeros(length(sinew),length(nitanew));
    local(iel).dispx=zeros(length(sinew),length(nitanew));
    local(iel).dispy=zeros(length(sinew),length(nitanew));
end


for iel= 1:nel                          % loop for total no. of element

    % for i= 1:length(local(iel).nd)
    %     xcoord(i)= gcoord(local(iel).nd(i), 1);
    %     ycoord(i)= gcoord(local(iel).nd(i), 2);
    % end

    xycoord = gcoord( local(iel).nd, :); % [x1 y1; x2 y2; ... ; xn yn]
    xcoord = xycoord(:, 1);
    ycoord = xycoord(:, 2);

    [index]= sbfeeldof( local(iel).nd, 2);
    % for i= 1:edof(iel)                       % using index to extract phi(edof*sdof) from phi1
    %     for j= 1:sdof
    %     local(iel).phi(i,j)= supelphi1(index(i),j);
    %     end
    % end
    local(iel).phi = supelphi1(index, :); % need comment for this

    coeffinestr = local(iel).phi;
    coefcoarsestr = [local(iel).phi(1:(edof(iel)-4),:);zeros(2,sdof);local(iel).phi((edof(iel)-1):edof(iel),:)];
    local(iel).coeffstr = coeffinestr;
    local(iel).coefcstr = coefcoarsestr;

    % set-up shape function matrix 
    N_mat = [];
    B1_mat = []; B2_mat = [];
    shapelq_mat = [];

    for i = 1:1:numel(nitanew)
        [N,shapelq,dhdnitalq,Nder]= lobatto(nitanew(i),p(iel));%row-vectors
        [modij, invmodij]=sbfemodijcirlipa(nnel,centrecoord(iel).coord,shapelq,dhdnitalq,xcoord,ycoord);   
        [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
        N_mat = [N_mat; N];
        shapelq_mat = [shapelq_mat; shapelq];
        B1_mat = [B1_mat; B1];
        B2_mat = [B2_mat; B2];
    end

    % store the shape matrix in each element
    local(iel).lobatto_N_mat = N_mat;
    local(iel).lobatto_shapelq_mat = shapelq_mat;
    local(iel).lobatto_B1_mat = B1_mat;
    local(iel).lobatto_B2_mat = B2_mat;    

    % better store sinew and nitanew to notify the size of shape matrix store in each element
    local(iel).sinew = sinew;
    local(iel).nitanew = nitanew;

    % calculate coordinator of each node in element = center origin + si * N * (xy - center)
    deltaCoord = xycoord - centrecoord(iel).coord;
    localXY_sii_1 = shapelq_mat * deltaCoord;

    % generate matrix of pointx = (i, j) => i is sii, j is nita
    local(iel).x = centrecoord(iel).coord(1, 1) + sinew(:) * localXY_sii_1(:, 1)';
    local(iel).y = centrecoord(iel).coord(1, 2) + sinew(:) * localXY_sii_1(:, 2)';

    % ======================================================
    % calculate pointdisp
    % start with sii^(-supelvalue) <= this one doesnt need to be in the loop, constant everywhere
    % should move them out a the loop later
    % sii^(-supelvalue) is a matrix A(i, j) => i is supelvalue, j is sii
    c_mat = repmat(supelc(:), [1 numel(sinew)]);
    si_mat = repmat(sinew(:)', [numel(supelc) 1]);
    lambda_mat = repmat(supelvalues(:), [1 numel(sinew)]);

    % c_sii_supelvalues_mat = repmat(supelc(:), [1 numel(sinew)]) .* (repmat(sinew(:)', [numel(supelc) 1]) .^ repmat(-supelvalues(:), [1 numel(sinew)]));
    c_sii_supelvalues_mat = c_mat .* (si_mat .^ -lambda_mat);
    % c_sii_supelvalues_minus1_mat = repmat(supelc(:), [1 numel(sinew)]) .* (repmat(sinew(:)', [numel(supelc) 1]) .^ repmat(-supelvalues(:), [1 numel(sinew)]) - 1);
    c_sii_supelvalues_minus1_mat = c_mat .* (si_mat .^ (-lambda_mat - 1));

    % ======================================================
    % below need to be in the loop
    % transpose matrix to get disp(i, j) => i is sii, j is nita
    pointdisp = real(N_mat * local(iel).phi * c_sii_supelvalues_mat)';

    % split the matrix to store in local
    local(iel).dispx = pointdisp(:, 1:2:end);
    local(iel).dispy = pointdisp(:, 2:2:end);

    % calculate stress without G*matmtrx, stress has size [3 * nita, si]
    stressXY = real(B1_mat * local(iel).phi * (c_sii_supelvalues_minus1_mat .* -lambda_mat) ...
                    + B2_mat * local(iel).phi * c_sii_supelvalues_minus1_mat);

    % reshape stress matrix to [3 newRow] to apply the G*matmtrx [3, 3]
    originalRow = size(stressXY, 1);
    originalCol = size(stressXY, 2);
    lengthOfStressXY = originalRow * originalCol / 3;
    stressXY = G * matmtrx * reshape(stressXY, [3 lengthOfStressXY]);
    % reshape the stress matrix back to normal to store the value
    stressXY = reshape(stressXY, [originalRow originalCol])';

    local(iel).strx = stressXY(:, 1:3:end);
    local(iel).stry = stressXY(:, 2:3:end);

%     for i=1:1:length(sinew)
%         for j=1:1:length(nitanew)
%             sii= sinew(i);
%             nitaa= nitanew(j);
%            % m= i; row
%            % n= j; col
%             [N,shapelq,dhdnitalq,Nder]= lobatto(nitaa,p(iel));%row-vectors
%             [modij, invmodij]=sbfemodijcirlipa(nnel,centrecoord(iel).coord,shapelq,dhdnitalq,xcoord,ycoord);   
%             [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
            
%             pointx=centrecoord(iel).coord(1,1);
%             pointy=centrecoord(iel).coord(1,2);
%             for ii=1:length(shapelq)
% %                 xcoord(ii)-centrecoord(iel).coord(ii,1)
%                 pointx= pointx + sii*shapelq(ii) *(xcoord(ii)-centrecoord(iel).coord(ii,1));
% %                 ycoord(ii)-centrecoord(iel).coord(ii,2)
%                 pointy= pointy + sii*shapelq(ii) *(ycoord(ii)-centrecoord(iel).coord(ii,2));
%             end
%             local(iel).x(i,j)=pointx; local(iel).y(i,j)=pointy;
            
%             pointdisp= zeros(2,1); pstr= zeros(3,1);
%             % IIfufpointdisp= zeros(2,1); IIcufpointdisp= zeros(2,1);
%             % IIfufpointstress= zeros(3,1); %N0pstr_1 + N1pstr1;
%             % IIcufpointstress= zeros(3,1); %N0pstr_1 + N1pstr1;


%             N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(2))];
%             si_star = 0;

%             for jj= 1:sdof
%                 pointdisp= real(pointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * local(iel).phi(:,jj));% 0^a= NaN
%                 pstr= real(pstr + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) *B1 +B2) * local(iel).phi(:,jj));% 0^a= NaN
%                 % IIfufpointdisp= real(IIfufpointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * coeffinestr(:,jj));
%                 % IIcufpointdisp= real(IIcufpointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * coefcoarsestr(:,jj));
%                 % IIfufpointstress= real(IIfufpointstress + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) * B1 + B2) * coeffinestr(:,jj));
%                 % IIcufpointstress= real(IIcufpointstress + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) * B1 + B2) * coefcoarsestr(:,jj));
%             end
                                
%             local(iel).dispx(i,j)=pointdisp(1); local(iel).dispy(i,j)=pointdisp(2);
%             local(iel).strx(i,j)= pstr(1); local(iel).stry(i,j)= pstr(2);
%             local(iel).strtxy(i, j) = pstr(3);
%             % local(iel).recovered_strx(i, j) = si_star(1); local(iel).recovered_stry(i, j) = si_star(2);
%            % local(iel).udisp(i,j).disp= pointdisp; 
%            % local(iel).pstress(i,j).pstr= pstr; 
          
%             % local(iel).IIfufpdx(i,j)= IIfufpointdisp(1); local(iel).IIcufpdx(i,j)= IIcufpointdisp(1);
%             % local(iel).IIfufpdy(i,j)= IIfufpointdisp(2); local(iel).IIcufpdy(i,j)= IIcufpointdisp(2);
%            % local(iel).IIfufpd(i,j).disp= IIfufpointdisp; local(iel).IIcufpd(i,j).disp= IIcufpointdisp;
            
%             % local(iel).IIfufpsx(i,j)= IIfufpointstress(1); local(iel).IIcufpsx(i,j)= IIcufpointstress(1);
%             % local(iel).IIfufpsy(i,j)= IIfufpointstress(2); local(iel).IIcufpsy(i,j)= IIcufpointstress(2);
%            % local(iel).IIfufps(i,j).stress= IIfufpointstress; local(iel).IIcufps(i,j).stress= IIcufpointstress;
%         end
%     end
end




