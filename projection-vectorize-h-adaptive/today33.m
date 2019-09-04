function [local]= sbfedstrecplot1(p,nnel,ndof,nnode,nel,nodes,gcoord,supelphi1,centrecoord,supelc,supelvalues,emodule,poisson)

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
sinew=[1:-0.1:0];%later should have bounded and unbounded as input, and if need should have different interval for dif ele to increase no of points for certain ele

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
     for i= 1:length(local(iel).nd)
        xcoord(i)= gcoord(local(iel).nd(i), 1);
        ycoord(i)= gcoord(local(iel).nd(i), 2);
    end
    [index]= sbfeeldof(local(iel).nd,2);
    for i= 1:edof(iel)                       % using index to extract phi(edof*sdof) from phi1
        for j= 1:sdof
        local(iel).phi(i,j)= supelphi1(index(i),j);
        end
    end
    coeffinestr = local(iel).phi;
    coefcoarsestr = [local(iel).phi(1:(edof(iel)-4),:);zeros(2,sdof);local(iel).phi((edof(iel)-1):edof(iel),:)];
    local(iel).coeffstr = coeffinestr;
    local(iel).coefcstr = coefcoarsestr;
    for i=1:1:length(sinew)
        for j=1:1:length(nitanew)
            sii= sinew(i);
            nitaa= nitanew(j);
           % m= i; row
           % n= j; col
            [N,shapelq,dhdnitalq,Nder]= lobatto(nitaa,p(iel));%row-vectors
            [modij, invmodij]=sbfemodijcirlipa(nnel,centrecoord(iel).coord,shapelq,dhdnitalq,xcoord,ycoord);   
            [b1, b2, B1, B2]=sbfederiv2(nnel,shapelq,dhdnitalq,invmodij);
            
            pointx=centrecoord(iel).coord(1,1);
            pointy=centrecoord(iel).coord(1,2);
            for ii=1:length(shapelq)
%                 xcoord(ii)-centrecoord(iel).coord(ii,1)
                pointx= pointx + sii*shapelq(ii) *(xcoord(ii)-centrecoord(iel).coord(ii,1));
%                 ycoord(ii)-centrecoord(iel).coord(ii,2)
                pointy= pointy + sii*shapelq(ii) *(ycoord(ii)-centrecoord(iel).coord(ii,2));
            end
            local(iel).x(i,j)=pointx; local(iel).y(i,j)=pointy;
            
            pointdisp= zeros(2,1); pstr= zeros(3,1);
            % IIfufpointdisp= zeros(2,1); IIcufpointdisp= zeros(2,1);
            % IIfufpointstress= zeros(3,1); %N0pstr_1 + N1pstr1;
            % IIcufpointstress= zeros(3,1); %N0pstr_1 + N1pstr1;


            N_for_sigma = [diag(ones(1, 3)*shapelq(1)), diag(ones(1, 3)*shapelq(2))];
            si_star = 0;

            for jj= 1:sdof
                pointdisp= real(pointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * local(iel).phi(:,jj));% 0^a= NaN
                pstr= real(pstr + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) *B1 +B2) * local(iel).phi(:,jj));% 0^a= NaN
                % IIfufpointdisp= real(IIfufpointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * coeffinestr(:,jj));
                % IIcufpointdisp= real(IIcufpointdisp + N * supelc(jj) * sii^(-supelvalues(jj)) * coefcoarsestr(:,jj));
                % IIfufpointstress= real(IIfufpointstress + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) * B1 + B2) * coeffinestr(:,jj));
                % IIcufpointstress= real(IIcufpointstress + G* matmtrx * supelc(jj) * sii^(-supelvalues(jj) -1) * (-supelvalues(jj) * B1 + B2) * coefcoarsestr(:,jj));
            end
                                
            local(iel).dispx(i,j)=pointdisp(1); local(iel).dispy(i,j)=pointdisp(2);
            local(iel).strx(i,j)= pstr(1); local(iel).stry(i,j)= pstr(2);
            local(iel).strtxy(i, j) = pstr(3);
            % local(iel).recovered_strx(i, j) = si_star(1); local(iel).recovered_stry(i, j) = si_star(2);
           % local(iel).udisp(i,j).disp= pointdisp; 
           % local(iel).pstress(i,j).pstr= pstr; 
          
            % local(iel).IIfufpdx(i,j)= IIfufpointdisp(1); local(iel).IIcufpdx(i,j)= IIcufpointdisp(1);
            % local(iel).IIfufpdy(i,j)= IIfufpointdisp(2); local(iel).IIcufpdy(i,j)= IIcufpointdisp(2);
           % local(iel).IIfufpd(i,j).disp= IIfufpointdisp; local(iel).IIcufpd(i,j).disp= IIcufpointdisp;
            
            % local(iel).IIfufpsx(i,j)= IIfufpointstress(1); local(iel).IIcufpsx(i,j)= IIcufpointstress(1);
            % local(iel).IIfufpsy(i,j)= IIfufpointstress(2); local(iel).IIcufpsy(i,j)= IIcufpointstress(2);
           % local(iel).IIfufps(i,j).stress= IIfufpointstress; local(iel).IIcufps(i,j).stress= IIcufpointstress;
        end
    end
end




