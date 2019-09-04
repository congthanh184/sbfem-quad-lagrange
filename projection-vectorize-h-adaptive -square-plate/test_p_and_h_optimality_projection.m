ngl = 64;
[p64, w64] = sbfeglqd1(ngl);
%% p-adaptive
% define u_fine
pf = [0.2 0.1 -0.44 -1.9];
pc = [0.2 0 0.8]; % project two end nodes, should be called u_project

% add random middle node from 0 -> 1
middleNode = -1:0.01:1;
pc = repmat(pc(:), 1, numel(middleNode)); 
pc(2, :) = middleNode;  % add middle node

N_p2 = [];
N_p3 = [];

dN_p2 = [];
dN_p3 = [];


for ii = 1:numel(p64)
    [~, N2, dN2] = lobatto( p64(ii), 2);
    [~, N3, dN3] = lobatto( p64(ii), 3);
    N_p2 = [N_p2; N2];
    N_p3 = [N_p3; N3];
    dN_p2 = [dN_p2; dN2];
    dN_p3 = [dN_p3; dN3];
end

u_coarse = N_p2 * pc;
u_fine = N_p3 * pf(:);

du_c = dN_p2 * pc;
du_f = dN_p3 * pf(:);
% calculate norm(u_coarse - u_fine)
delta_dn = du_c - repmat(du_f, 1, numel(middleNode));
delta_dn_sqr = delta_dn .^2;

% integral_delta = w64' * delta_dn_sqr;
integral_delta = sum(delta_dn_sqr, 1);

plot(middleNode, integral_delta, 'LineWidth', 3);

%% h-adaptive

% define u_fine 1 and u_fine 2
pf1 = [1 -1 1];
pf2 = [2 0.5 2];
pc = [0 0 0]; % should be called u_project
pc(1) = pf1(1);
pc(3) = pf2(3);

% add random middle node from -1 -> 1
middleNode = -1:0.01:1;
pc = repmat(pc(:), 1, numel(middleNode)); 
pc(2, :) = middleNode;  % add middle node

dN_fine = [];
dN_c_half1 = [];
dN_c_half2 = [];

% each fine element split into 64 points [-1, 1] , therefore
% each half coarse element also split into 64 points [-1, 0] and [0, 1]
nglob1 = [];
nglob2 = [];
for ii = 1:numel(p64)
    nita = p64(ii);
    [~, Nfine, dNfine] = lobatto( p64(ii), 2);
    nita_glob1 = 0.5 * (nita - 1);
    nita_glob2 = 0.5 * (nita + 1);
    [~, Nc1, dNc1] = lobatto( nita_glob1, 2);
    [~, Nc2, dNc2] = lobatto( nita_glob2, 2);
    
    dN_fine = [dN_fine; dNfine];
    dN_c_half1 = [dN_c_half1; dNc1];
    dN_c_half2 = [dN_c_half2; dNc2];
    nglob1 = [nglob1, nita_glob1];
    nglob2 = [nglob2, nita_glob2];    
end

% calculate derivative u_f1 and u_f2
du_f1 = dN_fine * pf1(:);
du_f2 = dN_fine * pf2(:);
du_c_h1 = dN_c_half1 * pc;
du_c_h2 = dN_c_half2 * pc;

% calculate delta derivative of uf - uc
delta_h1 = du_c_h1 - repmat(du_f1, 1, numel(middleNode));
delta_h2 = du_c_h2 - repmat(du_f2, 1, numel(middleNode));

delta_h1 = delta_h1 .^ 2;
delta_h2 = delta_h2 .^ 2;

integral_delta = (w64' * delta_h1) + (w64' * delta_h2);
plot(middleNode, integral_delta, 'LineWidth', 3);
ngl = 64;
[p64,w64] = sbfeglqd1(ngl);

% calculate expected middle node value
dNj = [];
dNhalf1 = [];
dNhalf2 = [];

for i = 1:ngl
    [N, nn, dnn] = lobatto( p64(i), 2);
    dNj = [dNj, dnn(:)];
    [N, nn, dnn] = lobatto( 0.5 * p64(i) - 0.5, 2);
    dNhalf1 = [dNhalf1; dnn(2)];
    [N, nn, dnn] = lobatto( 0.5 * p64(i) + 0.5, 2);
    dNhalf2 = [dNhalf2; dnn(2)];
end

k1 = dNj * (dNhalf1 .* w64(:)) / 2;
k2 = dNj * (dNhalf2 .* w64(:)) / 2;
fprintf(1, 'Calculation value %f \n', (pf1 * k1 + pf2 * k2));
% K = [k1.', k2.'];
