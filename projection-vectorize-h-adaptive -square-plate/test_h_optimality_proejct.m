ngl = 64;
[p64,w64] = sbfeglqd1(ngl);

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

k1 = dNj * (dNhalf1 .* w64(:));
k2 = dNj * (dNhalf2 .* w64(:));

%%
pf1 = [1 1 1];
pf2 = [2 1 2];
pc = [0 0 0];
pc(1) = pf1(1);
pc(3) = pf2(3);

% add random middle node from 0 -> 1
middleNode = -5:0.01:5;
pc = repmat(pc(:), 1, numel(middleNode)); 
pc(2, :) = middleNode;  % add middle node

dN_fine = [];
dN_c_half1 = [];
dN_c_half2 = [];

N_fine = [];
N_h1 = [];
N_h2 = [];
% each fine element split into 64 points [-1, 1] , therefore
% each half coarse element also split into 64 points [-1, 0] and [0, 1]
nglob1 = [];
nglob2 = [];
p64 = -1:0.01:1;
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
    
    N_fine = [N_fine; Nfine];
    N_h1 = [N_h1; Nc1];
    N_h2 = [N_h2; Nc2];
    
    nglob1 = [nglob1, nita_glob1];
    nglob2 = [nglob2, nita_glob2];    
end

f1 = dN_fine * pf1(:);
f2 = dN_fine * pf2(:);
ph1 = dN_c_half1 * pc;
ph2 = dN_c_half2 * pc;

% f1 = repmat(f1(:), 1, numel(midd

int_delta = sum((repmat(f1(:), 1, numel(middleNode)) - ph1).^2, 1) ...
            + sum((repmat(f2(:), 1, numel(middleNode)) - ph2).^2, 1);
plot(middleNode, int_delta')
pf1 * k1 + pf2 * k2
%%
plot(nglob1, f1); hold on
plot(nglob2, f2);
% plot(nglob1, ph1);
% plot(nglob2, ph2);