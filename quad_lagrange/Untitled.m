p3 = [1000 1500 2000; 0 250 500];

nitanew = [-1:0.1:1];

pc = [];
for i = 1:numel(nitanew)
    [N, NN] = lagrangian2(nitanew(i), 2);
    pc = [pc, p3 * NN'];
end