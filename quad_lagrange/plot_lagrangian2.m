collectN = [];

for i = -1:0.01:1
%     [~, temp] = lagrangian2(i, 2);
    temp = quadLagrangHierarchial(i, 2, 0);
    collectN = [collectN; temp.NN];    
end

temp = [];
for i = -1:0.01:0
    t = -1 + 2*(i--1);
    t = 1 - t^2;
    temp = [temp; t];
end

temp = [temp; temp(2:end)];
% plot(collectN)

xm1 = [1000; 2000];
x1 = [3000; 4000];
x0 = mean([xm1, x1], 2);
xm5 = mean([xm1, x0], 2);
x5 = mean([x0, x1], 2);

u = [xm1, x0, x1, xm5, x5];
result = [];
for i = -1:0.1:1
    temp = quadLagrangHierarchial(i, 1, 0);
    ui = u * temp.NN';
    result = [result, ui];
end

plot(result(1, :), result(2, :))