% process ET_save from h-adaptive p1 step 6
saveData = load('h_projection_fixed_p1_D.mat');

ETTable = saveData.ET_save{1, 6};

% delete noisy rows due to numerical sensitive
ETTable([13:17], :) = 0;

bot = sum( real(-ETTable(:, 1)));
top = sum( real(-(ETTable(:, 1) - ETTable(:, 2) - ETTable(:, 3) + ETTable(:, 4))));

100 * (top^0.5) / (bot^0.5)

deltaerror = sqrt(real(ETTable(:, 8)));
error_array = sqrt(ETTable(:, 5));
% error_array([13:17]) = [];
rhs65 = (1/3) * max(deltaerror);
aaa = error_array;
aaa([13:17]) = [];
rhs66 = 10 * mean(aaa);

refineIdx = (deltaerror > rhs65) | (error_array > rhs66);

sum(refineIdx)

%%
saveData = load('h_projection_fixed_p2_D.mat');

ETTable = saveData.ET_save{1, 5};

% delete noisy rows due to numerical sensitive
ETTable([4:7], :) = 0;

bot = sum( real(-ETTable(:, 1)));
top = sum( real(-(ETTable(:, 1) - ETTable(:, 2) - ETTable(:, 3) + ETTable(:, 4))));

100 * (top^0.5) / (bot^0.5)

deltaerror = sqrt(real(ETTable(:, 8)));
error_array = sqrt(ETTable(:, 5));
% error_array([13:17]) = [];
rhs65 = (1/3) * max(deltaerror);
aaa = error_array;
aaa([4:7]) = [];
rhs66 = 10 * mean(aaa);

refineIdx = (deltaerror > rhs65) | (error_array > rhs66);

sum(refineIdx)