%%
clear; close all; clc;
load('syt062.mat');
% denoising
R = ALRMA2(cube);

%% demonstration
idx = randperm(size(R,2), 9);
figure;

for i = 1:9
    subplot(3,3,i);
    plot(x, cube(:,idx(i)),'color','b','LineWidth', 1);
    hold on;
    plot(x, R(:,idx(i)),'color','r','LineWidth', 1);
    axis tight;
    title(num2str(i));
    legend('Noisy','Denoised');
end
