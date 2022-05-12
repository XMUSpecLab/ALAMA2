%% 计算每个patch和目标patch之间的距离,并且按距离从小到大重新排序
% 返回给定个patch,及其距离
function [sorted_g, I] = sort_patch(g, target_patch, sn)
%% 参数
% g: 光谱片段组，2D
% target_patch : 目标光谱片段
% theta: 硬阈值滤波阈值
% sn: 选择多少个patch
% sorted_g: 排序好的块
% I: 选择的序号
%% 计算所有片段和目标片段的距离
[~, num] = size(g);
dist = zeros(1, num);
% 对目标片段进行硬阈值滤波
for i = 1:num
    dist(i) = norm(target_patch-g(:,i));
end
[~, I] = sort(dist);
tmp = g(:,I);
if sn>num
    sn = num;
end
sorted_g = tmp(:,1:sn);
I = I(1:sn);
end