%% Grouping 分组并返回所有小片段
function G = grouping(search_window, window_size, step)
%% 参数设置
% search_window：整个搜索区间 
% window_size: 搜索框的尺寸半径
%% 执行列分组
[m,n] = size(search_window);
% 总共有多少个片段
% p 存储所有片段
G = [];
for i = 1:n
    for j=window_size+1:step:m-window_size
        tmp = search_window(j-window_size:j+window_size,i);
        G = [tmp,G];
    end
end
end