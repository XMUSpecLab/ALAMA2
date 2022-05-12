%% Block matching SVD
function Y = BM_SVD(I, ds, Ds, step, sn)
%I:含噪声光谱 1D or 2D
%ds:邻域窗口半径  
%Ds:搜索窗口半径  
%step: 搜素步长
%sn: Block个数  
%Y:去噪光谱 1D or 2D
%% Step1: 基础估计
[m,n]=size(I);
if n==1
    I = [I,I];
end
% 归一化
minima = min(I(:));
maxima = max(I(:));
I=(I - minima)/(maxima-minima);
Padded=[I(1:ds,:);I;I(end-ds+1:end,:)]; 
numerator = zeros(size(Padded));
denomerator = zeros(size(Padded));

for i=1:step:m
    for j=1:n
        i1=i+ds;
        % 邻域窗口1
        W1 = Padded(i1-ds:i1+ds,j);
        % 搜索窗口
        rmin = max(i1-Ds,1);
        rmax = min(i1+Ds,m+2*ds);
        search_window = Padded(rmin:rmax,:);
        % 对搜索窗口进行patch划分，并返回所有patch
        G = grouping(search_window, ds, 1);
        % 计算Patch 之间的距离,并按距离由小到大排列所有patch
        [SG, ~] = sort_patch(G, W1, sn);
        % 奇异值分解
        [u,s,v] = svd(SG,'econ'); 
        s=s(s~=0); v = v';
        % 去噪
        recon = s(1)*u(:,1)*v(1,:);
        NW = recon(:,1);
        % 权重
        wp = s(1);
        % 加权平均
        numerator(i1-ds:i1+ds,j) = ...
            numerator(i1-ds:i1+ds,j)  + wp*NW;
        denomerator(i1-ds:i1+ds,j) = ...
            denomerator(i1-ds:i1+ds,j) + wp;
    end
end
%基础估计
basic_result = numerator./denomerator;
basic = basic_result(ds+1:end-ds,:);
Y = basic*(maxima-minima)+minima;
end