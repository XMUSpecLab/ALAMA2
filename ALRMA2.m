%% An improved version of ALRMA
function [R, Sel, num] = ALRMA2(cube,type, num, thr)
%% parameters 
% data: 2D matrix of the signal 
% type: kmeans subclasses
% num: determine how many svd components are in considerate, 30 as default
% thr: determine which svd components is selected, 1e-2 as default
% Reconstructed: the denoised signal using svd denoising
% Index_selected: the selected SVD component index
% num: how many top svd components are examined in the current process
%% config the default parameters
switch nargin
    case 1
        type = 2;
        num = 30;
        thr = 1e-2;
    case 2
        num = 30;
        thr = 1e-2;
    case 4
        disp(' ');
    otherwise
        error("Need 4 parameters: data,columns,imaging region, window size, count, thr1, thr2, but only %d is given!\n", nargin);
end

%% K-means clustering
[class, ~, ~, D] = kmeans(cube', type);
sub_group = {};
sub_group_idx = {}; 
% 输出每个子类距离质心最近的sn条光谱，以及它们在原始数据中的索引
for i=1:type
    group_idx = find(class==i);
    d = D(group_idx, i);
    % 根据d,在每个子类中选取sn条光谱
    sn = 5;
    if length(group_idx)<5
         error("Error: only %d element in a subclass!\n", length(group_idx));
    end
    [~, I] = sort(d);
    sub_group_idx{i} = group_idx(I(1:sn));
    sub_group{i}= cube(:,I(1:sn));
end

ref_spec = []; ref_idx = [];
for i=1:type
    ref_idx = [sub_group_idx{i};ref_idx];
    out = BM_SVD(sub_group{i}, 500, 600, 50, 25);
    ref_spec = [out, ref_spec];
end

%% SVD and stastical analysis
% determine how many top SVD components are in consideration
[U,S,V] = svd(cube,'econ');V = V';s = S(S~=0);
snrs = []; recon = 0; 
for i = 1:num
    tmp = s(i)*U(:,i)*V(i,:);
    recon = recon + tmp;
    for k = 1:length(ref_idx)
        newspec = recon(1:end-50,ref_idx(k));
        residual = ref_spec(1:end-50,k) - newspec;
        snrs(k,i) = snr(ref_spec(1:end-50,k), residual);
    end
end
%%
% determine which SVD componets is selected to reconstruct the signal
% according to the mean SNRS matrix
pad = zeros(length(ref_idx),1);
snrs = [pad, snrs]; dfsnr = [];
for i = 1:length(ref_idx)
    dfsnr(i,:) = diff(snrs(i,:));
end
mdfsnr = mean(dfsnr); [Sel,~]=find(mdfsnr'>=thr);
if nargout ~= 1
    figure;
    subplot 211; plot(mdfsnr,'b*-'); xlabel('SVs'); ylabel('SNR Contribution/dB'); title(['Mean curves of ', num2str(type), ' regions'])
    subplot 212; plot(dfsnr'); xlabel('SVs'); ylabel('SNR Contribution/dB'); 
end
ss = zeros(size(S)); 
for j = 1:length(Sel)
    ss(Sel(j),Sel(j))= s(Sel(j));
end
R = U*ss*V;
end
