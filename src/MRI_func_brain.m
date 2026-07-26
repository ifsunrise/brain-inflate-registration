function [varargout] = MRI_func_brain(sMRI,func,varargin)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

sMRI = MRI_read_Data(sMRI,'verbose',0);
switch func
    case 'seg_kmeans'
        output = mri_seg_kmeans(sMRI,varargin{:});
    case 'rand_crop'
        output = mri_randcrop(sMRI,varargin{:});
    case 'bw_seedgrowth'
        output = mri_bw_seedgrowth(sMRI,varargin{:});
    case 'directional_smooth'
        output = mri_directional_smooth(sMRI,varargin{:});
    case 'denoise'
        output = mri_denoise(sMRI);
    case 'pvcnorm'
        output = mri_pvcnorm(sMRI,varargin{:});
end

if ~iscell(output)
    output = {output};
end
varargout = output;
end

function sImg = mri_seg_kmeans(sMRI,varargin)
options_default = structure('N',3,'mask',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

temp = single(sMRI.data);

if ~isempty(mask)
    [L,L_cen] = imsegkmeans3(temp(mask>0),3);
    L_cen = sort(L_cen);
else
    [L,L_cen] = imsegkmeans3(temp,N);
    [L_cen,idx] = sort(L_cen);
    [~,idx] = sort(idx);
    L = reshape(idx(L),size(L));
end

for i = 1:3
    L_frac{i} = power(exp(1),-2*abs((temp-L_cen(i))));
end

L_sum = L_frac{1}+L_frac{2}+L_frac{3};
for i = 1:3
    temp = L_frac{i}./L_sum;
    csf_mean = mean(temp(mask==0),"omitnan");
    temp = (temp - csf_mean)./(1-csf_mean); temp(temp<0) = 0;
    L_frac{i} = temp;
end
end

function  [out] = mri_randcrop(img,varargin)
options_default = structure('mask',[],'sz',50);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if isstruct(img)
    img = img.data;
end

if length(sz) == 1
    sz = [sz,sz];
end

mask(1:round(sz(1)/2),:) = 0; mask(:,1:round(sz(2)/2)) = 0;
mask(end-round(sz(1)/2):end,:) = 0; mask(:,end-round(sz(2)/2):end) = 0;

idx = find(mask>0);
crop_idx = idx(randi(length(idx),1));

[y,x] = ind2sub(size(img),crop_idx);
crop_loc = [y,x] - round(sz/2)-1;

img_crop = img(crop_loc(1)+[1:sz(1)],crop_loc(2)+[1:sz(2)]);

out = {img_crop,crop_loc,crop_idx};
end

function bw_grow = mri_bw_seedgrowth(bw,varargin)
% 从一个种子点寻找邻域（如沿着灰白质边界走
options_default = structure('seed',[],'N',5);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

temp = single(bw);
temp(seed(1),seed(2)) = 100; % 设为一个大值
bw_grow = 0*single(bw);

for i = 1:N
    temp = imfilter(temp,ones(3),'same').*bw;
    temp(temp>20) = 100; temp(temp>0 & temp<20) = 1;
    bw_grow(temp>1) = bw_grow(temp>1)+1;
end
end

function img_filt = mri_directional_smooth(img,Gx0,Gy0,sm_mask,N)
img_filt = img;
for i = 1:N
    img_smx = imfilter(img_filt,[1/2 -1 1/2],'replicate');
    img_smy = imfilter(img_filt,[1/2;-1;1/2],'replicate');
    temp = img_smx.*Gx0+img_smy.*Gy0;
    img_filt = img_filt + temp.*sm_mask;
end
end

function output = mri_denoise(sMRI)
% from cat12: cat_vol_sanlm.m
% nonlocal means for denoising
job.intlim = [1,1]-0.0001;
job.rician = 0;

src = single(sMRI.data);

% histogram limit
[src,srcth] = cat_stat_histth(src,job.intlim);

% use intensity normalisation because cat_sanlm did not filter values below ~0.01
th  = max( cat_stat_nanmean( src(src(:)>cat_stat_nanmean(src(src(:)>0))) ) , ...
    abs(cat_stat_nanmean( src(src(:)<abs(cat_stat_nanmean(src(src(:)<0)))))) );

src = (src / th) * 100;
src = (src - srcth(1)); % avoid negative values!
cat_sanlm(src,3,1,job.rician);
src = src + srcth(1);  % restore original intensity range
src = (src / 100) * th;

sMRI_filt = sMRI;
sMRI_filt.data = src;

output = {sMRI_filt};
end

function output = mri_pvcnorm(sMRI,varargin)
options_default = structure('label',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

% from % [Ym] = cat_main_gintnorm1639(single(T1_raw.data),Tth);
Ysrc = single(sMRI.data);
if isequal(sMRI.size,size(label))
    T3th = 0; T3thx = 0;
    for i = 1:3
        T3th(i+1) = median(Ysrc(abs(label - i) < 0.3));
        T3thx(i+1) = i;
    end
else
    T3th = label.T3th;
    T3thx = label.T3thx;
end

Ym = 0*Ysrc;
for i=2:4
    M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
end
M  = Ysrc>=T3th(end);
Ym(M(:)) = T3thx(i) + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));
Ym = Ym / 3;

sMRI_norm = sMRI;
sMRI_norm.data = Ym;

output = {sMRI_norm,T3th};
end