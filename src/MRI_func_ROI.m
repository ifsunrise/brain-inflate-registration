
function [varargout] = MRI_func_ROI(sMRI,func,varargin)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if ~isstruct(sMRI)
    sMRI = MRI_read_Data(sMRI);
end

switch func
    case 'slice'
        output = slicing_mri(sMRI,varargin{:});
    case 'crop'
        output = crop_mri(sMRI,varargin{:});
    case 'excrop'
        output = excrop_mri(sMRI,varargin{:});
    case 'decrop'
        output = decrop_mri(sMRI,varargin{:});
    case 'repos'
        output = repos_mri(sMRI,varargin{:});
    case 'update_Tmat'
        output = update_Tmat(sMRI,varargin{:});
    case 'imresize3'
        output = imresize3_mri(sMRI,varargin{:});
    case 'rescale'
        output = rescale_mri(sMRI,varargin{:});
    case 'sample'
        output = sample_mri(sMRI,varargin{:});
    case 'cat'
        output = cat_mri(sMRI,varargin{:});
    case 'permute'
        output = permute_mri(sMRI,varargin{:});
    case 'make_grid'
        output = make_grid(sMRI,varargin{:});
    case 'make_bubbles'
        output = make_bubbles(sMRI,varargin{:});
    case 'make_subcort'
        output = make_subcort(sMRI,varargin{:});
    case 'dist_map'
        output = dist_map(sMRI,varargin{:});
    case 'smooth_field'
        output = smooth_field(sMRI,varargin{:});
    case 'fuse_rgb'
        output = fuse_rgb(sMRI,varargin{:});
end

if ~iscell(output)
    output = {output};
end
varargout = output;
end

function sImg = slicing_mri(sMRI,varargin)
options_default = struct('dim','z','index',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

switch dim
    case 'x'
    case 'y'
    case 'z'
        sImg.data = sMRI.data(:,:,index);
        sImg.pixsize = sMRI.voxsize([1,2]);
end
end

function sROI = crop_mri(sMRI, varargin)
% 裁剪 MRI 体积，同时更新 T_mat 以保证世界坐标对齐
%
% 输入：
%   sMRI   - 符合自定义结构（至少包含 data, T_mat, voxsize）的 MRI 结构体
%   roi    - 裁剪范围，可以是：
%            * 2x3 数值矩阵 [x1, y1, z1; x2, y2, z2]（闭区间）
%            * 1x3 cell 数组 {x_idx, y_idx, z_idx}，每个元素为索引向量
% 输出：
%   sROI   - 裁剪后的 MRI 结构体，其 T_mat 已更新，世界坐标与原体积一致

options_default = struct('roi', []);
[options, eval_str] = resolve_input(options_default, varargin);
eval(eval_str);

% 解析裁剪范围
if isnumeric(roi)
    % roi 格式：[x1 y1 z1; x2 y2 z2]
    y_ind = roi(1,2) : roi(2,2);   % data 第一维为 y
    x_ind = roi(1,1) : roi(2,1);   % data 第二维为 x
    z_ind = roi(1,3) : roi(2,3);   % data 第三维为 z
else
    % roi 格式：{x_idx, y_idx, z_idx}
    x_ind = roi{1};
    y_ind = roi{2};
    z_ind = roi{3};
end

% 裁剪起始索引（1‑based）
y1 = y_ind(1);
x1 = x_ind(1);
z1 = z_ind(1);

% 复制结构体并替换数据
sROI = sMRI;
sROI.data = sMRI.data(y_ind, x_ind, z_ind, :);

% 更新仿射变换矩阵，保持世界坐标不变
% 原约定：行向量 [y, x, z, 1] 右乘 T_mat 得到世界坐标
if isfield(sROI, 'T_mat') && iscell(sROI.T_mat) && size(sROI.T_mat,2) >= 2
    T = sROI.T_mat{1,2};   % 4x4 矩阵
    % 更新平移向量（第四行前三列）
    if ~isempty(T)
        T(4,1) = T(4,1) + (y1 - 1) * T(1,1);
        T(4,2) = T(4,2) + (x1 - 1) * T(2,2);
        T(4,3) = T(4,3) + (z1 - 1) * T(3,3);
        sROI.T_mat{1,2} = T;
    end
end

% 保留原始体素尺寸
sROI.voxsize = sMRI.voxsize;

% 调用 MRI_read_Data 对 struct 类型重新计算 size、shape 等辅助字段
sROI = MRI_read_Data(sROI);
end

function sROI = excrop_mri(sMRI, varargin)
% 扩展裁剪（Extendable Crop）：允许裁剪矩形的起止索引落在原始体积之外。
% 超出范围的体素用 val 填充，同时更新 T_mat 保持世界坐标不变。
%
% 输入：
%   sMRI - mri 结构体（需含 data, T_mat, voxsize）
%   roi  - 2x3 矩阵，格式 [x1, y1, z1; x2, y2, z2]（闭区间）
%          x1,y1,z1 可以是任意整数（包括 0、负值或大于原始尺寸）
%   val  - 填充值，默认 0
%
% 输出：
%   sROI - 裁剪/扩展后的结构体，data 尺寸为 [ny, nx, nz, nt]
%          其 T_mat 已调整，保证体素的世界坐标与原体积对应位置完全一致

% 解析输入
options_default = struct('roi', [], 'val', 0);
[options, eval_str] = resolve_input(options_default, varargin);
eval(eval_str);

if isempty(roi) || ~isnumeric(roi) || ~isequal(size(roi), [2,3])
    error('roi 必须为 2x3 数值矩阵，格式 [x1,y1,z1; x2,y2,z2]。');
end

% 起始索引（y,x,z 顺序──data 第一维是 y，第二维是 x）
y1 = roi(1,2); x1 = roi(1,1); z1 = roi(1,3);
y2 = roi(2,2); x2 = roi(2,1); z2 = roi(2,3);

% 原始尺寸（y, x, z, t）
orig_sz = size(sMRI.data);
if length(orig_sz) < 4
    orig_sz(4) = 1;
end
ny_orig = orig_sz(1);
nx_orig = orig_sz(2);
nz_orig = orig_sz(3);

% 新体积尺寸
ny_new = y2 - y1 + 1;
nx_new = x2 - x1 + 1;
nz_new = z2 - z1 + 1;
if ny_new <= 0 || nx_new <= 0 || nz_new <= 0
    error('新体积尺寸必须为正。');
end

% 创建新数据，全部预填 val
if isfloat(sMRI.data)
    new_data = val * ones(ny_new, nx_new, nz_new, orig_sz(4), 'like', sMRI.data);
else
    new_data = cast(val, 'like', sMRI.data) * ...
               ones(ny_new, nx_new, nz_new, orig_sz(4), 'like', sMRI.data);
end

% 确定原数据落在新体积中的有效索引，并复制数据
new_y_vec = y1 : y2;
y_valid = new_y_vec >= 1 & new_y_vec <= ny_orig;
new_y_idx = find(y_valid);
orig_y_idx = new_y_vec(y_valid);

new_x_vec = x1 : x2;
x_valid = new_x_vec >= 1 & new_x_vec <= nx_orig;
new_x_idx = find(x_valid);
orig_x_idx = new_x_vec(x_valid);

new_z_vec = z1 : z2;
z_valid = new_z_vec >= 1 & new_z_vec <= nz_orig;
new_z_idx = find(z_valid);
orig_z_idx = new_z_vec(z_valid);

new_data(new_y_idx, new_x_idx, new_z_idx, :) = ...
    sMRI.data(orig_y_idx, orig_x_idx, orig_z_idx, :);

% 构建输出结构体
sROI = sMRI;
sROI.data = new_data;
sROI.voxsize = sMRI.voxsize;

% 更新仿射矩阵：使新体积的 (1,1,1) 体素对应原体积中 (x1, y1, z1) 的世界坐标
% 注：T_mat{1,2} 为行向量前乘形式，平移在第四行，对角线依次对应 x,y,z 方向缩放
if isfield(sROI, 'T_mat') && iscell(sROI.T_mat) && size(sROI.T_mat,2) >= 2
    T = sROI.T_mat{1,2};
    if ~isempty(T)
        T(4,1) = T(4,1) + (x1 - 1) * T(1,1);   % x 方向
        T(4,2) = T(4,2) + (y1 - 1) * T(2,2);   % y 方向
        T(4,3) = T(4,3) + (z1 - 1) * T(3,3);   % z 方向
        sROI.T_mat{1,2} = T;
    end
end

% 调用 MRI_read_Data 刷新 size, shape 等辅助字段
sROI = MRI_read_Data(sROI);
end

function sVOL_reo = repos_mri(sVOL,sRef,varargin)
options_default = struct('verbose',1,'interp_type','linear','fill_val',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

temp = sVOL.data; sVOL_reo = sRef;
sVOL_reo.data = cast(fill_val + zeros([sRef.size,size(sVOL.data,4)]),'like',temp);

phyShape_ref = sRef.shape.*sRef.voxshape;
phyShape_vol = sVOL.shape.*sVOL.voxshape;

T = sRef.T_mat{1,2};
pos_ref = T(4,1:3).*sign(T([1,6,11]));
phyOri_ref = phyShape_ref+pos_ref;

T = sVOL.T_mat{1,2};
pos_vol = T(4,1:3).*sign(T([1,6,11]));
phyOri_vol = phyShape_vol+pos_vol;

dOri_ref = (phyShape_ref-phyOri_ref) - (phyShape_vol-phyOri_vol) + 1;
dShape_ref = round(phyShape_vol./sRef.voxshape);
lim_ref = round(dOri_ref./sRef.voxshape)+[0 0 0; dShape_ref-1];
over_ref = (lim_ref - [1 1 1; sRef.shape]).*[-1;1];
over_ref(over_ref<0) = 0;
lim_ref = round(lim_ref - over_ref.*[-1;1]);

lim_vol = [1 1 1; sVOL.shape]  - over_ref.*[-1;1].*sRef.voxshape./sVOL.voxshape;

temp = lim_ref(2,:)-lim_ref(1,:)+1;
sVOL = MRI_func_ROI(sVOL,'crop',round(lim_vol));
sVOL = MRI_func_ROI(sVOL,'imresize3',temp([2,1,3]) ,'interp_type',interp_type);
sVOL_reo.data(lim_ref(1,2):lim_ref(2,2),lim_ref(1,1):lim_ref(2,1),lim_ref(1,3):lim_ref(2,3),:) = sVOL.data;

if verbose
    disp(strcat(num2str(round(lim_vol)),' => ',num2str(round(lim_ref))))
end
end

% deepseek, 有问题
function sVOL_reo = repos_mri_new(sVOL, sRef, varargin)
% REPOS_MRI  将 sVOL 重采样到 sRef 的网格空间，保持世界坐标对齐
%
%   输出与 sRef 具有相同的尺寸、voxsize 和 T_mat，数据来自对 sVOL 的插值
%
% 输入：
%   sVOL, sRef : MRI 结构体，需包含 data 和 T_mat{1,2}
% 可选参数：
%   'verbose'    : 是否显示信息 (0/1, 默认 1)
%   'interp_type': 插值类型 ('linear','nearest','cubic', 默认 'linear')
%   'fill_val'   : 边界外填充值 (默认 0)
% 输出：
%   sVOL_reo     : 重采样后的 MRI 结构体，与 sRef 同空间属性

options_default = struct('verbose',1, 'interp_type','linear', 'fill_val',0);
[options, eval_str] = resolve_input(options_default, varargin);
eval(eval_str);

% 检查 T_mat 存在性
if ~isfield(sVOL,'T_mat') || isempty(sVOL.T_mat) || ...
   ~isfield(sRef,'T_mat') || isempty(sRef.T_mat)
    error('repos_mri 要求输入的 MRI 结构体包含有效的 T_mat{1,2} (4x4 仿射矩阵)');
end

T_vol = sVOL.T_mat{1,2};   % 将体素索引 [x;y;z;1] 映射到世界坐标
T_ref = sRef.T_mat{1,2};

% 参考网格的尺寸 (y, x, z)
ny = sRef.size(1);
nx = sRef.size(2);
nz = sRef.size(3);

% 1. 计算 sRef 网格覆盖的世界坐标范围（体素角点）
%    取 8 个角点 (1 和 size+1)，经 T_ref 映射到世界坐标
[cx, cy, cz] = ndgrid([1, nx+1], [1, ny+1], [1, nz+1]);
corners = [cx(:), cy(:), cz(:), ones(numel(cx),1)]';
worldCorners = T_ref * corners;   % 4 x 8
xWorldLim = [min(worldCorners(1,:)), max(worldCorners(1,:))];
yWorldLim = [min(worldCorners(2,:)), max(worldCorners(2,:))];
zWorldLim = [min(worldCorners(3,:)), max(worldCorners(3,:))];

% 2. 构建 sRef 的输出空间参考
Rout = imref3d([ny, nx, nz], xWorldLim, yWorldLim, zWorldLim);

% 3. 计算 sVOL 体素索引 → sRef 体素索引 的仿射变换矩阵
%    由 w = T_vol * X_vol = T_ref * X_ref  →  X_ref = T_ref \ T_vol * X_vol
M = T_ref \ T_vol;               % 4x4
tform = affine3d(M);             % imwarp 所需的变换对象

% 4. 执行重采样（自动处理多通道，保持数据类型）
sVOL_reo = sRef;
sVOL_reo.data = imwarp(sVOL.data, tform, ...
    'OutputView', Rout, ...
    'FillValues',  fill_val, ...
    'Interp',     interp_type);

if verbose
    fprintf('repos: %s → %s\n', mat2str(sVOL.size), mat2str(sRef.size));
end
end

function sMRI_new = update_Tmat(sMRI, sRef, T_new)
% 当 sRef 的 T_mat 变为 T_new 时，同步更新 sMRI2 的 T_mat，
% 以保持两者在世界坐标中的相对位置不变。
%
% 输入：
%   sMRI1, sMRI2 : MRI 结构体，需包含 T_mat{1,2}（原始对齐状态下）
%   T_new        : 4x4 矩阵，sMRI1 将要使用的新的体素→世界变换
%
% 输出：
%   sMRI2_new    : 更新 T_mat 后的 sMRI2，data 和 voxsize 不变，
%                  仅 T_mat 和自动计算的 size/shape 等字段更新
%
% 原理：
%   原始关系：世界坐标 w = T1 * X1 = T2 * X2
%   得相对变换：X2 = T2 \ T1 * X1 = M * X1,  M = T2 \ T1
%   新关系要求：T2_new \ T_new * X1 = M * X1
%   => T2_new = T_new * M = T_new * (T2 \ T1)
%   由于 M 更好表示为 T1 \ T2，计算 T2_new = T_new * (T1 \ T2)

if ~isfield(sRef, 'T_mat') || isempty(sRef.T_mat) || ...
   ~isfield(sMRI, 'T_mat') || isempty(sMRI.T_mat)
    error('输入结构体必须包含有效的 T_mat{1,2}');
end

T1 = sRef.T_mat{1,2};
T2 = sMRI.T_mat{1,2};

% 相对变换：从 mri1 索引到 mri2 索引的映射
M = T1 \ T2;          % 等价于 inv(T1) * T2

% 新的 mri2 变换
T2_new = T_new * M;

% 复制结构体并更新
sMRI_new = sMRI;
sMRI_new.T_mat{1,2} = T2_new;

% 调用 MRI_read_Data 重新计算 size, shape, voxsize 等辅助字段
sMRI_new = MRI_read_Data(sMRI_new);
end

function nMRI = imresize3_mri(sMRI,varargin)
options_default = struct('new_size',1,'interp_type','linear');
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if isequal(new_size,1)
    nMRI = sMRI;
else
    raw_sz = size(sMRI.data);
    if length(new_size) == 1
        new_size = new_size.*raw_sz(1:3);
    end

    if length(raw_sz)==3
        raw_sz(4) = 1;
    end

    % if length(new_size) == 1 || ~isequal(new_size,round(new_size)) || mean(new_size)<5
    %     new_size = round(new_size.*raw_sz(1:3));
    % end

    new_size = round(new_size);
    nMRI = zeros([new_size,raw_sz(4)],'single');
    for i = 1:raw_sz(4)
        nMRI(:,:,:,i) = imresize3(sMRI.data(:,:,:,i),new_size,interp_type);
    end
    nMRI = MRI_read_Data(nMRI);

    new_sz = size(nMRI.data);
    if length(new_sz)==2
        new_sz(3) = 1;
    end
    nMRI.voxsize(1:3) = sMRI.voxsize(1:3).*raw_sz(1:3)./new_sz(1:3);
    nMRI = MRI_read_Data(nMRI);

    res = raw_sz(1:3)./new_sz(1:3);

    res = res([2,1,3]);
    % nMRI.T_mat{1,2}(1:3,4) = sMRI.T_mat{1,2}(1:3,4)./res';
    % try
    %     for i = 1:3
    %         nMRI.T_mat{1,2}(i,i) = sMRI.T_mat{1,2}(i,i)*res(i);
    %     end
    % end

    try
        temp = sMRI.T_mat{1,2}*diag([res,1]);
        nMRI.T_mat{1,2} = temp;
        nMRI.T_mat{1,2}(4,:) = sMRI.T_mat{1,2}(4,:);
    end
end
end

function nMRI = rescale_mri(sMRI, varargin)
% 对 MRI 体积进行各向同性或各向异性物理缩放，保持体素矩阵大小不变，
% 锚定原点（第一个体素角的世界坐标）不变。
%
% 输入：
%   sMRI  - MRI 结构体（需包含 data, voxsize, T_mat）
%   scale - 缩放因子，标量或 [sx, sy, sz]（分别对应 x,y,z 方向）
%            sx>0 表示 x 方向物理尺寸变为原来的 sx 倍
% 输出：
%   nMRI - 缩放后的结构体（data 不变，voxsize 和 T_mat 已更新）

options_default = struct('scale', 1);
[options, eval_str] = resolve_input(options_default, varargin);
eval(eval_str);

nMRI = sMRI;   % 复制原始数据（体素矩阵不变）

% 参数解析
if isscalar(scale)
    scale = [scale, scale, scale];
end
if numel(scale) ~= 3
    error('scale 必须为标量或包含 3 个元素对应于 x,y,z 方向。');
end
sx = scale(1);   % x 方向缩放因子
sy = scale(2);   % y 方向
sz = scale(3);   % z 方向

% 1. 更新体素尺寸（voxsize 顺序为 [dy, dx, dz]）
nMRI.voxsize(1) = sMRI.voxsize(1) * sy;   % dy
nMRI.voxsize(2) = sMRI.voxsize(2) * sx;   % dx
nMRI.voxsize(3) = sMRI.voxsize(3) * sz;   % dz

% 2. 更新仿射矩阵：右乘缩放矩阵，保持平移不变（固定原点）
if isfield(nMRI, 'T_mat') && iscell(nMRI.T_mat) && size(nMRI.T_mat,2) >= 2
    T = nMRI.T_mat{1,2};
    if ~isempty(T)
        % T 的列对应 x,y,z，右乘 diag([sx, sy, sz, 1]) 仅缩放对角线，
        % 平移 T(4,1:3) 不受影响
        nMRI.T_mat{1,2} = T * diag([sx, sy, sz, 1]);
    end
end

% 3. 刷新 size, shape 等辅助字段
nMRI = MRI_read_Data(nMRI);
end

function output = sample_mri(sMRI,varargin)
options_default = struct('pos',[],'method',1,'refine',1,'tex_sz',[],'scale',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if ~ischar(method)
    switch method
        case 0
            method = 'nearest';
        case 1
            method = 'linear';
        case 3
            method = 'cubic';
        otherwise
            method = 'linear';
    end
end

if scale ~= 1
    cen = mean(pos,1,'omitnan');
    pos = cen + scale.*(pos - cen);
end

if refine>1
    if isempty(tex_sz)
        N = sqrt(size(pos,1)); M = N;
    else
        N = tex_sz(1); M = tex_sz(2);
    end
    Xs = reshape(pos(:,1),N,M); Ys = reshape(pos(:,2),N,M); Zs = reshape(pos(:,3),N,M);
    % Xsk = imresize(Xs,refine); Ysk = imresize(Ys,refine); Zsk = imresize(Zs,refine);
    Xsk = imresize_grid(Xs,refine); Ysk = imresize_grid(Ys,refine); Zsk = imresize_grid(Zs,refine);
    pos = [Xsk(:),Ysk(:),Zsk(:)];
end

pos = pos./sMRI.voxshape;
index = interp3(single(sMRI.data),pos(:,1),pos(:,2),pos(:,3),method,0);

output = {index,pos};
end

function nMRI = cat_mri(sMRI,varargin)
options_default = struct('dim',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

nMRI = sMRI{1};
mri = fcellfun(@(x) x.data, sMRI);
nMRI.data = cat(dim,mri{:});
nMRI = MRI_read_Data(nMRI);
end

function nMRI = permute_mri(sMRI,varargin)
options_default = struct('dim',[3,1,2]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

nMRI = sMRI;
nMRI.data = permute(sMRI.data,[dim,4]);
nMRI.voxsize = sMRI.voxsize(dim);
T_mat = nMRI.T_mat{2};
T_mat = T_mat([dim,4],[dim,4]);
nMRI.T_mat{2} = T_mat;
nMRI = MRI_read_Data(nMRI);
end

function nMRI = make_grid(sMRI,varargin)
mode = '';
options_default = struct('mode','block','d',10);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

sz = sMRI.size;
[X,Y,Z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
A = mod(mod(floor(X/d),2)+mod(floor(Y/d),2)+mod(floor(Z/d),2),2);
switch mode
    case 'block'

    case 'line'
        A = edge3(A,'approxcanny',0.6);
end
nMRI = sMRI; nMRI.data = A;
end

function nMRI = make_bubbles(sMRI,varargin)
options_default = struct('pos',[],'r',1,'bin',1,'refine',1,'rm_border',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

nMRI = sMRI;
if refine>1
    nMRI = MRI_func_ROI(nMRI,'imresize3',refine);
end
V = single(0*nMRI.data);
pos = round(pos*refine);

if length(bin) == 1
    pos = pos(1:bin:end,:);
else
    N = sqrt(size(pos,1)); % 需要是矩形纹理采样
    if length(bin)  == 2
        bin([3,4]) = [1,1];
    end
    if length(bin)  == 3
        [bin(3),bin(4)] = ind2sub(bin([1,2]),bin(3));
    end
    idx = reshape(1:N^2,N,N);

    % 去除边缘的重叠
    idx([1:rm_border,end-rm_border+1:end],:) = 0;
    idx(:,[1:rm_border,end-rm_border+1:end],:) = 0;

    idx = idx(bin(3):bin:end,bin(4):bin:end);
    pos = pos(nonzeros(idx(:)),:);
end

for i = 1:size(pos,1)
    V(pos(i,2),pos(i,1),pos(i,3)) = 1;
end
H = fspecial3('ellipsoid',r*refine); H = size(H,1)^2*H;
V = imfilter(V,H,'replicate');
nMRI.data = V;
end

function T_in = smooth_field(T_in,varargin)
options_default = struct('smoothness',[],'smooth_mask',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if smoothness == 0
    smoothness = [];
end

% smooth_func = @(x,sm) smooth3(x,'gaussian',3*sm,sm);
smooth_func = @(x,sm) imgaussfilt3(x,sm);

if ~isempty(smoothness)
    if length(smoothness)>1
        smooth_num = smoothness(2);
    else
        smooth_num = 1;
    end
    smoothness = smoothness(1);
    if ~isempty(smooth_mask)
        if ndims(smooth_mask) < 3
            smooth_mask = T_in.data>0;
        end
        smooth_mask = double(smooth_mask);
        for i_dim = 1:3
            smooth_mask(:,:,:,i_dim) = smooth_func(smooth_mask(:,:,:,i_dim),smoothness);
        end
    else
        smooth_mask = 1;
    end
    for i = 1:smooth_num
        for i_dim = 1:3
            T_in.data(:,:,:,i_dim) = smooth_func(T_in.data(:,:,:,i_dim),smoothness);
        end
        if mod(i,2) == 1
            T_in.data = T_in.data./(1e-6+smooth_mask);
        end
    end
end
end

function inner_GM = make_subcort(sMRI,varargin)
options_default = structure('mask',[],'GM',[],'sharpen',[5,2],'fuse_rate',0.3);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

vol = amp(sMRI.data,[0 1e-3]);
inner_gray = 1 - (1-vol).*double(mask);

if ~isempty(GM)
    GM = amp(GM,1e-3);
    inner_GM = (1 - GM.*double(mask)).*inner_gray;
    inner_GM = fuse_rate*inner_GM + (1-fuse_rate)*inner_gray;
else
    inner_GM = inner_gray;
end

for i = 1:size(inner_GM,3)
    inner_GM(:,:,i) = imsharpen(inner_GM(:,:,i),'Radius',sharpen(1),'Amount',sharpen(2));
end
% inner_GM = imresize3(inner_GM,2);
inner_GM(inner_GM<0) = 0; inner_GM(inner_GM>1) = 1;
end

function [dist_map,bnd_box] = dist_map(in_mask,varargin)
options_default = structure('if_inv',0,'ext',-1,'method','bwdist');
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

dist_map = in_mask;
in_mask = in_mask.data;

if ext < 0
    out_box = in_mask;
else
    bnd_box = zeros(1,6);
    for i = 1:3
        temp = sum(1-in_mask,setdiff(1:3,i));
        temp = find(temp>0);
        bnd_box(i) = temp(1); bnd_box(i+3) = temp(end)-temp(1)+1;
    end

    if isempty(ext)
        if if_inv
            ext = 5;
        else
            ext = 0;
        end
    end
    bnd_box(1:3) = bnd_box(1:3) - ext;
    bnd_box(4:6) = bnd_box(4:6) + 2*ext;
    out_box = in_mask(bnd_box(1)+[1:bnd_box(4)],bnd_box(2)+[1:bnd_box(5)],bnd_box(3)+[1:bnd_box(6)]);
end

if if_inv
    out_box = ~out_box;
end

switch method
    case 'bwdist'
        temp = bwdist(out_box);
    case 'bwdistsc'
        temp = bwdistsc(out_box);
    case 'knnsearch'
        [temp,idx] = bwdist_knnsearch(out_box);
end

if ext >= 0
    out_dist = 0*in_mask;
    out_dist(bnd_box(1)+[1:bnd_box(4)],bnd_box(2)+[1:bnd_box(5)],bnd_box(3)+[1:bnd_box(6)]) = temp;
else
    out_dist = temp;
end

dist_map.data = out_dist;
end

function [distmap,idx] = bwdist_knnsearch(mask)
sz = size(mask);
[X,Y,Z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
Vall = [X(:),Y(:),Z(:)];

mask = single(mask);
mask_edge = mask.*(imgaussfilt3(mask,0.25)<(1-1e-6));
Vtarg = Vall(mask_edge>0,:);
Vq = Vall(mask==0,:);

[idx,dist] = knnsearch(Vtarg,Vq);

distmap = 0*mask;
distmap(mask==0) = dist;
end

function sMRI_fuse = fuse_rgb(sMRI,varargin)
alpha = 0;
options_default = structure('heat',[],'heat_cmap',jet,'alpha',0.2,'heat_clim',[0 1],'base_cmap',[],'base_clim',[0 1]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

heat_rgb = heat2rgb(heat,heat_clim,heat_cmap);
if isempty(base_cmap)
    mri_rgb = amp_auto(sMRI.data,base_clim);
else
    mri_rgb = heat2rgb(sMRI.data,base_clim,base_cmap);
end

sMRI_fuse = sMRI;
sMRI_fuse.data = single(mri_rgb).*(1-alpha) + single(heat_rgb).*alpha;
end

function heat_rgb = heat2rgb(heat,heat_clim,heat_cmap)
heat_rgb = round(amp_auto(heat,heat_clim).*size(heat_cmap,1));
heat_rgb = ind2rgb(heat_rgb(:),heat_cmap);
heat_rgb = reshape(heat_rgb,[size(heat),3]);
end

function B = amp_auto(A,lim)
if max(lim(:))<=1
    B = amp(A,lim);
else
    if length(lim)==1
        lim(2) = max(A(:));
    end
    B = (double(A)-lim(1))/((lim(2)-lim(1)));
    B(B<0) = 0; B(B>1) = 1;
end
end