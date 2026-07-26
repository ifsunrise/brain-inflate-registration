function [sReg,trans_params,def_field,inv_field,sMov,sFix] = MRI_register(sMov,sFix,mode,varargin)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
all_mode = {'elastix','demons'};
if ~contains(all_mode,mode)
    varargin = [{mode},varargin];
    mode = 'elastix';
end

if ~isstruct(sMov)
    sMov = MRI_read_Data(sMov);
end
if ~isstruct(sFix)
    sFix = MRI_read_Data(sFix);
end

switch mode
    case 'elastix'
        [sReg,trans_params,def_field,inv_field,sMov,sFix] = MRI_elastix(sMov,sFix,varargin{:});
    case 'demons'
        [sReg,def_field] = MRI_demons(sMov,sFix,varargin{:});
        trans_params = [];
end
end

function [sReg,trans_params,def_field,inv_field,sMov,sFix] = MRI_elastix(sMov,sFix,varargin)
options_default = structure('reg_type','nonrigid','auto_params',1,'reg_params',[],'movingPts',[],'fixedPts',[],'ptsWeight',[],'remake_sz',0,'format','',...
    'get_field',[0,0],'field_type','xyz','warp_byField',0,'get_inv',0,'field_fix_bnd',0,'interp_type','linear','smoothness',0,'verbose',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

elastix_root = replace(which('Melastix'),'Melastix_code\Melastix.m','');

if isempty(format)
    if ismatrix(sMov.data) || ~isempty(movingPts)
        format = '.tif';
    else
        format = '.nii.gz';
    end
end

remake_size = [];
sFix_raw = sFix;
if remake_sz
    sMov_raw = sMov;
    sFix_zoom = MRI_func_ROI(sFix,'imresize3',remake_sz,'interp_type',interp_type); % 防止高分辨率的sMov损失
    remake_size = round(sMov.size.*sMov.voxsize/min(sMov.voxsize));

    sMov = remake_mri(sMov,remake_size); sFix = remake_mri(sFix,remake_size);
end
movFile = [elastix_root,'\Images\mov',format]; MRI_save_File(sMov,movFile)
fixFile = [elastix_root,'\Images\fix',format]; MRI_save_File(sFix,fixFile)

% E:\work_Code\Matlab\ZYS\DIY_App\Elastix_toolbox\Parameter
% (FixedInternalImagePixelType "short")
% (MovingInternalImagePixelType "short")
% (FixedImageDimension 3)
% (MovingImageDimension 3)
% (UseDirectionCosines "true")
% (Registration "MultiResolutionRegistration")
% (Interpolator "BSplineInterpolator")
% (ResampleInterpolator "FinalBSplineInterpolator")
% (OpenCLDeviceID "0")
% (OpenCLDeviceType "GPU")
% (Resampler "OpenCLResampler")
% (OpenCLResamplerUseOpenCL "true")
% //(FixedImagePyramid "OpenCLFixedGenericImagePyramid")
% //(OpenCLFixedGenericImagePyramidUseOpenCL "true")
% //(MovingImagePyramid "OpenCLMovingGenericImagePyramid")
% //(OpenCLMovingGenericImagePyramidUseOpenCL "true")
% (FixedImagePyramid "FixedRecursiveImagePyramid")
% (MovingImagePyramid "MovingRecursiveImagePyramid")
% (Optimizer "AdaptiveStochasticGradientDescent")
% (Transform "BSplineTransform")
% (Metric "DisplacementMagnitudePenalty")
% (AutomaticScalesEstimation "true")
% (AutomaticTransformInitialization "true")
% (HowToCombineTransforms "Compose")
% (UseFastAndLowMemoryVersion "true")
% (NumberOfHistogramBins 32)
% (ErodeMask "false")
% (NumberOfResolutions 7)
% (FinalGridSpacingInPhysicalUnits 19 19 19)
% (MaximumNumberOfIterations 300)
% (NumberOfSpatialSamples 1000)
% (SampleRegionSize 19 19 19)
% (NewSamplesEveryIteration "true")
% (ImageSampler "Random")
% (BSplineInterpolationOrder 3)
% (FinalBSplineInterpolationOrder 3)
% (DefaultPixelValue 0)
% (WriteResultImage "false")
% (WriteTransformParametersEachIteration "false")
% (WriteTransformParametersEachResolution "false")
% (ShowMetricValue "false")
% (ShowExactMetricValue "false" "false" "false" "false")
% (WriteResultImageAfterEachResolution "false")
% (ResultImagePixelType "short")
% (ResultImageFormat "nii.gz")
% (Metric0Weight 1)

if auto_params
    if auto_params == 1
        base_sz = 20;
    else
        base_sz = auto_params;
    end
    grid_size = round(mean(sMov.size)/base_sz)*[1,1,1];
    %NumberOfSpatialSamples =
    samp_size = round(mean(sMov.size)/base_sz)*[1,1,1];

    if iscell(reg_params)
        temp = [];
        for i = 1:2:length(reg_params)
            temp.(reg_params{i}) = reg_params{i+1};
        end
        reg_params = temp;
    end
    params = reg_params;

    reg_params.FinalGridSpacingInPhysicalUnits = grid_size;
    reg_params.SampleRegionSize = samp_size;
    reg_params.ResultImageFormat = format(2:end);

    if isstruct(params)
        fields = fieldnames(params);
        for i = 1:length(fields)
            field = fields{i};
            reg_params.(field) = params.(field);
        end
    end

    if ~isempty(movingPts) && ~isempty(ptsWeight)
        reg_params.Metric0Weight = ptsWeight(1);
        reg_params.Metric1Weight = ptsWeight(2);
    end

    if verbose
        disp(reg_params)
    end
end

[mov_reg,trans_params,def_field] = Melastix(movFile,fixFile,[],'reg_type',reg_type,'reg_params',reg_params,...
    'movingPts',movingPts,'fixedPts',fixedPts,'dim',ndims(sMov.data),'get_field',get_field,'get_inv',get_inv,'verbose',verbose);
verbose = verbose>=1;
%% 将相对位移的变形场变为绝对的坐标图
if ~isempty(def_field)
    switch field_type
        case 'xyz'
            % 1. 将原始 def_field 数据转为坐标场（基于配准时网格）
            temp = def_field.data;
            if remake_sz
                sFix_used = sFix;      % 重采样后的固定图像（remake_size）
            else
                sFix_used = sFix_raw;  % 原始固定图像
            end
            def_field = sFix_used;
            def_field.data = temp;
            [def_field] = Melastix_getfield(def_field, sMov, sFix_used, sMov, field_fix_bnd);

            % 2. 若启用了 remake_sz，直接将坐标场放大并映射到原始移动图像
            if remake_sz
                % 尺寸
                sz_zoom = sFix_zoom.size;      % 目标固定图像尺寸 [ny,nx,nz]
                
                % 对坐标场的三个通道分别放大到 sz_zoom
                vy = imresize3(def_field.data(:,:,:,1), sz_zoom, 'linear');
                vx = imresize3(def_field.data(:,:,:,2), sz_zoom, 'linear');
                vz = imresize3(def_field.data(:,:,:,3), sz_zoom, 'linear');
                
                % 移动图像缩放因子（从配准尺寸到原始高分辨率尺寸）
                scale_mov = sMov_raw.size ./ sMov.size;   % [ny,nx,nz]
                
                % 坐标值映射回原始移动图像坐标系
                def_field = sFix_zoom;
                def_field.data = cat(4, vy * scale_mov(1), ...
                                        vx * scale_mov(2), ...
                                        vz * scale_mov(3));
            end

            def_field.remake_size = remake_size;
    end
end

if warp_byField
    mov_reg = MRI_warp(def_field,sMov,'BST','ref',sFix_raw,'interp_type',interp_type,'smoothness',smoothness); mov_reg = mov_reg.data;
end

inv_field = [];
if get_inv
    inv_parafile = trans_params.inv_parafile;
    inv_field = Mtransformix_pts('all', inv_parafile,'verbose',verbose);
    [inv_field] = Melastix_getfield(inv_field,sFix_raw,sMov);

    sTran = MRI_warp(inv_field,sFix_raw,'BST','ref',sMov,'interp_type',interp_type);
end

if remake_sz
    sReg = sFix_zoom; sReg.data = imresize3(mov_reg,sFix_zoom.size);
else
    sReg = sFix; sReg.data = mov_reg;
end
end

function sMRI = remake_mri(sMRI,remake_sz)
sMRI.data = imresize3(sMRI.data,remake_sz);
sMRI.voxsize = [1,1,1]; sMRI.T_mat = [];
sMRI = MRI_read_Data(sMRI);
end

function inv_def()
map.img = double(def_field.data);

src_mask.img = invmap_mask;
tar_mask.img = imgaussfilt3(sMask.data,3)>0;

size_tar=size(tar_mask.img);
size_src=size(src_mask.img); size_src=size_src(1:3);
tar_mask_ind=find(tar_mask.img>0);src_mask_ind=find(src_mask.img>0);

mapx=squeeze(map.img(:,:,:,1));mapx=mapx(src_mask_ind);
mapy=squeeze(map.img(:,:,:,2));mapy=mapy(src_mask_ind);
mapz=squeeze(map.img(:,:,:,3));mapz=mapz(src_mask_ind);

[s{1},s{2},s{3}]=ind2sub(size_src,src_mask_ind);
[tX,tY,tZ]=ind2sub(size_tar,tar_mask_ind);

inv_map{1}=single(zeros(size_tar(1),size_tar(2),size_tar(3)));
inv_map{2}=inv_map{1};inv_map{3}=inv_map{1};
clear src_mask src_mask_ind
for jj=1:3
    inv_map{jj}(tar_mask_ind)=mygriddata3_IntNearestHack(mapx,mapy,mapz,s{jj},tX,tY,tZ);
end
end

function [sReg,map] = MRI_demons(sMov,sFix,varargin)
box = [];
options_default = structure('iter',20,'lv',1,'box',[],'smooth_lv',1,'use_gpu',1,'data_type','single','bin',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if isempty(box)
    Im = sMov.data;
    If = sFix.data;
else
    Im = sMov.data(box(1):box(4),box(2):box(5),box(3):box(6));
    If = sFix.data(box(1):box(4),box(2):box(5),box(3):box(6));
end

Im(isnan(Im)) = 0;
If(isnan(If)) = 0;

switch data_type
    case 'single'
        Im = single(Im);
        If = single(If);
    case 'uint16'
        Im = uint16(Im);
        If = uint16(If);
end

if (bin ~= 1)
    Im = imresize3(Im,1/bin);
    If = imresize3(If,1/bin);
end

if (use_gpu>0)
    [D,moving_reg] = imregdemons(gpuArray(Im),gpuArray(If),iter,'PyramidLevels',lv,'AccumulatedFieldSmoothing',smooth_lv);
    moving_reg = gather(moving_reg);
    D = gather(D);
else
    [D,moving_reg] = imregdemons(Im,If,iter,'PyramidLevels',lv,'AccumulatedFieldSmoothing',smooth_lv);
end

sReg = sMov;
if ~isempty(box)
    sReg.data(box(1):box(4),box(2):box(5),box(3):box(6)) = moving_reg;
else
    sReg.data = moving_reg;
end

% sReg.data = imwarp(Im,D);

sz = sFix.size;
[X,Y,Z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
D0 = cat(4,X,Y,Z);

if (bin ~= 1)
    Drsz = D0*0;
    for i = 1:3
        Drsz(:,:,:,i) = imresize3(D(:,:,:,i), sz(1:3));
    end
    D = Drsz*bin;
end

map = D + D0;

if (bin ~= 1)
    sReg = MRI_warp(map,sMov,'BST','ref',sFix);
end

map = MRI_read_Data(map);
map.data = single(map.data);
end