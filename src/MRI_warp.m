function [vout,info] = MRI_warp(vMap,vsubj,mode,varargin)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin<3
    mode = 'auto';
end

all_mode = {'BST','elastix','bspline'};
if ~any(contains(all_mode,mode))
    if ~isequal(mode,'auto')
        varargin = [{mode},varargin];
    end
    
    if isstruct(vMap) && isfield(vMap,'type')
        mode = vMap.type;
    else
        if isstruct(vMap) && isfield(vMap,'data')
            mode = 'BST';
        end
        if isstruct(vMap) && isfield(vMap,'TransformParameters')
            mode = 'elastix';
        end
        if ischar(vMap) && contains('TransformParameters')
            mode = 'elastix';
        end
    end

    if isnumeric(vMap) && ndims(vMap)>2
        mode = 'BST';
    end
end

info = [];
vsubj = MRI_read_Data(vsubj,'verbose',0);
vout = vsubj;
switch mode
    case 'affine'
        tform = affinetform3d(vMap);
        Rout = affineOutputView(size(vsubj.data),tform,"BoundsStyle","SameAsInput");
        vout.data = imwarp(vsubj.data,tform,'OutputView',Rout);
    case 'BST'
        [vout,info] = svreg_apply_map(vMap,vsubj,varargin{:});
    case 'elastix'
        vout = MRI_transformix(vMap,vsubj,varargin{:});
    case 'bspline'
        [vout,info.T] = bspline_transform(vMap.T,vsubj.data,vMap.spacing,1);

        % vout = vsubj;
        % fBox = round(vMap.Vm_min - vMap.ext);
        % fBox = fBox + [0,0,0; vMap.box-1];
        % 
        % crop_data = vsubj.data(fBox(1,2):fBox(2,2),fBox(1,1):fBox(2,1),fBox(1,3):fBox(2,3));
        % 
        % % [warp_data,info.T] = bspline_transform(vMap.T,crop_data,vMap.spacing,1);
        % T = permute(vMap.T,[2,1,3,4]);
        % for i = 1:3
        %     Tbig(:,:,:,i) = imresize3(T(:,:,:,i),size(T,1:3).*vMap.spacing);
        % end
        % warp_data = interp3(single(crop_data),Tbig(:,:,:,1),Tbig(:,:,:,2),Tbig(:,:,:,3));
        % 
        % mBox = round(vMap.Vf_min - vMap.ext);
        % mBox = mBox + [0,0,0; size(Tbig,[2,1,3])-1];
        % 
        % vout = vsubj;
        % vout.data(mBox(1,2):mBox(2,2),mBox(1,1):mBox(2,1),mBox(1,3):mBox(2,3)) = warp_data;

    otherwise
        error
end
end

function [vout,vMap] = svreg_apply_map(vMap,vsubj,varargin)
% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy
% raw params: map_file,data_file,out_file,target_file,smoothness,datatype,bitpix,interp_type
ref = '';
options_default = structure('ref',[],'smoothness',0,'smooth_mask',[],'interp_type','linear','remove_out',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if ~isstruct(vMap)
    vMap = MRI_read_Data(vMap);
end
if ~isstruct(vsubj)
    vsubj = MRI_read_Data(vsubj);
end
if isempty(ref)
    % ref = vsubj;
    ref = vMap.data(:,:,:,1);
end
if ~isempty(ref) && ~isstruct(ref)
    ref = MRI_read_Data(ref);
end

if ~ischar(interp_type)
    switch interp_type
        case 0
            interp_type = 'nearest';
        case 1
            interp_type = 'linear';
        case 3
            interp_type = 'cubic';
        otherwise
            interp_type = 'linear';
    end
end

vMap.data = squeeze(vMap.data);
use_2d = ismatrix(vsubj.data) || size(vMap.data,3)==2;

if ~use_2d
    vMap = MRI_func_ROI(vMap,'smooth_field',smoothness,smooth_mask);
    vMap.data = single(vMap.data);
    
    vt = ref;
    vout = ref;
    vout.data = single(vout.data);
    for i = 1:size(vsubj.data,4)
        vout.data(:,:,:,i) = apply_warp(single(vsubj.data(:,:,:,i)),vMap,interp_type);
    end
else
    vt=ref;
    vout = ref;
    vout.data = single(vout.data);
    for i = 1:size(vsubj.data,3)
        vout.data(:,:,i) = interp2(vsubj.data(:,:,i),vMap.data(:,:,1),vMap.data(:,:,2),interp_type);
    end
end

if ~isempty(smoothness) && remove_out
    vout.data=double(vout.data).*double(vt.data>0);
end
end

function out = apply_warp(data,vMap,interp_type)
    sz = size(vMap.data);
    if sz(3)<1000
        out = interp3(data,vMap.data(:,:,:,1),vMap.data(:,:,:,2),vMap.data(:,:,:,3),interp_type);
    else
        z = round(sz(3)/2);
        out1 = interp3(data,vMap.data(:,:,1:z,1),vMap.data(:,:,1:z,2),vMap.data(:,:,1:z,3),interp_type);
        out2 = interp3(data,vMap.data(:,:,z+1:end,1),vMap.data(:,:,z+1:end,2),vMap.data(:,:,z+1:end,3),interp_type);
        out = cat(3,out1,out2);
    end
end

function svreg_apply_map_forDTI(map_file,data_file,out_file,target_file,smoothness,datatype,bitpix,interp_type)
% source code: svreg_apply_map.m
in_base=tempname;%
%   subbasename=data_file(1:strfind(data_file,'.dwi.')-1);
out_base=[in_base,'.atlas'];
%  out_base_orig=[subbasename,'.atlas'];
eig2nifti(data_file,in_base);
warp_DTI(map_file, [in_base '.L1.nii.gz'], [in_base '.L2.nii.gz'], [in_base '.L3.nii.gz'], [in_base '.V1.nii.gz'], [in_base '.V2.nii.gz'], [in_base '.V3.nii.gz'], out_base);
nifti2eig([out_base '.L1.nii.gz'], [out_base '.L2.nii.gz'], [out_base '.L3.nii.gz'], [out_base '.V1.nii.gz'], [out_base '.V2.nii.gz'], [out_base '.V3.nii.gz'], out_file);
%     v=load_untouch_eig_gz(data_file);
%     [interpL1, interpL2, interpL3, interpV1, interpV2, interpV3] = interp_dti_data(vMap.data, v.data(:,:,:,1), v.data(:,:,:,2), v.data(:,:,:,3), v.data(:,:,:,4:6), v.data(:,:,:,7:9), v.data(:,:,:,10:12), false);
%     regMap=load_nii_BIG_Lab(map_file);
%     [J1,J2,J3] = jacobian_nii(regMap, 'no_smooth');
%     [W1,W2,W3] = PPD(interpV1, interpV2, interpV3, J1, J2, J3);
%     generate_eig_file(interpL1, interpL2, interpL3, W1, W2, W3, out_file);
end

%%
function [sTrans] = MRI_transformix(trans_params,trans_file,varargin)
ref = '';
options_default = structure('mov',[],'ref',[],'remake_sz',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if remake_sz
    remake_size = round(mov.size.*mov.voxsize/min(mov.voxsize));
    ref_zoom = MRI_func_ROI(ref,'imresize3',remake_sz);
    trans_file = remake_mri(trans_file,remake_size);
end

if isstruct(trans_params)
    trans_params = trans_params.reg_parafile;
end

[trans_result] = Mtransformix(trans_file, trans_params);

if remake_sz
    sTrans = ref_zoom; sTrans.data = imresize3(trans_result,ref_zoom.size);
else
    sTrans = sFix; sTrans.data = trans_result;
end
end

function sMRI = remake_mri(sMRI,remake_sz)
sMRI.data = imresize3(sMRI.data,remake_sz); sMRI.voxsize = [1,1,1]; sMRI.T_mat = [];
sMRI = MRI_read_Data(sMRI);
end