function [sMRI] = MRI_read_Data(mriFile,varargin)
%UNTITLED 此处显示有关此函数的摘要
%   MRI .nii 文件要素: T_mat， 对角线尺寸参数需要与分辨率对应，最右边平移参数是物理坐标
options_default = structure('mode','bst','voxsize',[],'val_offset',0,'T_serial',[],'verbose',1,'if_trans',0,'if_permute',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if ischar(mriFile) | isstring(mriFile)
    if ~isfile(mriFile)
        disp(['error: file ',mriFile,' not exist!'])
        error
    end

    mriFile = char(mriFile);
    [~,~,file_format] = fileparts(mriFile);
else
    if isnumeric(mriFile) | islogical(mriFile)
        file_format = 'numeric';
    else
        file_format = class(mriFile);
    end
end
if strcmp(file_format,'.gz')
    file_format = '.nii';
    if ~isfile(mriFile)
        mriFile = mriFile(1:end-3);
    end
end
if strcmp(file_format,'.nii')
    if ~isfile(mriFile)
        mriFile = [mriFile,'.gz'];
    end
end

mode = EEG_replace_acronyms(mode);
T_one = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

sMRI = [];
switch file_format
    case '.tif'
        sMRI.data = loadtiff(mriFile);
        tiff_info = imfinfo(mriFile);
        tiff_info = tiff_info(1);
        pixsize = 1./[tiff_info.YResolution,tiff_info.XResolution];
        if isempty(pixsize)
            pixsize = [1,1];
        end
        
        try
            eval_str = splitlines(upper(tiff_info.ImageDescription));
            for i = 1:length(eval_str)
                try
                    eval([eval_str{i},';']);
                end
            end
        end
        
        try
            spacing_z = SPACING;
        catch
            spacing_z = 1;
        end
        
        sMRI.voxsize = 0.001*[pixsize,spacing_z];
    case 'numeric'
        if isempty(voxsize)
            voxsize = [1,1,1];
        end

        sMRI.data = mriFile;
        sMRI.name = 'numeric';
%         sMRI.size = size(sMRI.data);
        sMRI.voxsize = voxsize;
        T_mat = T_one;
        T_mat([1,6,11]) = voxsize([2,1,3]);
        T_mat(4,1:3) = -0.5*(size(sMRI.data,[2,1,3]).*voxsize([2,1,3]))';
        sMRI.T_mat{2} = T_mat;
        sMRI.val_offset = 0;
        val_offset = 0;
        
        verbose = 0;
    case 'struct'
        sMRI = mriFile;
        val_offset = 0;
        verbose = 0;
    otherwise
        try
            BST_toolbox_path = which('bst_get');
            if strcmp(mode,'bst') && isempty(BST_toolbox_path)
                mode = 'native';
            end
            switch mode
                case 'bst'
                    BST_toolbox_path = replace(BST_toolbox_path,'core\bst_get.m','');
                    BST_io_path = [BST_toolbox_path,'io'];
                    addpath(BST_io_path)
                    
                    iSubject = 0;
                    FileFormat = 'ALL';
                    isInteractive = 0;
                    isAutoAdjust = 0;
                    Comment = [];
                    Labels = [];
                    clog = evalc('[sMRI_BST] = import_mri_(iSubject, mriFile, FileFormat, isInteractive, isAutoAdjust, Comment, Labels);');
                    
                    sMRI.data = sMRI_BST.Cube;
                    sMRI.name = sMRI_BST.Comment;
                    sMRI.voxsize = sMRI_BST.Voxsize;
                    sMRI.T_mat = sMRI_BST.InitTransf;
                    sMRI.hdr = sMRI_BST.Header;
                    
                    % sMRI.val_offset = sMRI_BST.Header.dim.funused2;

                    % reorient
                    % if size(sMRI.T_mat,1)==2 && strcmp(sMRI.T_mat{2,1},'reorient')
                    %     [~,reo_dims] = max(abs(sMRI.T_mat{2,2}),[],1);
                    %     reo_dims(end) = [];
                    %     voxshape_raw = sMRI.hdr.dim.pixdim(2:4);
                    %     sMRI.voxshape = voxshape_raw(reo_dims);
                    % end

                    if size(sMRI.T_mat,1)>1
                        sMRI.T_mat{2,2} = sMRI.T_mat{2,2}';
                    end

                case 'native'
                    error
                otherwise
                    error
            end
        catch read_error
            if strcmp(file_format,'.nii')
                if ~strcmp(mode, 'native') * verbose==1
                    disp(['error [',mode,'] mode, try [nifttiread] => ', read_error.message])
                end
                %             if strcmp(read_error.identifier,'MATLAB:structRefFromNonStruct')
                try
                    % !!!! native方法不能执行T_mat !!!!
                    nii_info = niftiinfo(mriFile);
                    sMRI.data = squeeze(niftiread(mriFile));
                    sMRI.voxsize = squeeze(nii_info.PixelDimensions);
                    % sMRI.voxsize = sMRI.voxsize([2,1,3]);
                    sMRI.T_mat{2} = nii_info.Transform.T;
                    sMRI.val_offset = 0;
                    sMRI.nifti = nii_info;

                    pixdim = nii_info.raw.pixdim;
                    % if pixdim(1) == -1 && if_permute
                    %     dims = [1,3,2];
                    %     sMRI.data = permute(sMRI.data,dims);
                    %     sMRI.voxsize(1:3) = sMRI.voxsize(dims);
                    %     sMRI.T_mat{2} = sMRI.T_mat{2}([dims,4],[dims,4]);
                    % end
                end
            end
        end
        
        if isequal(T_serial,[])
            T_serial = {'dim(213)'};
        end
end

if ischar(mriFile) && (endsWith(mriFile,'.nii.gz') || endsWith(mriFile,'.nii'))
    sMRI.nifti = niftiinfo(mriFile);
    ext = MRI_read_ext(mriFile);
    if ~isempty(ext)
        sMRI.nifti.ext = ext;
    end
end

%%
% if isempty(sMRI) && isfile(mriFile)
%     [mri_dir,mri_name] = fileparts(mriFile);
%     [V,info] = ReadData3D(mriFile);
%     sMRI.data = V(:,end:-1:1,:);
%     sMRI.voxsize = info.PixelDimensions;
%     sMRI.name = mri_name;
%     sMRI.T_mat = {'vox2ras',[reshape(info.TransformMatrix,3,3),[0;0;0]; 0 0 0 1]};
% end

if ~isempty(voxsize)
    sMRI.voxsize = voxsize;
end
%%

% if isfield(sMRI,'val_offset') && (sMRI.val_offset==0) && size(sMRI.data,4)==1 && ~isa(sMRI.data,'int16')
%     if isempty(val_offset)
%         baseline = min(sMRI.data(1:666:end));
%         if baseline ~= 0
%             disp("apply data = data - "+baseline)
%         end
%     else
%         baseline = val_offset;
%     end
%     try
%         sMRI.data = sMRI.data-baseline;
%     end
% end

% sMRI.hdr.datatype = class(sMRI.data);
if isa(sMRI.data,'int16')
    % if mean(sMRI.data(:))<-1000
    %     sMRI.data = single(sMRI.data)+32768;
    % end
    sMRI.data = uint16(sMRI.data);
end
if isa(sMRI.data,'int8')
    sMRI.data = uint8(sMRI.data);
end

if isfield(sMRI, 'nifti')
    % scl_slope = sMRI.hdr.nifti.scl_slope;
    % scl_inter = sMRI.hdr.nifti.scl_inter;

    scl_slope = sMRI.nifti.MultiplicativeScaling;
    scl_inter = sMRI.nifti.AdditiveOffset;

    if (scl_slope == 0)
        scl_slope = 1;
    end

    sMRI.data = single(sMRI.data) * scl_slope + scl_inter;

    % sMRI.hdr.nifti.scl_slope = 1;
    % sMRI.hdr.nifti.scl_inter = 0;

    sMRI.nifti.MultiplicativeScaling = 1;
    sMRI.nifti.AdditiveOffset = 0;
end

if isequal(T_serial,[])
    T_serial = {};
end
sMRI.T_serial = T_serial;
if iscell(T_serial)
    sMRI.data = serial_eval_func(sMRI.data,T_serial,{@dim});
else
    if isempty(T_serial)
        pycode = @serial_resolve_method;
        this = sMRI;
        code = [T_serial];
        eval(pycode(code));
        sMRI = this;
    end
end
sMRI.size = size(sMRI.data);
if length(sMRI.size) == 2
    sMRI.size = [sMRI.size,1];
end
%% Transform matrix
if ~isfield(sMRI,'T_mat')
    sMRI.T_mat = [];
end

if isempty(sMRI.T_mat)
    T_ref = zeros(4);
    T_ref([1,6,11]) = sMRI.voxsize;
    T_ref(4,4) = 1;
    temp = sMRI.voxsize.*size(sMRI.data,[1:3]);
    T_ref(4,1:3) = -0.5*temp([2,1,3]);
    sMRI.T_mat = {'vox2ras',T_ref};
end

T_ref = sMRI.T_mat{1,2};

if all(T_ref(4,1:3) == 0)
    T_ref = T_ref';
end
% sMRI.T_mat = sMRI.T_mat(1,:);

T_translate = -0.5*(size(sMRI.data,[2,1,3]))';
if isfield(sMRI,'T_mat') && ~isempty(sMRI.T_mat) && mean2((sMRI.T_mat{1,2} - T_one).^2) < 0.001
    T_ref(4,1:3) = T_translate;
end

if if_trans
    % 负的变换系数会使图像翻转
    % for i_dim = 1:3
    %     if T_ref(i_dim,i_dim)<0
    %         switch i_dim
    %             case 1
    %                 sMRI.data = sMRI.data(:,end:-1:1,:,:);
    %             case 2
    %                 sMRI.data = sMRI.data(end:-1:1,:,:,:);
    %             case 3
    %                 sMRI.data = sMRI.data(:,:,end:-1:1,:);
    %         end
    %         T_ref(i_dim,i_dim) = -T_ref(i_dim,i_dim);
    %         T_ref(4,i_dim) = -T_ref(4,i_dim);
    %     end
    % end

    tform = affinetform3d(sMRI.T_mat{2}');
    sMRI.data = imwarp(sMRI.data,tform);
    sMRI.voxsize = [1,1,1];
    T_ref = 0*T_ref;
    T_ref([1,6,11,16]) = 1;
    T_ref(4,1:3) = -0.5*(size(sMRI.data,[2,1,3]))';
end

sMRI.T_mat{1,2} = T_ref;

%%
if (isempty(T_serial) == 0) && isfile('myDiaryFile')
    log = importdata('myDiaryFile');
    log = strjoin(log,' ');
    eval(log);
    delete('myDiaryFile')
end

% y,x,z --> x,y,z
sMRI.shape = sMRI.size([2,1,3]);
if isempty(sMRI.voxsize)
    sMRI.voxsize = 1+0*sMRI.size;
end
sMRI.voxshape = sMRI.voxsize([2,1,3]);

if verbose
    fprintf("sMRI "+mode+" "+mat2STR(sMRI.size)+" res:"+mat2STR(sMRI.voxsize));
    fprintf('\n')
end
end

function str = mat2STR(mat)
str = join(split(mat2str(round(mat,2))),',');
end

function mri = dim(mri,dims)
if ischar(dims)==0
    dims = num2str(dims);
end

mri = permute(mri,[str2double(dims(1)),str2double(dims(2)),str2double(dims(3)),4]);

%fprintf(['=> apply transform: ', 'dims permute [',dims,']; '])

eval_code = ['sMRI.voxsize = sMRI.voxsize([',dims(1),',',dims(2),',',dims(3),']);'];
% disp(['#MRI_read_Data# ',eval_code])

fid = fopen('myDiaryFile', 'a+');
fprintf(fid, '%s\n', eval_code);
fclose(fid);
end

function [sMri] = import_mri_(iSubject, MriFile, FileFormat, isInteractive, isAutoAdjust, Comment, Labels)
% IMPORT_MRI: Import a MRI file in a subject of the Brainstorm database
% 
% USAGE: [BstMriFile, sMri] = import_mri(iSubject, MriFile, FileFormat='ALL', isInteractive=0, isAutoAdjust=1, Comment=[], Labels=[])
%               BstMriFiles = import_mri(iSubject, MriFiles, ...)   % Import multiple volumes at once
%
% INPUT:
%    - iSubject  : Indice of the subject where to import the MRI
%                  If iSubject=0 : import MRI in default subject
%    - MriFile   : Full filename of the MRI to import (format is autodetected)
%                  => if not specified : file to import is asked to the user
%    - FileFormat : String, one on the file formats in in_mri
%    - isInteractive : If 1, importation will be interactive (MRI is displayed after loading)
%    - isAutoAdjust  : If isInteractive=0 and isAutoAdjust=1, relice/resample automatically without user confirmation
%    - Comment       : Comment of the output file
%    - Labels        : Labels attached to this file (cell array {Nlabels x 3}: {index, text, RGB})
% OUTPUT:
%    - BstMriFile : Full path to the new file if success, [] if error
%%
% Initialize returned variables
BstMriFile = [];
sMri = [];
% Get Protocol information
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
	sSubject = ProtocolSubjects.DefaultSubject;
% Normal subject 
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end
%% ===== DICOM CONVERTER =====
if strcmpi(FileFormat, 'DICOM-SPM')
    % Convert DICOM to NII
    DicomFiles = MriFile;
    MriFile = in_mri_dicom_spm(DicomFiles, bst_get('BrainstormTmpDir'), isInteractive);
    if isempty(MriFile)
        return;
    end
    FileFormat = 'Nifti1';
end

%% ===== LOOP ON MULTIPLE MRI =====
if iscell(MriFile) && (length(MriFile) == 1)
    MriFile = MriFile{1};
elseif iscell(MriFile) && ~strcmpi(FileFormat, 'SPM-TPM')
    % Only allow multiple import if there is already a MRI
    if isempty(sSubject.Anatomy)
        error(['You must import the first MRI in the subject folder separately.' 10 'Please select only one volume at a time.']);
    end
    % Initialize returned values
    nFiles = length(MriFile);
    BstMriFile = cell(1, nFiles);
    sMri = cell(1, nFiles);
    % Import all volumes without supervision
    for i = 1:nFiles
        [BstMriFile{i}, sMri{i}] = import_mri(iSubject, MriFile{i}, FileFormat, isInteractive, isAutoAdjust);
    end
    % All the files are imported: exit
    return;
end

%% ===== LOAD MRI FILE =====
isProgress = bst_progress('isVisible');
if ~isProgress
    bst_progress('start', 'Import MRI', 'Loading MRI file...');
end
% MNI / Atlas?
isMni = ismember(FileFormat, {'ALL-MNI', 'ALL-MNI-ATLAS'});
isAtlas = ismember(FileFormat, {'ALL-ATLAS', 'ALL-MNI-ATLAS', 'SPM-TPM'});
% Load MRI
isNormalize = 0;
sMri = in_mri(MriFile, FileFormat, isInteractive && ~isMni, isNormalize);
if isempty(sMri)
    bst_progress('stop');
    return
end
%% ===== GET ATLAS LABELS =====
% Try to get associated labels
if isempty(Labels) && ~iscell(MriFile)
    Labels = mri_getlabels(MriFile, sMri, isAtlas);
end
% Save labels in the file structure
if ~isempty(Labels)   % Labels were found in the input folder
    sMri.Labels = Labels;
    tagAtlas = '_volatlas';
    isAtlas = 1;
elseif isAtlas    % Volume was explicitly imported as an atlas
    tagAtlas = '_volatlas';
else
    tagAtlas = '';
end
% Get atlas comment
if isAtlas && isempty(Comment) && ~iscell(MriFile)
    [fPath, fBase, fExt] = bst_fileparts(MriFile);
    switch (fBase)
        case 'aseg'
            Comment = 'ASEG';
        case 'aparc+aseg'
            Comment = 'Deskian-Killiany';
        case 'aparc.a2009s+aseg'
            Comment = 'Destrieux';
        case 'aparc.DKTatlas+aseg'
            Comment = 'DKT';
    end
end


%% ===== MANAGE MULTIPLE MRI =====
fileTag = '';
% Add new anatomy
iAnatomy = length(sSubject.Anatomy) + 1;   
% If add an extra MRI: read the first one to check that they are compatible
if (iAnatomy > 1) && (isInteractive || isAutoAdjust)
    % Load the reference MRI (the first one)
    refMriFile = sSubject.Anatomy(1).FileName;
    sMriRef = in_mri_bst(refMriFile);
    % Adding an MNI volume to an existing subject
    if isMni
        sMri = mri_reslice_mni(sMri, sMriRef, isAtlas);
        isSameSize = 1;
        errMsg = '';
    % Regular coregistration options between volumes
    else
        % If some transformation where made to the intial volume: apply them to the new one ?
        if isfield(sMriRef, 'InitTransf') && ~isempty(sMriRef.InitTransf) && any(ismember(sMriRef.InitTransf(:,1), {'permute', 'flipdim'}))
            if ~isInteractive || java_dialog('confirm', ['A transformation was applied to the reference MRI.' 10 10 'Do you want to apply the same transformation to this new volume?' 10 10], 'Import MRI')
                % Apply step by step all the transformations that have been applied to the original MRI
                for it = 1:size(sMriRef.InitTransf,1)
                    ttype = sMriRef.InitTransf{it,1};
                    val   = sMriRef.InitTransf{it,2};
                    switch (ttype)
                        case 'permute'
                            sMri.Cube = permute(sMri.Cube, [val, 4]);
                            sMri.Voxsize = sMri.Voxsize(val);
                        case 'flipdim'
                            sMri.Cube = bst_flip(sMri.Cube, val(1));
                    end
                end
                % Modifying the volume disables the option "Reslice"
                isResliceDisabled = 1;
            else
                isResliceDisabled = 0;
            end
        else
            isResliceDisabled = 0;
        end

        % === ASK REGISTRATION METHOD ===
        % Get volumes dimensions
        refSize = size(sMriRef.Cube(:,:,:,1));
        newSize = size(sMri.Cube(:,:,:,1));
        isSameSize = all(refSize == newSize) && all(round(sMriRef.Voxsize(1:3) .* 1000) == round(sMri.Voxsize(1:3) .* 1000));
        % Ask what operation to perform with this MRI
        if isInteractive
            % Initialize list of options to register this new MRI with the existing one
            strOptions = '<HTML>How to register the new volume with the reference image?<BR>';
            cellOptions = {};
            % Register with the SPM
            strOptions = [strOptions, '<BR>- <U><B>SPM</B></U>:&nbsp;&nbsp;&nbsp;Coregister the two volumes with SPM (requires SPM toolbox).'];
            cellOptions{end+1} = 'SPM';
            % Register with the MNI transformation
            strOptions = [strOptions, '<BR>- <U><B>MNI</B></U>:&nbsp;&nbsp;&nbsp;Compute the MNI transformation for both volumes (inaccurate).'];
            cellOptions{end+1} = 'MNI';
            % Skip registration
            strOptions = [strOptions, '<BR>- <U><B>Ignore</B></U>:&nbsp;&nbsp;&nbsp;The two volumes are already registered.'];
            cellOptions{end+1} = 'Ignore';
            % Ask user to make a choice
            RegMethod = java_dialog('question', [strOptions '<BR><BR></HTML>'], 'Import MRI', [], cellOptions, 'Reg+reslice');
        % In non-interactive mode: ignore if possible, or use the first option available
        else
            RegMethod = 'Ignore';
        end
        % User aborted the import
        if isempty(RegMethod)
            sMri = [];
            bst_progress('stop');
            return;
        end

        % === ASK RESLICE ===
        if isInteractive && (~strcmpi(RegMethod, 'Ignore') || ...
            (isfield(sMriRef, 'InitTransf') && ~isempty(sMriRef.InitTransf) && any(ismember(sMriRef.InitTransf(:,1), 'vox2ras')) && ...
             isfield(sMri,    'InitTransf') && ~isempty(sMri.InitTransf)    && any(ismember(sMri.InitTransf(:,1),    'vox2ras')) && ...
             ~isResliceDisabled))
            % If the volumes don't have the same size, add a warning
            if ~isSameSize
                strSizeWarn = '<BR>The two volumes have different sizes: if you answer no here, <BR>you will not be able to overlay them in the same figure.';
            else
                strSizeWarn = [];
            end
            % Ask to reslice
            isReslice = java_dialog('confirm', [...
                '<HTML><B>Reslice the volume?</B><BR><BR>' ...
                'This operation rewrites the new MRI to match the alignment, <BR>size and resolution of the original volume.' ...
                strSizeWarn ...
                '<BR><BR></HTML>'], 'Import MRI');
        % In non-interactive mode: never reslice
        else
            isReslice = 0;
        end

        % === REGISTRATION ===
        switch (RegMethod)
            case 'MNI'
                % Register the new MRI on the existing one using the MNI transformation (+ RESLICE)
                [sMri, errMsg, fileTag] = mri_coregister(sMri, sMriRef, 'mni', isReslice, isAtlas);
            case 'SPM'
                % Register the new MRI on the existing one using SPM + RESLICE
                [sMri, errMsg, fileTag] = mri_coregister(sMri, sMriRef, 'spm', isReslice, isAtlas);
            case 'Ignore'
                if isReslice
                    % Register the new MRI on the existing one using the transformation in the input files (files already registered)
                    [sMri, errMsg, fileTag] = mri_reslice(sMri, sMriRef, 'vox2ras', 'vox2ras', isAtlas);
                else
                    % Just copy the fiducials from the reference MRI
                    [sMri, errMsg, fileTag] = mri_coregister(sMri, sMriRef, 'vox2ras', isReslice, isAtlas);
                    % Transform error in warning
                    if ~isempty(errMsg) && ~isempty(sMri) && isSameSize && ~isReslice
                        disp(['BST> Warning: ' errMsg]);
                        errMsg = [];
                    end
                end
                % Copy the old SCS and NCS fields to the new file (only if registered)
                if isSameSize || isReslice
                    sMri.SCS = sMriRef.SCS;
                    %sMri.NCS = sMriRef.NCS;
                end
        end
    end
    % Stop in case of error
    if ~isempty(errMsg)
        if isInteractive
            bst_error(errMsg, [RegMethod ' MRI'], 0);
            sMri = [];
            bst_progress('stop');
            return;
        else
            error(errMsg);
        end
    end
end
%% ===== SAVE MRI IN BRAINSTORM FORMAT =====
% Add a Comment field in MRI structure, if it does not exist yet
if ~isempty(Comment)
    sMri.Comment = Comment;
    importedBaseName = file_standardize(Comment);
else
    if ~isfield(sMri, 'Comment') || isempty(sMri.Comment)
        sMri.Comment = 'MRI';
    end
    % Use filename as comment
    if (iAnatomy > 1) || isInteractive || ~isAutoAdjust
        [fPath, fBase, fExt] = bst_fileparts(MriFile);
        fBase = strrep(fBase, '.nii', '');
        if isMni
            sMri.Comment = file_unique(fBase, {sSubject.Anatomy.Comment});
        else
            sMri.Comment = file_unique([fBase, fileTag], {sSubject.Anatomy.Comment});
        end
    end
    % Add MNI tag
    if isMni
        if isfield(sMri, 'NCS') && isfield(sMri.NCS, 'y_method') && ~isempty(sMri.NCS.y_method)
            sMri.Comment = [sMri.Comment ' (MNI-' sMri.NCS.y_method ')'];
        elseif isfield(sMri, 'NCS') && isfield(sMri.NCS, 'y') && isfield(sMri.NCS, 'iy') && ~isempty(sMri.NCS.y) && ~isempty(sMri.NCS.iy)
            sMri.Comment = [sMri.Comment ' (MNI-nonlin)'];
        elseif isfield(sMri, 'NCS') && isfield(sMri.NCS, 'R') && isfield(sMri.NCS, 'T') && ~isempty(sMri.NCS.R) && ~isempty(sMri.NCS.T)
            sMri.Comment = [sMri.Comment ' (MNI-linear)'];
        else
            sMri.Comment = [sMri.Comment ' (MNI)'];
        end
    end
    % Get imported base name
    [tmp__, importedBaseName] = bst_fileparts(MriFile);
    importedBaseName = strrep(importedBaseName, 'subjectimage_', '');
    importedBaseName = strrep(importedBaseName, '_subjectimage', '');
    importedBaseName = strrep(importedBaseName, '.nii', '');
end
%% ===== MRI VIEWER =====
if isInteractive
    % First MRI: Edit fiducials
    if (iAnatomy == 1)
        % MRI Visualization and selection of fiducials (in order to align surfaces/MRI)
        hFig = view_mri(BstMriFile, 'EditMri');
        drawnow;
        bst_progress('stop');
        % Wait for the MRI Viewer to be closed
        if ishandle(hFig)
            waitfor(hFig);
        end
    % Other volumes: Display registration
    else
        % If volumes are registered
        if isSameSize || isReslice
            % Open the second volume as an overlay of the first one
            hFig = view_mri(refMriFile, BstMriFile);
            % Set the amplitude threshold to 30%
            if ~isAtlas
                panel_surface('SetDataThreshold', hFig, 1, 0.3);
            end
        else
            hFig = view_mri(BstMriFile);
        end
    end
else
    if ~isProgress
        bst_progress('stop');
    end
end
end