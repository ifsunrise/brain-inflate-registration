
global if_show
if_show = 0;
ext_sm = 6;
if_demons = 1;

atlas_surferdir = 'E:\work_Code\Matlab\ZYS\DIY_App\SEEG_Mapper\Atlas\HCPex\surfer_zys\';
disp_paths(atlas_surferdir)

mniF_mask = MRI_read_Data([atlas_surferdir,'mniF_cortex_mask.nii.gz']);
mniF_mask = mniF_mask.data>0;
mni_inflat = MRI_read_Data([atlas_surferdir,'inflat_atlas\iter0\mniF_T1.mean.nii.gz']);
mni_inflat = mni_inflat.data;

sMNIlh = load([atlas_surferdir,'mniL_pialL_bin2.mat']);
sMNIrh = load([atlas_surferdir,'mniR_pialR_bin2.mat']);
sMNIlh_flat = load([atlas_surferdir,'mniL_pialL_bin2_flat.mat']);
sMNIrh_flat = load([atlas_surferdir,'mniR_pialR_bin2_flat.mat']);
load([atlas_surferdir,'mni_frac_inflat.mat'],'mniL_inflat','mniR_inflat')

surf_mniR = Surf_read_Data([atlas_surferdir,'mniR_pialR_bin2.~']);
surf_mniL = Surf_read_Data([atlas_surferdir,'mniL_pialL_bin2.~']);
cmap_FSLre = surf_mniR.cmap_FSLre;

mnisz = size(mni_inflat);
[X,Y,Z] = meshgrid(1:mnisz(2),1:mnisz(1),1:mnisz(3));
Vall = [X(:),Y(:),Z(:)];

% root_dir = 'D:\File\work_file\Dataset\Brain\openFCD\';
root_dir = FCD_rootdir;

sub_dirs = dir(root_dir); sub_dirs(1:2) = [];
sub_dirs = sub_dirs([sub_dirs.isdir]);
sub_num = length(sub_dirs);

i_list = 1:sub_num;
%%

for i_sub = 1:length(i_list)
    disp("------------ sub "+i_sub+": "+i_list(i_sub)+" ------------")
    sub_dir = sub_dirs(i_list(i_sub)).name;

    if ~if_demons
        surf_dir = [root_dir,sub_dir,'\surf\'];
    else
        surf_dir = [root_dir,sub_dir,'\surf2\'];
    end

    tic
    disp('reading, improving and warping inner surf...')

    if ~isfile([surf_dir,'/subjL_warp.mat'])
        continue
    end

    subjLreg = load([surf_dir,'/subjL_warp.mat']);
    subjRreg = load([surf_dir,'/subjR_warp.mat']);
    mri_subj = MRI_read_Data([surf_dir,'\..\T1_rb\T1.reg.norm.nii.gz'],'verbose',0);
    mri_subj.data = mat_cut(3*mri_subj.data-1,0);
    % mask_subj = MRI_read_Data([surf_dir,'\..\T1_rb\T1.reg.norm.nii.gz'],'verbose',0);

    subjL_warp = subjLreg.subjL_warp;
    subjR_warp = subjRreg.subjR_warp;

    if if_show
        figure('Name','sphere reg'); lay_sub(1,4); colormap jet
        next; plot_surf(sMNIlh.mniL_sph,'colors',sMNIlh_flat.mniLK3);
        next; plot_surf(subjL_warp,'colors',subjLreg.subjLK3);
        next; plot_surf(sMNIrh.mniR_sph,'colors',sMNIrh_flat.mniRK3);
        next; plot_surf(subjR_warp,'colors',subjRreg.subjRK3);
        lay_fig([12,3]); lay_fig('func','legend off'); lay_linkview(gcf,'view');
    end
    %% attach surface to gray-white matter border
    subjL_file = [root_dir,sub_dir,'\surf\T1_reg.DL.L.inner.ply'];
    surf_subjL = Surf_read_Data(subjL_file);
    subjR_file = [root_dir,sub_dir,'\surf\T1_reg.DL.R.inner.ply'];
    surf_subjR = Surf_read_Data(subjR_file);

    % subjL_file = [surf_dir,'\T1_reg.DL.right.inner.ply'];
    % surf_subjL_raw = Surf_read_Data(subjL_file);
    % subjR_file = [surf_dir,'\T1_reg.DL.left.inner.ply'];
    % surf_subjR_raw = Surf_read_Data(subjR_file);

    % v_shift = [2,2,2];
    % surf_subjL_raw.Vertices = surf_subjL_raw.Vertices + v_shift;
    % surf_subjR_raw.Vertices = surf_subjR_raw.Vertices + v_shift;
    % 
    % if if_show
    %     figure('Name','attached surf');
    %     next; slice_mri(mri_subj,[50,80,100],'add_surf',{surf_subjL_raw,surf_subjR_raw}); view(150,30); legend off
    % end
    % 
    % [surf_subjL,surf_subjR] = surfs_attach_grad(mri_subj,surf_subjL_raw,surf_subjR_raw);
    %%
    surfL_warp = surf_subjL;
    surfL_warp.Vertices = interp_sph(sMNIlh.mniL_sph.Vertices,double(mniL_inflat.Vertices),subjL_warp.Vertices);

    surfR_warp = surf_subjR;
    surfR_warp.Vertices = interp_sph(sMNIrh.mniR_sph.Vertices,double(mniR_inflat.Vertices),subjR_warp.Vertices);

    if if_show
        figure('Name','surf warp');
        next; colormapa(cmap_FSLre); plot_surf(surf_mniL,'colors',surf_mniL.Labels_FSL); view(-90,0)
        next; plot_surf(surfL_warp); view(-90,0)
        next; colormapa(cmap_FSLre); plot_surf(surf_mniR,'colors',surf_mniR.Labels_FSL); view(90,0)
        next; plot_surf(surfR_warp); view(90,0)
        lay_fig([12,3]); lay_fig('func','legend off; axis off');
    end
    %%
    disp("using "+round(toc,1)+"s")
    disp('generating multi-layers pials...')
    tic

    pl_d = [8,8];
    [surf_subjL_pl,surf_subjR_pl,Vm,Vf,sksubj,sksubj_dist,sksubj_warp,sksubj_warp_indist] = surfs_multiPeel(mri_subj,mni_inflat,surf_subjL,surf_subjR,surfL_warp,surfR_warp);
    sksubj_pl_fill = Surf_voxelise([surf_subjL_pl,surf_subjR_pl],'mode',1,'ref',mri_subj);
    %%
    disp("using "+round(toc,1)+"s")
    disp('warping volume (mni_flat => subj) via interpolation...')
    tic

    opts = {'linear','none'};
    interp_h = scatteredInterpolant(Vm,Vf,opts{:});

    sk_mask_sm = imgaussfilt3(single(sksubj.data),ext_sm) + (sksubj_dist>2)>0;
    % sum(sk_mask_sm(:))

    if if_show
        mask = mniF_mask;
        Vout = round(interp_h(Vall(mask>0,:)));
        idx = sub2ind(mnisz,Vout(:,2),Vout(:,1),Vout(:,3));
        mask_new = mni_inflat*0;
        mask_new(idx) = 1;

        sksubj = Surf_voxelise({surf_subjL,surf_subjR},'mode',0,'ref',mri_subj);
        mask = sksubj.data;

        figure('Name','surf inflatReg');
        next; plot_mri({mask+mri_subj.data,mask_new+mni_inflat},'off','enhance',10)
        lay_fig([12,4])
    end
    %%
    Vm_ext = Vall(sk_mask_sm>0,:);
    Vout_ext = interp_h(Vm_ext);
    idx_ext = sk_mask_sm>0;
    %%

    mask_warp = single(sksubj_pl_fill.data>0);
    Vout_extsm = Vout_ext;
    Vout_extsm(isnan(Vout_extsm)) = 0;

    sm_sz = 3; sm_std = 0.65;
    if sm_sz>0
        [XYZ_mni2subj,V_extsm,mask_warpsm] = smooth_XYZ(mnisz,Vout_extsm,idx_ext,mask_warp,sm_sz,sm_std);
    else
        Vall_new = Vall*0; Vall_new(idx_ext,:) = Vout_extsm;
        XYZ_mni2subj = reshape(Vall_new,[mnisz,3]);
    end
    XYZ_mni2subj = single(XYZ_mni2subj);

    % mri_warpsm = getVol_vals(mri_inflat,Vout_extsm,idx_ext);
    %
    % mri_warpsm.data = mri_warpsm.data.*imgaussfilt3(single(mask_warp),0.75);
    % mri_warpsm.data(mri_warpsm.data<0 | isnan(mri_warpsm.data)) = 0;
    % mri_warpsm.data(sksubj_dist>pl_d(1)) = 2;
    %
    % if if_show
    %     mri_warp = getVol_vals(mri_inflat,Vout_ext,idx_ext);
    %     mri_warp.data(mri_warp.data<0) = 0;
    %
    %     figure; lay_sub(1,4)
    %     next; plot_mri({mask_warp,mri_warp.data,mri_warpsm,mri_subj},'off','enhance',10);
    %     lay_fig([12,3]); lay_linkview
    % end

    %%
    disp("using "+round(toc,1)+"s")
    disp('inverse warping volume (subj => mni_inflat) via interpolation...')
    tic

    opts = {'linear','linear'};
    interpinv_h = scatteredInterpolant(Vf,Vm,opts{:});

    maskinv = sksubj_warp.data;
    sk_maskinv_sm = imgaussfilt3(single(maskinv),ext_sm*1) + (sksubj_warp_indist.data>4) >0;
    % sum(sk_maskinv_sm(:));

    Vm_extinv = Vall(sk_maskinv_sm>0,:);
    Vout_extinv = interpinv_h(Vm_extinv);
    idx_extinv = sk_maskinv_sm>0;

    Vall_new = Vall*0; Vall_new(idx_extinv,:) = Vout_extinv;
    XYZ_subj2mni = reshape(Vall_new,[mnisz,3]);
    XYZ_subj2mni = single(XYZ_subj2mni.*mniF_mask);

    % mri_warpinv = getVol_vals(mri_subj,Vout_extinv,idx_extinv);
    % mri_warpinv.data(mri_warpinv.data<0 | isnan(mri_warpinv.data)) = 0;
    % mri_warpinv.data = mri_warpinv.data.*mniF_mask;

    disp("using "+round(toc,1)+"s")
    %%
    disp('saving results...')
    tic

    if ~if_demons
        saveFile_dir = [root_dir,sub_dir,'\save\inflat\'];
    else
        saveFile_dir = [root_dir,sub_dir,'\save\inflat2\'];
    end
    mkdir(saveFile_dir)

    mri_warpsm = MRI_warp(XYZ_mni2subj,mni_inflat,'BST');
    MRI_save_File(mri_warpsm.data,[saveFile_dir,'mni_inflatReg.nii.gz'],'info',mri_subj)

    figure('Visible',if_show); plot_mri(mri_warpsm,'off','enhance',10)
    exportgraphics(gcf,[saveFile_dir,'mni_inflatReg.png']); close(gcf)
    
    % MRI_save_File(mri_warpinv.data,[saveFile_dir,'T1.norm.reg_inflat.nii.gz'])

    MRI_save_File(XYZ_mni2subj,[saveFile_dir,'XYZ_mniF2subj.gz'],'info',mri_subj)
    MRI_save_File(XYZ_subj2mni,[saveFile_dir,'XYZ_subj2mniF.nii.gz'])
    disp("using "+round(toc,1)+"s")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function val_new = interp_sph(v,val_raw,v_q)
[TH,PHI] = cart2sph(v(:,1),v(:,2),v(:,3));
[TH_q,PHI_q] = cart2sph(v_q(:,1),v_q(:,2),v_q(:,3));

val_new = [];
for i = 1:2
    rot_TH = pi*(i-1);
    th = mod(TH+rot_TH,2*pi);
    th_q = mod(TH_q+rot_TH,2*pi);

    vv_in = [th,PHI];
    vv_q = [th_q,PHI_q];

    for j = 1:size(val_raw,2)
        interp_h = scatteredInterpolant(vv_in,val_raw(:,j));
        val_new(:,j,i) = interp_h(vv_q);
    end
end
w = 1 - abs(mod(TH_q,2*pi)/pi - 1);
w = permute(w,[1,3,2]);
val_new = val_new(:,:,1).*w + val_new(:,:,2).*(1-w);
end

function [XYZ_new,Vout_new,mask_warpsm] = smooth_XYZ(mnisz,Vout_extsm,idx_ext,mask_warp,sm_sz,sm_std)
Vall_new = zeros([prod(mnisz),3]);
Vall_new(idx_ext,:) = Vout_extsm;
XYZ_new = reshape(Vall_new,[mnisz,3]);
mask_warpsm = smooth3(mask_warp,"gaussian",sm_sz);
for i = 1:3
    XYZ_new(:,:,:,i) = smooth3(XYZ_new(:,:,:,i).*mask_warp,"gaussian",sm_sz,sm_std)./(1e-6+mask_warpsm);
end
Vall_new = reshape(XYZ_new,[],3);
Vout_new = Vall_new(idx_ext,:);
end

function mri_warp = getVol_vals(mri_src,V,idx)
if isnumeric(mri_src)
    mri_src = single(mri_src);
end
mri_src = MRI_read_Data(mri_src);
temp = MRI_func_ROI(mri_src,'sample',V,'method',1);
mri_warp = mri_src; mri_warp.data(:) = 0;
mri_warp.data(idx) = temp;
end

function [surf_subjL,surf_subjR] = surfs_attach_grad(mri_subj,surf_subjL_raw,surf_subjR_raw)
LAB = mri_subj.data;
WG = LAB;
WG(WG<1) = 1;
WG = imgradient3(WG);

opts ={'attach_grad','map',WG,'d',1,'iter',12,'d_decay',1.25};
[surf_subjL,Gmean,d] = Surf_func_ROI(surf_subjL_raw,opts{:});
% Gmean(round(linspace(1,end,10)))
[surf_subjR,Gmean,d] = Surf_func_ROI(surf_subjR_raw,opts{:});
% Gmean(round(linspace(1,end,10)))

global if_show
if if_show
    Gmean(round(linspace(1,end,10)))
    sksubj_raw = Surf_voxelise({surf_subjL_raw,surf_subjR_raw},'mode',0,'ref',mri_subj);
    sksubj = Surf_voxelise({surf_subjL,surf_subjR},'mode',0,'ref',mri_subj);

    figure('Name','attach grad');
    % next; plot_mri(WG,'off');
    next; slice_mri(WG,[60,90,110],'add_surf',{surf_subjL_raw,surf_subjR_raw}); view(90,0)
    next; slice_mri(WG,[60,90,110],'add_surf',{surf_subjL,surf_subjR}); view(90,0)
    lay_linkview;
    next; plot_mri({mri_subj.data+sksubj_raw.data, mri_subj.data+sksubj.data})
    lay_fig([12,3]); lay_fig('func','legend off; axis off');
end
end

function [surf_subjL_pl,surf_subjR_pl,Vm,Vf,sksubj,sksubj_dist,sksubj_warp,sksubj_warp_indist] = surfs_multiPeel(mri_subj,mni_inflat,surf_subjL,surf_subjR,surfL_warp,surfR_warp)

pl_d = 8*[1,1];
pl_iter = 8;
pl_d_sm = pl_d*0.75;

%%
sksubj = Surf_voxelise({surf_subjL,surf_subjR},'mode',0,'ref',mri_subj);
sksubj_fill = Surf_voxelise({surf_subjL,surf_subjR},'mode',1,'ref',mri_subj);
sksubj_indist = MRI_func_ROI(sksubj_fill.data==0,'dist_map');
sksubj_outdist = MRI_func_ROI(sksubj_fill.data==0,'dist_map','if_inv',1,'ext',-1);
sksubj_dist = sksubj_indist.data - sksubj_outdist.data;

global if_show
if if_show
    figure('Name','skull distance'); lay_sub(1,4); next; plot_mri({mri_subj,mri_subj.data+sksubj.data,sksubj_fill,amp(sksubj_dist)},'off');
    lay_fig([12,3]);
end
%%
sksubj_warp = Surf_voxelise({surfL_warp,surfR_warp},'mode',0,'ref',mri_subj);
sksubj_warp_fill = Surf_voxelise({surfL_warp,surfR_warp},'mode',1,'ref',mri_subj);
sksubj_warp_indist = MRI_func_ROI(sksubj_warp_fill.data==0,'dist_map');
sksubj_warp_outdist = MRI_func_ROI(sksubj_warp_fill.data==0,'dist_map','if_inv',1,'ext',-1');

if if_show
    mri = mni_inflat;
    sksubj_warp_dist = sksubj_warp_indist.data - sksubj_warp_outdist.data;

    figure('Name','warp distance'); lay_sub(1,4); next; plot_mri({mri,sksubj_warp.data+mri,sksubj_warp_fill,amp(sksubj_warp_dist)},'off','enhance',10);
    lay_fig([12,3]);
end
%%

distmap_subj = {sksubj_outdist,sksubj_indist};
distmap_subj_warp = {sksubj_warp_outdist,sksubj_warp_indist};

surf_subjL_pl = []; surf_subjL_plwarp = [];
for i = 1:2
    [surf_subjL_pl{i},~,vL_pl_iters{i}] = Surf_func_ROI(surf_subjL,'peel_distgrad','map',distmap_subj{i},'d',pl_d(i),'iter',pl_iter);
    d_real = squeeze(vecnorm(diff(cat(3,surf_subjL.Vertices,vL_pl_iters{i}),1,3),2,2));
    % mean(d_real(:))
    [surf_subjL_plwarp{i},~,vL_plwarp_iters{i}] = Surf_func_ROI(surfL_warp,'peel_distgrad','map',distmap_subj_warp{i},'d',pl_d_sm(i),'iter',pl_iter);
end
surf_subjR_pl = []; surf_subjR_plwarp = [];
for i = 1:2
    [surf_subjR_pl{i},~,vR_pl_iters{i}] = Surf_func_ROI(surf_subjR,'peel_distgrad','map',distmap_subj{i},'d',pl_d(i),'iter',pl_iter);
    d_real = squeeze(vecnorm(diff(cat(3,surf_subjR.Vertices,vR_pl_iters{i}),1,3),2,2));
    % mean(d_real(:))
    [surf_subjR_plwarp{i},~,vR_plwarp_iters{i}] = Surf_func_ROI(surfR_warp,'peel_distgrad','map',distmap_subj_warp{i},'d',pl_d_sm(i),'iter',pl_iter);
end

%%
Vm = [vL_pl_iters{1}; vL_pl_iters{2}; vR_pl_iters{1}; vR_pl_iters{2}];
Vf = [vL_plwarp_iters{1}; vL_plwarp_iters{2}; vR_plwarp_iters{1}; vR_plwarp_iters{2}];
Vm = reshape(permute(Vm,[1,3,2]),[],3);
Vf = reshape(permute(Vf,[1,3,2]),[],3);

Vm = [Vm; surf_subjL.Vertices; surf_subjR.Vertices];
Vf = [Vf; surfL_warp.Vertices; surfR_warp.Vertices];
Vm = double(Vm); Vf = double(Vf);

if if_show
    figure('Name','peel warp');
    next; plot_surf(surf_subjL_pl,[1,1,0; 0,0,1],0.3); view(-90,0)
    next; plot_surf(surf_subjL_plwarp,[1,1,0; 0,0,1],0.3); view(-90,0)
    next; plot_surf(surf_subjR_pl,[1,1,0; 0,0,1],0.3); view(90,0)
    next; plot_surf(surf_subjR_plwarp,[1,1,0; 0,0,1],0.3); view(90,0)
    lay_fig([12,3]); lay_fig('func','legend off; axis off');
end

if if_show
    sksubj_pl = Surf_voxelise([surf_subjL_pl,surf_subjR_pl],'mode',0,'ref',mri_subj);
    sksubj_plwarp = Surf_voxelise([surf_subjL_plwarp,surf_subjR_plwarp],'mode',0,'ref',mri_subj);
    figure('Name','peel voxelize'); lay_sub(1,4); next; plot_mri({mri_subj.data+sksubj.data,mri_subj.data+ismember(sksubj_pl.data,[1,3]), ...
        mri_subj.data+ismember(sksubj_pl.data,[2,4]),mri+(sksubj_plwarp.data>0)},'off','enhance',10);
    lay_fig([12,3]); lay_linkview
end
end