atlas_surferdir = 'E:\work_Code\Matlab\ZYS\DIY_App\SEEG_Mapper\Atlas\HCPex\surfer_zys\';
disp_paths(atlas_surferdir)

% souredata producted by SEEG_Process_Atlas_surfer and SEEG_Process_Atlas_surferFlat
sMNIlh = load([atlas_surferdir,'mniL_pialL_bin2.mat']);
sMNIlh_attach = load([atlas_surferdir,'subjL_sphreg.mat']);
sMNIlh_flat = load([atlas_surferdir,'mniL_pialL_bin2_flat.mat']);
sMNIlh_reg = load([atlas_surferdir,'mniL_pialL_bin2_subjLreg.mat']);

cmap_sign = jetw(255,'half');
%%
load('E:\work_Code\Matlab\Toolbox\brain_science\Registration\SphericalDemons_SDv1.5.1\zys_DIY\SphericalDemons\freesurfer\lh.DW.Atlas1to39.2.10.mat')
%%
sm_iters = [50,100,400];
mniL_raw = sMNIlh.mniL_pialL_bin2;
mniL_sm1 = Surf_smooth(mniL_raw,'iter',sm_iters(1),'type',3);
mniL_sm2 = Surf_smooth(mniL_sm1,'iter',sm_iters(2),'type',3);
mniL_sm3 = Surf_smooth(mniL_sm2,'iter',sm_iters(1),'type',3);

fig = figure('Visible',1);
next; colormapa(cmap_sign); plot_surf(mniL_sm1,'colors',sMNIlh_flat.mniLK3);
next; colormapa(cmap_sign); plot_surf(mniL_sm1,'colors',sMNIlh_flat.mniLK3);
lay_fig([12,5]); lay_fig('func','legend off; axis off'); lay_linkview(gcf,'view'); view(-90,0)
next; colormapa(cmap_sign); plot_surf(sMNIlh_flat.mniL_sph,'colors',sMNIlh_flat.mniLK3);
next; colormapa(cmap_sign); plot_surf(sMNIlh.mniL_sph,'colors',sMNIlh_flat.mniLK2);

%%
sMNIrh = load([atlas_surferdir,'mniR_pialR_bin2.mat']);
sMNIrh_attach = load([atlas_surferdir,'subjR_sphreg.mat']);
sMNIrh_flat = load([atlas_surferdir,'mniR_pialR_bin2_flat.mat']);
sMNIrh_reg = load([atlas_surferdir,'mniR_pialR_bin2_subjRreg.mat']);
%%
sm_iters = [50,100,400];
mniR_raw = sMNIrh.mniR_pialR_bin2;
mniR_sm1 = Surf_smooth(mniR_raw,'iter',sm_iters(1),'type',3);
mniR_sm2 = Surf_smooth(mniR_sm1,'iter',sm_iters(2),'type',3);
mniR_sm3 = Surf_smooth(mniR_sm2,'iter',sm_iters(1),'type',3);

fig = figure('Visible',1);
next; colormapa(cmap_sign); plot_surf(mniR_sm1,'colors',sMNIrh_flat.mniRK3);
next; colormapa(cmap_sign); plot_surf(mniR_sm1,'colors',sMNIrh_flat.mniRK3);
lay_fig([12,5]); lay_fig('func','legend off; axis off'); lay_linkview(gcf,'view'); view(90,0)
next; colormapa(cmap_sign); plot_surf(sMNIrh_flat.mniR_sph,'colors',sMNIrh_flat.mniRK3);
next; colormapa(cmap_sign); plot_surf(sMNIrh.mniR_sph,'colors',sMNIrh_flat.mniRK2);

%%
map_subjL = sMNIlh_reg.subjL_flat.Vertices(:,[1,2]);
map_subjR = sMNIrh_reg.subjR_flat.Vertices(:,[1,2]);

def_weight = [];

subjL_square_corners = 40;
area_eps = 1e-5;
square_ratio = 2;
temp = size(sMNIlh_flat.img_mniL_raw);
texNy = temp(1); texNx = temp(2);

surf_inL = sMNIlh_attach.subjL_attachL;
[surf_div,~,idx] = Surf_func_ROI(surf_inL,'divide',sMNIlh_reg.vregidxL_iswhite==0); surf_divL = surf_div{2};
vsubjL_idx = idx{2};

subjR_square_corners = 40;
% area_eps = 1e-5;
% square_ratio = 2;
% temp = size(sMNIrh_flat.img_mniR_raw);
% texNy = temp(1); texNx = temp(2);

surf_inR = sMNIrh_attach.subjR_attachR;
[surf_div,~,idx] = Surf_func_ROI(surf_inR,'divide',sMNIrh_reg.vregidxR_iswhite==0); surf_divR = surf_div{2};
vsubjR_idx = idx{2};
%%
% root_dir = 'D:\File\work_file\Dataset\Brain\openFCD\';
pialL_tail = '\surf\T1_reg.DL.left.inner.ply';
pialR_tail = '\surf\T1_reg.DL.right.inner.ply';

% root_dir = FCD_rootdir;
root_dir = 'F:\Dataset\Brain\openFCD\reg_test1\';

if_plot = 1;
if_visible = 1;

sub_dirs = dir(root_dir); sub_dirs(1:2) = [];
sub_dirs = sub_dirs([sub_dirs.isdir]);
sub_num = length(sub_dirs);

i_list = 1:sub_num;
for i_sub = 1:length(i_list)

    disp("sub "+i_sub+": "+i_list(i_sub)+" => ")

    sub_dir = sub_dirs(i_list(i_sub)).name;
    saveFile_dir = [root_dir,sub_dir,'/surf/'];

    tsubj_tic = tic;
    %%
    subjL_file = [root_dir,sub_dir,pialL_tail];
    surf_raw = Surf_read_Data(subjL_file);
    surf_subjL = surf_raw;

    [subjL_sph] = Surf_func_ROI(surf_raw,'mapping_sphere','fine',0,'get_field',0,'keep_area',0);
    % subjL_sph.Vertices = 1*mni_sph.Vertices+0*subjL_sph.Vertices;
    surf_sm1 = Surf_smooth(surf_raw,'iter',sm_iters(1),'type',3);

    subjL_curv1 = Surf_func_Val(surf_sm1,'curvature','fine',0);
    subjLK1 = mat_cut(subjL_curv1.k1,'10','90');
    subjLK1re = subjLK1;
    subjLK1re(subjLK1re<0) = 4*subjLK1re(subjLK1re<0);

    [surf_sm2] = Surf_smooth(surf_sm1,'iter',sm_iters(2),'type',3,'calc_dist',0);
    subjLK2 = get_reduce_curvK1(surf_sm2,16,subjL_sph.Vertices);
    subjLK2 = mat_cut(subjLK2,'10','90');

    [surf_sm3] = Surf_smooth(surf_sm2,'iter',sm_iters(3),'type',3,'calc_dist',0);
    subjL_sm3 = surf_sm3;
    subjLK3 = get_reduce_curvK1(surf_sm3,32,subjL_sph.Vertices);
    subjLK3 = mat_cut(subjLK3,'10','90');

    val = [subjLK3,subjLK2,subjLK1re];
    val = surf_valsmooth(surf_raw,val,50);
    subjLK3 = val(:,1); subjLK2 = val(:,2); subjLK1sm = val(:,3);
    subjLK1sm = mat_cut(subjLK1sm,'10','90');

    % temp = min((subjLK1sm-mean(subjLK1sm))./std(subjLK1sm),(subjLK3-mean(subjLK3))./std(subjLK3));
    subjLK3 = mat_cut(subjLK1sm,'10','90');

    % [subjL_sph] = Surf_func_ROI(surf_sm3,'mapping_sphere','fine',0,'get_field',0,'keep_area',0);
    subjL_sph = sMNIlh_reg.subjL_sphreg;
    subjL_sph = Surf_smooth(subjL_sph,'iter',3,'type',3);
    subjL_sph.Vertices = subjL_sph.Vertices./vecnorm(subjL_sph.Vertices,2,2)*SD_atlas.radius;

    if if_plot
        fig = figure('Visible',if_visible);
        next; colormapa(cmap_sign); plot_surf(surf_sm3,'colors',subjLK3);
        next; colormapa(cmap_sign); plot_surf(surf_sm2,'colors',subjLK2);
        next; colormapa(cmap_sign); plot_surf(surf_sm2,'colors',subjLK1); % caxis([-1,1]/2);
        lay_fig([12,5]); lay_fig('func','legend off; axis off'); lay_linkview(gcf,'view'); view(-90,0)
        next; colormapa(cmap_sign); plot_surf(subjL_sph,'colors',subjLK3);
        next; colormapa(cmap_sign); plot_surf(subjL_sph,'colors',subjLK2);
        next; colormapa(cmap_sign); plot_surf(subjL_sph,'colors',subjLK1); % caxis([-1,1]/2);

        % exportgraphics(fig,[saveFile_dir,'fig_pialL.png']); close(fig)
    end
    %%
    % V = sMNIlh_reg.subjL_sphreg.Vertices(vsubjL_idx,:); % !!!!! sphreg with sm
    % F = surf_div.Faces;
    % sph_subjLdiv = struct('Vertices',V,'Faces',F);
    % [subjL_flat,map_subjL] = Surf_func_GEO(sph_subjLdiv,'mapping_square','corners',subjL_square_corners,'area_eps',area_eps);
    %
    % map_subjL(:,1) = square_ratio*map_subjL(:,1)/max(map_subjL(:,1));
    % subjL_flat.Vertices(:,[1,2]) = map_subjL;
    subjL_flat = sMNIlh_reg.subjL_flat;

    subjLK3_div = subjLK1(vsubjL_idx); val = subjLK3_div;

    img_subjL_raw = scatteredInterpolant_img(map_subjL,val,[1,square_ratio],[texNy,texNx]);
    img_subjL = amp(img_subjL_raw);

    if if_plot
        fig = figure('Visible',if_visible); colormap jet
        next; plot_surf(surf_inL,'colors',double(sMNIlh_reg.vregidxL_iswhite));
        next; plot_surf(surf_divL,'colors',val);
        % next; plot_surf(subjL_sphreg,'colors',subjLK3);
        lay_fig([7,6]); lay_fig('func','axis off; legend off'); lay_linkview(gcf,'view'); view(90,0)
        next; plot_surf(subjL_flat,'colors',val); lighting none; axis ij
        next; showc(img_subjL)

        exportgraphics(fig,[saveFile_dir,'fig_flatL.png']); close(fig)
    end
    %%
    im_fix_cat = amp(sMNIlh_reg.img_mniL_raw);
    im_mov_cat = amp(img_subjL);

    bound_sz = size(im_mov_cat);
    pts_sz = [5,4];
    x = linspace(1,bound_sz(2),pts_sz(2)); y = linspace(1,bound_sz(1),pts_sz(1));
    movingPts = [x,y*0+bound_sz(2),fliplr(x),y*0; x*0,y,x*0+bound_sz(1),fliplr(y)];
    fixedPts = movingPts;
    reg_params = struct('Metric0Weight',1,'Metric1Weight',0.2,'NumberOfResolutions',3,'NumberOfSpatialSamples',2000);
    params = {};
    % params = {'movingPts',movingPts,'fixedPts',fixedPts,'reg_params',reg_params};

    t_tic = tic;
    [sReg,trans_params,def_field] = MRI_register(im2uint8(im_fix_cat),im2uint8(im_mov_cat),'reg_type','nonrigid', ...
        'auto_params',20,'get_field',1,'field_type','xyz',params{:});
    subjL_reg = sReg.data;

    if if_plot
        fig = figure('Visible',if_visible); lay_sub(1,3); colormap jet
        next; showc(im_fix_cat);
        next; showc(im_mov_cat);
        next; showc(subjL_reg);
        lay_fig([12,3]); lay_fig('func','axis off; legend off; axis off');
        disp_toc(t_tic)

        exportgraphics(fig,[saveFile_dir,'fig_regL.png']); close(fig)
    end
    %%
    fix_bnd_sz = 20;
    def_sz = size(def_field.data);
    if isempty(def_weight)
        def_weight = make_maskWeight(def_sz,fix_bnd_sz,2);
    end

    Vm = pts_warp(map_subjL*texNy,def_field,'field');
    Vwarp = [];
    for i = 1:3
        Vwarp(:,i) = interp2(sMNIlh_flat.imgXYZ_mniL(:,:,i),Vm(:,1),Vm(:,2));
    end

    Vraw = sMNIlh_reg.sph_subjLdiv.Vertices;
    weight = interp2(def_weight,Vm(:,1),Vm(:,2));
    weight(isnan(weight)) = 0;
    Vwarp(isnan(Vwarp)) = 0;
    Vwarp = Vwarp.*weight + Vraw.*(1-weight);
    Vwarp = Vwarp./vecnorm(Vwarp,2,2)*SD_atlas.radius;

    subjL_warp = sMNIlh_reg.subjL_sphreg; subjL_warp.Vertices(vsubjL_idx,:) = Vwarp;
    [subjL_warp] = Surf_smooth(subjL_warp,'iter',30,'type',3);
    subjL_warp.Vertices = subjL_warp.Vertices./vecnorm(subjL_warp.Vertices,2,2)*SD_atlas.radius;

    if if_plot
        fig = figure('Visible',if_visible); lay_sub(2,2); colormap jet
        next; plot_scatter(Vm,1,subjLK3_div); colormap jet; axis ij; axis tight
        next; plot_surf(sMNIlh.mniL_sph,'colors',sMNIlh_flat.mniLK3);
        next; plot_surf(sMNIlh_reg.subjL_sphreg,'colors',subjLK3);
        next; plot_surf(subjL_warp,'colors',subjLK3);
        lay_fig([7,6]); lay_fig('func','axis off; legend off'); lay_linkview(gcf,'view');

        exportgraphics(fig,[saveFile_dir,'fig_warpL.png']); close(fig)
    end

    %%
    save([saveFile_dir,'subjL_warp.mat'],'subjL_warp','subjLK3','subjLK2','subjLK1','def_field','img_subjL')
    %%
    subjR_file = [root_dir,sub_dir,pialR_tail];
    surf_raw = Surf_read_Data(subjR_file);
    surf_subjR = surf_raw;

    [subjR_sph] = Surf_func_ROI(surf_raw,'mapping_sphere','fine',0,'get_field',0,'keep_area',0);
    % subjR_sph.Vertices = 1*mni_sph.Vertices+0*subjR_sph.Vertices;
    surf_sm1 = Surf_smooth(surf_raw,'iter',sm_iters(1),'type',3);

    subjR_curv1 = Surf_func_Val(surf_sm1,'curvature','fine',0);
    subjRK1 = mat_cut(subjR_curv1.k1,'10','90');
    subjRK1re = subjRK1;
    subjRK1re(subjRK1re<0) = 4*subjRK1re(subjRK1re<0);

    [surf_sm2] = Surf_smooth(surf_sm1,'iter',sm_iters(2),'type',3,'calc_dist',0);
    subjRK2 = get_reduce_curvK1(surf_sm2,16,subjR_sph.Vertices);
    subjRK2 = mat_cut(subjRK2,'10','90');

    [surf_sm3] = Surf_smooth(surf_sm2,'iter',sm_iters(3),'type',3,'calc_dist',0);
    subjR_sm3 = surf_sm3;
    subjRK3 = get_reduce_curvK1(surf_sm3,32,subjR_sph.Vertices);
    subjRK3 = mat_cut(subjRK3,'10','90');

    val = [subjRK3,subjRK2,subjRK1re];
    val = surf_valsmooth(surf_raw,val,50);
    subjRK3 = val(:,1); subjRK2 = val(:,2); subjRK1sm = val(:,3);
    subjRK1sm = mat_cut(subjRK1sm,'10','90');

    % temp = min((subjRK1sm-mean(subjRK1sm))./std(subjRK1sm),(subjRK3-mean(subjRK3))./std(subjRK3));
    subjRK3 = mat_cut(subjRK1sm,'10','90');

    % [subjR_sph] = Surf_func_ROI(surf_sm3,'mapping_sphere','fine',0,'get_field',0,'keep_area',0);
    subjR_sph = sMNIrh_reg.subjR_sphreg;
    subjR_sph = Surf_smooth(subjR_sph,'iter',3,'type',3);
    subjR_sph.Vertices = subjR_sph.Vertices./vecnorm(subjR_sph.Vertices,2,2)*SD_atlas.radius;

    if if_plot
        fig = figure('Visible',if_visible);
        next; colormapa(cmap_sign); plot_surf(surf_sm3,'colors',subjRK3);
        next; colormapa(cmap_sign); plot_surf(surf_sm2,'colors',subjRK2);
        next; colormapa(cmap_sign); plot_surf(surf_sm2,'colors',subjRK1); % caxis([-1,1]/2);
        lay_fig([12,5]); lay_fig('func','legend off; axis off'); lay_linkview(gcf,'view'); view(90,0)
        next; colormapa(cmap_sign); plot_surf(subjR_sph,'colors',subjRK3);
        next; colormapa(cmap_sign); plot_surf(subjR_sph,'colors',subjRK2);
        next; colormapa(cmap_sign); plot_surf(subjR_sph,'colors',subjRK1); % caxis([-1,1]/2);

        exportgraphics(fig,[saveFile_dir,'fig_pialR.png']); close(fig)
    end
    %%
    % V = sMNIlh_reg.subjR_sphreg.Vertices(vsubjR_idx,:); % !!!!! sphreg with sm
    % F = surf_div.Faces;
    % sph_subjRdiv = struct('Vertices',V,'Faces',F);
    % [subjR_flat,map_subjR] = Surf_func_GEO(sph_subjRdiv,'mapping_square','corners',subjR_square_corners,'area_eps',area_eps);
    %
    % map_subjR(:,1) = square_ratio*map_subjR(:,1)/max(map_subjR(:,1));
    % subjR_flat.Vertices(:,[1,2]) = map_subjR;
    subjR_flat = sMNIrh_reg.subjR_flat;

    subjRK3_div = subjRK1(vsubjR_idx); val = subjRK3_div;

    img_subjR_raw = scatteredInterpolant_img(map_subjR,val,[1,square_ratio],[texNy,texNx]);
    img_subjR = amp(img_subjR_raw);

    if if_plot
        fig = figure('Visible',if_visible); colormap jet
        next; plot_surf(surf_inR,'colors',double(sMNIrh_reg.vregidxR_iswhite));
        next; plot_surf(surf_divR,'colors',val);
        % next; plot_surf(subjR_sphreg,'colors',subjRK3);
        lay_fig([7,6]); lay_fig('func','axis off; legend off'); lay_linkview(gcf,'view'); view(-90,0)
        next; plot_surf(subjR_flat,'colors',val); lighting none; axis ij
        next; showc(img_subjR)

        exportgraphics(fig,[saveFile_dir,'fig_flatR.png']); close(fig)
    end
    %%
    im_fix_cat = amp(sMNIrh_reg.img_mniR_raw);
    im_mov_cat = amp(img_subjR);

    bound_sz = size(im_mov_cat);
    pts_sz = [5,4];
    x = linspace(1,bound_sz(2),pts_sz(2)); y = linspace(1,bound_sz(1),pts_sz(1));
    movingPts = [x,y*0+bound_sz(2),fliplr(x),y*0; x*0,y,x*0+bound_sz(1),fliplr(y)];
    fixedPts = movingPts;
    reg_params = struct('Metric0Weight',1,'Metric1Weight',0.2,'NumberOfResolutions',3,'NumberOfSpatialSamples',2000);
    params = {};
    % params = {'movingPts',movingPts,'fixedPts',fixedPts,'reg_params',reg_params};

    t_tic = tic;
    [sReg,trans_params,def_field] = MRI_register(im2uint8(im_fix_cat),im2uint8(im_mov_cat),'reg_type','nonrigid', ...
        'auto_params',20,'get_field',1,'field_type','xyz',params{:});
    subjR_reg = sReg.data;

    if if_plot
        fig = figure('Visible',if_visible); lay_sub(1,3); colormap jet
        next; showc(im_fix_cat);
        next; showc(im_mov_cat);
        next; showc(subjR_reg);
        lay_fig([12,3]); lay_fig('func','legend off; axis off');
        disp_toc(t_tic)

        exportgraphics(fig,[saveFile_dir,'fig_regR.png']); close(fig)
    end
    %%
    fix_bnd_sz = 20;
    def_sz = size(def_field.data);
    if isempty(def_weight)
        def_weight = make_maskWeight(def_sz,fix_bnd_sz,2);
    end

    Vm = pts_warp(map_subjR*texNy,def_field,'field');
    Vwarp = [];
    for i = 1:3
        Vwarp(:,i) = interp2(sMNIrh_flat.imgXYZ_mniR(:,:,i),Vm(:,1),Vm(:,2));
    end

    Vraw = sMNIrh_reg.sph_subjRdiv.Vertices;
    weight = interp2(def_weight,Vm(:,1),Vm(:,2));
    weight(isnan(weight)) = 0;
    Vwarp(isnan(Vwarp)) = 0;
    Vwarp = Vwarp.*weight + Vraw.*(1-weight);
    Vwarp = Vwarp./vecnorm(Vwarp,2,2)*SD_atlas.radius;

    subjR_warp = sMNIrh_reg.subjR_sphreg; subjR_warp.Vertices(vsubjR_idx,:) = Vwarp;
    [subjR_warp] = Surf_smooth(subjR_warp,'iter',30,'type',3);
    subjR_warp.Vertices = subjR_warp.Vertices./vecnorm(subjR_warp.Vertices,2,2)*SD_atlas.radius;

    if if_plot
        fig = figure('Visible',if_visible); lay_sub(2,2); colormap jet
        next; plot_scatter(Vm,1,subjRK3_div); colormap jet; axis ij; axis tight
        next; plot_surf(sMNIrh.mniR_sph,'colors',sMNIrh_flat.mniRK3);
        next; plot_surf(sMNIrh_reg.subjR_sphreg,'colors',subjRK3);
        next; plot_surf(subjR_warp,'colors',subjRK3);
        lay_fig([7,6]); lay_fig('func','axis off; legend off'); lay_linkview(gcf,'view');

        exportgraphics(fig,[saveFile_dir,'fig_warpR.png']); close(fig)
    end
    %%
    save([saveFile_dir,'subjR_warp.mat'],'subjR_warp','subjRK3','subjRK2','subjRK1','def_field','img_subjR')
    %%
    disp_toc(tsubj_tic)

    if mod(i_sub,5)
        clc
        % close all hidden
    end
end

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

function [K1] = get_reduce_curvK1(surf_raw,bin,sph_v)
[surf_rd1, idx_rd] = Surf_func_ROI(surf_raw,'reduce','bin',bin,'get_idx',1);
surf_rd2 = surf_rd1; surf_rd2.Vertices = surf_raw.Vertices(idx_rd,:);
val_curv = Surf_func_Val(surf_rd2,'curvature','fine',0);
v_q = sph_v; v = v_q(idx_rd,:);
% K1 = val_curv.k1; K2 = val_curv.k2;
% K = interp_sph(v,[K1,K2],v_q);
% K1 = K(:,1); K2 = K(:,2);
K1 = interp_sph(v,val_curv.k1,v_q);
end

function val = surf_valsmooth(surf_in,val,sm_iter)
temp = surf_in; temp.Vertices(:,1:size(val,2)) = val;
temp = Surf_smooth(temp,'iter',sm_iter,'type',3);
val = temp.Vertices(:,1:size(val,2));
end

function eigenmodes = get_surfEigmode(surf_raw,num_modes)
PWD = pwd;
brainEigmode_dir = 'E:\work_Code\Matlab\from_github\BrainEigenmodes-main\';
cd(brainEigmode_dir)
data_dir = 'data\test_data_zys/';
surface_interest = 'zys';
hemisphere = 'lh';
mesh_interest = 'white';

% num_modes = 10;

surface_input_filename = sprintf('%s/%s_%s-%s.vtk', data_dir,surface_interest, mesh_interest, hemisphere);
% mask_filename = sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere);
output_eval_filename = sprintf('%s/%s_%s-%s_eval_%i.txt',data_dir,surface_interest, mesh_interest, hemisphere, num_modes);
output_emode_filename = sprintf('%s/%s_%s-%s_emode_%i.txt',data_dir,surface_interest, mesh_interest, hemisphere, num_modes);

Surf_save_File(surf_raw,surface_input_filename);

tic
cmd_str = [];
cmd_str{1} = 'CALL conda.bat activate pytorch3d';
cmd_str{end+1,1} = ['python surface_eigenmodes.py ',surface_input_filename,' ',output_eval_filename,' ',output_emode_filename,' -save_cut 0 -N ',num2str(num_modes),' -is_mask 0'];
bat_name = 'test1.bat';
writecell(cmd_str,bat_name,'FileType','text')
system(bat_name);
toc

eigenmodes = dlmread(output_emode_filename);
cd(PWD)
end

function [disk_img,Vnew] = sph2disk(cur_sph,val_disk,texN,if_plot)
surf_disk = cur_sph;
v = surf_disk.Vertices;
[TH,PHI] = cart2sph(v(:,3),v(:,2),v(:,1));
PHI = PHI+pi/2;
Vnew = [cos(TH).*PHI,sin(TH).*PHI];

x = linspace(-pi,pi,texN);
[Y,X] = ndgrid(x);

interp_h = scatteredInterpolant(Vnew,val_disk);
val_new = interp_h([X(:),Y(:)]);
disk_img = reshape(val_new,texN,[]);

if if_plot
    Vplot = Vnew;
    Vplot(PHI>pi*0.95,:) = NaN; Vplot(:,3) = 0;
    surf_disk.Vertices = Vplot;

    fig = figure('Visible',if_visible);;
    next; plot_surf(surf_disk,'colors',val_disk); next; showc(disk_img); axis xy
    lay_fig([12,3]); colormap jet
end
end

function img = scatteredInterpolant_img(v,val,lim_max,img_sz)
x = linspace(0,lim_max(2),img_sz(2)); y = linspace(0,lim_max(1),img_sz(1));
[X,Y] = meshgrid(x,y);
interp_fmni = scatteredInterpolant(v,val);
img = interp_fmni([X(:),Y(:)]);
img = reshape(img,img_sz);
end

function def_weight = make_maskWeight(def_sz,fix_bnd_sz,gamma)
def_weight = ones(round(def_sz([1,2])/4));
def_weight = imfilter(def_weight,fspecial('average',fix_bnd_sz));
def_weight = imresize(def_weight,def_sz([1,2]));
def_weight = amp(def_weight,0).^gamma;
end