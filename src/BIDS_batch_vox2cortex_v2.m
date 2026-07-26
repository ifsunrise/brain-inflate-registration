% FCD_rootdir = FCD_rootdir;

mri_tag = 'T1_reg';
voxin_sz = [218 182 182];

vox2cortex_root = 'E:\work_Code\Python\MiniProjects\from_github\Vox2Cortex-main\';
vox2cortex_dataroot = [vox2cortex_root,'dataset_root\dataset_raw\'];
vox2cortex_meshroot = [vox2cortex_root,'experiments\public_model\test_template_42016_DATASET_NAME\meshes'];

sub_dirs = dir(FCD_rootdir); sub_dirs(1:2) = [];
sub_dirs = sub_dirs([sub_dirs.isdir]);
sub_num = length(sub_dirs);

i_list = 1:sub_num;
%%

for i_sub = 1:length(i_list)
    tic
    % disp_wait(i_sub,length(i_list))
    fprintf("sub "+i_sub+": "+i_list(i_sub)+" => ")
    sub_dir = sub_dirs(i_list(i_sub)).name;

    T1rb_dir = [FCD_rootdir,sub_dir,'\T1_rb\'];
    T1 = MRI_read_Data([T1rb_dir,mri_tag,'.nii.gz'],'verbose',0);
    mri = T1.data;

    mri_rsz = imresize3(mri,voxin_sz);
    mri_rsz = MRI_read_Data(mri_rsz,'voxsize',[1,1,1],'verbose',0); mri_rsz.T_mat = [];

    save_dir = [vox2cortex_dataroot,'\',sub_dir,'\'];
    mkdir_silent(save_dir)
    MRI_save_File(mri_rsz,[save_dir,'mri.nii.gz']);
    t_use = toc;
    fprintf(round(t_use,1)+"s")
    fprintf('\n')
end

%%
sub_names = {sub_dirs(i_list).name}';
writecell(sub_names,[vox2cortex_dataroot,'test_ids.txt'])

winopen(vox2cortex_root)
% system([vox2cortex_root,'run.bat'])
%%

mesh_tags = {'_epoch76_class1_voxelpred','_epoch76_struc0_meshpred','_epoch76_struc1_meshpred','_epoch76_struc2_meshpred','_epoch76_struc3_meshpred'};
% mesh_names = {'bndpial','left.inner','right.inner','left.pial','right.pial'};
ply_names = {'bndpial','right.inner','left.inner','right.pial','left.pial'};
plyat_names = {'R.inner','L.inner','R.pial','L.pial'};
v_shift = [2,2,2];
attach_opts = {'d',0.75,'iter',12,'d_decay',1.25,'sm_iter',1};
if_show = 1;
% plot_pos = [60,90,110];
plot_pos = {[],[],100};

if if_show
    fig = figure('Name','attach grad');
end
for i_sub = 1:length(i_list)
    tic
    % disp_wait(i_sub,length(i_list))
    fprintf("sub "+i_sub+": "+i_list(i_sub)+" => ")
    sub_dir = sub_dirs(i_list(i_sub)).name;

    surf_dir = [FCD_rootdir,sub_dir,'\surf\'];
    mkdir_silent(surf_dir)

    T1rb_dir = [FCD_rootdir,sub_dir,'\T1_rb\'];

    for i_mesh = 1:length(mesh_tags)
        mesh_file  = [vox2cortex_meshroot,'\',sub_dir,mesh_tags{i_mesh},'.ply'];
        surf_ply = Surf_read_Data(mesh_file);
        surf_ply.Vertices = (surf_ply.Vertices./voxin_sz([2,1,3]) + 1/2).*mni_frac.shape + v_shift;
        surf_ply.vColor = [];
        save_file = [surf_dir,mri_tag,'.DL.',ply_names{i_mesh},'.ply'];
        Surf_save_File(surf_ply,save_file,'verbose',0)

        surf_plys{i_mesh} = surf_ply;
    end
    surf_plys(1) = [];

    T1_reg = MRI_read_Data([T1rb_dir,'\T1.reg.norm.nii.gz'],'verbose',0);
    LAB = mat_cut(3*T1_reg.data-1,0);

    [surf_plys_at,dWG,dGM,Gmean,d] = surfs_attach_grad(LAB,surf_plys,attach_opts);

    for i_mesh = 1:4
        save_file = [surf_dir,mri_tag,'.DL.',plyat_names{i_mesh},'.ply'];
        Surf_save_File(surf_plys_at{i_mesh},save_file,'verbose',0)
    end

    if if_show
        clf; lay_sub(2,2)
        
        next; slice_mri(dWG,plot_pos,'add_surf',surf_plys(1:2),'use_alpha',0); 
        view(2); legend off; axis off; title("d = "+Gmean{1}(1))
        next; slice_mri(dWG,plot_pos,'add_surf',surf_plys_at(1:2),'use_alpha',0);
        view(2); legend off; axis off; title("d = "+Gmean{1}(end))
        next; slice_mri(dGM,plot_pos,'add_surf',surf_plys(3:4),'use_alpha',0);
        view(2); legend off; axis off; title("d = "+Gmean{3}(1))
        next; slice_mri(dGM,plot_pos,'add_surf',surf_plys_at(3:4),'use_alpha',0);
        view(2); legend off; axis off; title("d = "+Gmean{3}(end))
        % lay_linkview;
        lay_fig([6,7]); sgtitle(sub_dir,'interpreter','none')

        exportgraphics(fig,[surf_dir,'/surfs_attach.png'])
    end

    t_use = toc;
    fprintf(round(t_use,1)+"s")
    fprintf('\n')
end
%% 
sk_vals = [0,4,5,6,7];

for i_sub = 1:length(i_list)
    t_tic = tic;
    % disp_wait(i_sub,length(i_list))
    sub_dir = sub_dirs(i_list(i_sub)).name;
    fprintf(i_sub+": sub"+i_list(i_sub)+" "+sub_dir+" => ")

    T1rb_dir = [FCD_rootdir,sub_dir,'\T1_rb\'];
    T1_reg = MRI_read_Data([T1rb_dir,'T1.reg.norm.nii.gz'],'verbose',0);
    T1_reg.data = mat_cut(3*T1_reg.data-1,0);
    T1_sk = T1_reg;

    surf_dir = [FCD_rootdir,sub_dir,'\surf\'];
    for i_mesh = 2:length(mesh_tags)
        surf_file = [surf_dir,mri_tag,'.DL.',plyat_names{i_mesh},'.ply'];
        surf_obj = Surf_read_Data(surf_file);

        skPial = Surf_voxelise(surf_obj,'mode',0,'ref',T1_reg);
        T1_sk.data(skPial.data == 1) = sk_vals(i_mesh);
    end
    MRI_save_File(T1_sk,[surf_dir,'T1_reg.pial_label.nii.gz']);
    t_use = toc(t_tic);
    fprintf(round(t_use,1)+"s")
    fprintf('\n')
end

function [surf_plys_at,dWG,dGM,Gmean,d] = surfs_attach_grad(LAB,surf_plys,opts_raw)
surf_plys_at = surf_plys;

% inner surface
WG = LAB;
WG(WG<1) = 1;
dWG = imgradient3(WG);

opts = [{'attach_grad','map',dWG},opts_raw];
[surf_plys_at{1},Gmean{1},d{1}] = Surf_func_ROI(surf_plys{1},opts{:});
[surf_plys_at{2},Gmean{2},d{2}] = Surf_func_ROI(surf_plys{2},opts{:});

% out pial
GM = LAB;
GM(GM>1) = 1;
dGM = imgradient3(GM);

opts = [{'attach_grad','map',dGM},opts_raw];
[surf_plys_at{3},Gmean{3},d{3}] = Surf_func_ROI(surf_plys{3},opts{:});
[surf_plys_at{4},Gmean{4},d{4}] = Surf_func_ROI(surf_plys{4},opts{:});
end