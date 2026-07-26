function [varargout] = Surf_func_ROI(Surface,func,varargin)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
func_short = '';
switch func
    case 'single'
        output = surf_numType(Surface,'single',varargin{:});
        func_short = '';
    case 'double'
        output = surf_numType(Surface,'double',varargin{:});
        func_short = '';
    case 'coloring'
        output = surf_coloring(Surface,varargin{:});
        func_short = '';
    case 'combine'
        output = surf_combine(Surface,varargin{:});
        func_short = 'comb';
    case 'join'
        output = surf_join(Surface,varargin{:});
        func_short = 'join';
    case 'divide'
        output = surf_divide(Surface,varargin{:});
        func_short = 'div';
    case 'split'
        output = surf_split(Surface,varargin{:});
        func_short = 'split';
    case 'intersect'
        surf_ind = varargin{1};
        [Ind] = SurfacesIntersect(Surface.Vertices,Surface.Faces,surf_ind.Vertices,surf_ind.Faces);
        output = Ind(:,[3:5]);
    case 'remove_isolated'
        output = surf_rm_isolated(Surface,varargin{:});
        func_short = 'rmiso';
    case 'remove_nan'
        output = surf_rm_nan(Surface,varargin{:});
        func_short = 'rmnan';
    case 'remove_unused'
        output = surf_rm_unuse(Surface,varargin{:}); % 去除F中没有的V
        func_short = 'rmuu';
    case 'reduce'
        output = surf_reduce(Surface,varargin{:});
        func_short = 'rdu';
    case 'resize8c'
        output = surf_resize8c(Surface,varargin{:});
        func_short = 'rdu8';
    case 'close8c'
        output = surf_close8c(Surface,varargin{:});
        func_short = 'clo8';
    case 'rotate'
        output = surf_rotate(Surface,varargin{:});
        func_short = 'rot';
    case 'refine'
        output = surf_refine(Surface,varargin{:});
        func_short = 'fine';
    case 'remesh'
         output = Surf_remesh(Surface,varargin{:});
         func_short = 're';
    case 'unify_normals'
        % 大的mesh极慢，可将mesh的F分为多块处理
        output = surf_unify_normals(Surface,varargin{:});
        func_short = 'uni';
    case 'smooth'
        output = Surf_smooth(Surface,varargin{:});
        func_short = 'sm';
    case 'boundary_smooth'
        output = surf_boundary_smooth(Surface,varargin{:});
        func_short = 'bsm';
    case 'peel_normal'
        output = surf_peel_normal(Surface,varargin{:});
        func_short = 'peel';
    case 'peel_lap'
        output = surf_peel_lap(Surface,varargin{:});
        func_short = 'peel';
    case 'peel_sphere'
        output = surf_peel_sphere(Surface,varargin{:});
    case 'peel_distgrad'
        output = surf_peel_distgrad(Surface,varargin{:});
        func_short = 'peel';
    case 'mapping_sphere'
        output = surf_mapping_sphere(Surface,varargin{:});
        func_short = 'sph';
    case 'mapping_rect'
        output = surf_mapping_rect(Surface,varargin{:});
        func_short = 'rect';
    case 'inBox'
        output = surf_inBox(Surface,varargin{:});
        func_short = 'box';
    case 'fill_hole'
        output = surf_fill_hole(Surface,varargin{:});
        func_short = 'fill';
    case 'fill_hole8c'
        output = surf_fill_hole8c(Surface,varargin{:});
        func_short = 'fill8c';
    case 'attach'
        output = surf_attach(Surface,varargin{:});
        func_short = 'at';
    case 'attach_edge'
        output = surf_attach_edge(Surface,varargin{:});
        func_short = 'ate';
    case 'attach_grad'
        output = surf_attach_grad(Surface,varargin{:});
    case 'uniform_sampling'
        output = surf_uniform_sampling(Surface,varargin{:});
        func_short = 'us';
    case 'grid_sampling'
        [output] = surf_space_sampling(Surface,varargin{:});
        func_short = 'gs';
    case 'is_inside'
        output = surf_inside(Surface,varargin{:});
        func_short = 'in';
    case 'sym_icp' % use icp to find a symmetrical, non-slanting posture
        output = surf_sym(Surface,'mode','icp',varargin{:});
        func_short = 'symicp';
    case 'sym_cpd' % use cpd to become symmetrical
        output = surf_sym(Surface,'mode','cpd',varargin{:});
        func_short = 'symcpd';
    case 'sph_alignAxis' % 将球脑中轴面对齐到YZ平面
        output = surf_sph_alignAxis(Surface,varargin{:});
    case 'partition' % 分块染色
        output = surf_partition(Surface,varargin{:});
        func_short = 'part';
    case 'slice'
        output = surf_slice(Surface,varargin{:});
        func_short = 'sli';
    otherwise
        disp('Process type error, must be ["combine \ divide \ split \ reduce \ refine \ smooth \ boundary_smooth \ peel \ peal_sphere \ inBox \ fill_hole \ attach]')
        return
end

% add comments, record the process history
if iscell(Surface) % if is cell, depend on the 1th one
    Surface = Surface{1};
end
if isfield(Surface, 'Comment')
    Comment = Surface.Comment;
else
    Comment = inputname(1);
end
varargout = Surf_comment_add(output,Comment,func_short);
end

function Surface = surf_numType(Surface,type,varargin)
fields = fieldnames(Surface);

eval(['type = @',type,';']);

for i = 1:length(fields)
    if isnumeric(Surface.(fields{i}))
        Surface.(fields{i}) = type(Surface.(fields{i}));
    end
end
end

function output = surf_coloring(Surface,varargin)
mode = '';
options_default = structure('colors','auto','cmap_fix',[],'label_start',1,'mode','make');
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if ismember(colors,{'label','l','L','v','V'})
    colors = 'labels';
end
if ismember(colors,{'face','f','fL','F'})
    colors = 'faceLabels';
end

surf_h = struct();
for iter = 1:2
    switch colors
        case 'auto'
            if isfield(Surface,'vColor') && ~isempty(Surface.vColor)
                if size(Surface.vColor,1)>1
                    surf_h.FaceVertexCData = Surface.vColor;
                else
                    surf_h.FaceColor = Surface.vColor;
                end
                colors = 'pass';
                continue
            end
            if isfield(Surface,'faceLabels') && ~isempty(Surface.faceLabels)
                colors = 'faceLabels';
                continue
            end
            if isfield(Surface,'Labels') && ~isempty(Surface.Labels)
                colors = 'labels';
                continue
            end
        case 'labels'
            if isfield(Surface,'Labels')
                if isempty(cmap_fix)
                    surf_h.FaceVertexCData = Surface.Labels;
                else
                    Labels = uint32(double(Surface.Labels) - label_start + 1); %??? => 对于uint32的0/1，+double的1后会产生计算错误
                    surf_h.FaceVertexCData = get_cdata(Surface,Labels,cmap_fix);
                end
            end
            colors = 'pass';
        case 'faceLabels'
            % shading flat
            faceLabels = uint32(double(Surface.faceLabels) - label_start + 1);
            surf_h.FaceVertexCData = get_cdata(Surface,faceLabels,cmap_fix);
        case 'pass'
        otherwise
            surf_h.FaceColor = colors;
            colors = 'pass';
    end
end

switch mode
    case 'plot'
        output = surf_h;
    case 'make'
        if isfield(surf_h,'FaceVertexCData') && ~isempty(surf_h.FaceVertexCData)
            Surface.vColor = surf_h.FaceVertexCData;
        end
        if ~isempty(cmap_fix)
            Surface.colormap = cmap_fix;
        end
        output = Surface;
end
end

function cdata = get_cdata(Surface,label_field,cmap_fix)
if isfield(Surface,'colormap')
    cmap = Surface.colormap;
else
    ax = gca;
    cmap = ax.Colormap;
end

if ischar(label_field)
    labels = Surface.(label_field);
else
    labels = label_field;
end

if isempty(cmap_fix)
    if isfield(Surface,'colormap')
        cmap_fix = 'auto0';
    else
        cmap_fix = 'amp';
    end
end

if ischar(cmap_fix)
    switch cmap_fix
        case 'auto0'
            cmap_fix = cmap; %cmap 从0开始
        case 'auto1'
            cmap_fix = [0,0,0;cmap]; %cmap 从1开始，需要补0
        case 'abp'
            cmap_fix = cmap;
            labels = amp(double(labels),0)*length(cmap);
        case 'amp'
            cmap_fix = cmap;
            labels = amp(double(labels))*length(cmap);
        otherwise
            eval(['cmap_fix = ',cmap_fix,'(',num2str(length(cmap)),');']);
    end
end

% 如果将 X 指定为整数数据类型的数组，则值 0 对应于颜色图 map 中的第一种颜色。
% 对于包含 c 种颜色的颜色图，图像 X 的值会被裁剪到范围 [0, c-1] 内。
% 如果将 X 指定为 single 或 double 数据类型的数组，则值 1 对应于颜色图中的第一种颜色。
% 对于包含 c 种颜色的颜色图，图像 X 的值会被裁剪到范围 [1, c] 内。
cdata = squeeze(ind2rgb(round(labels),cmap_fix)); % must be double
cdata(isnan(labels),:) = NaN;
end

function Surface_all = surf_combine(Surface,varargin)
options_default = struct('cmap','turbo','update_label','');
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

surf_num = length(Surface);
if surf_num<2
    Surface_all = Surface;
    return
end

if ischar(cmap)
    eval(['cmap = ',cmap,'(',num2str(surf_num),');'])
end

update_L = 0;
update_fL = 0;
if ischar(update_label)
    switch update_label
        case 'L'
            update_L = 1;
        case 'F'
            update_fL = 1;
        case 'all'
            update_L = 1;
            update_fL = 1;
    end
end

Surface_all = struct('Vertices',[],'Faces',[],'vColor',[],'Labels',[],'faceLabels',[]);
for i = 1:surf_num
    Surf = Surface{i};
    Nv = size(Surf.Vertices,1);
    Nf = size(Surf.Faces,1);
    
    if isempty(Surf)
        continue
    end
    
    Surface_all.Faces = [Surface_all.Faces;Surf.Faces+size(Surface_all.Vertices,1)];
    Surface_all.Vertices = [Surface_all.Vertices;Surf.Vertices];
    
    if ~isfield(Surf,'vColor')
        Surf.vColor = repmat(cmap(i,:),[Nv,1]);
    else
        if size(Surf.vColor,1) < Nv
            Surf.vColor = repmat(Surf.vColor(1,:),[Nv,1]);
        end
    end
    Surface_all.vColor = [Surface_all.vColor;Surf.vColor];
    
    if ~isfield(Surf,'Labels')
        if i == 1
            L = 1;
        else
            L = 1+max(Surface_all.Labels);
        end
        Surf.Labels = L;
    end

    if length(Surf.Labels) < Nv
        Surf.Labels = repmat(Surf.Labels(1),[Nv,1]);
    end

    if update_L
        Surf.Labels = ones(Nv,1)*i;
    end
    Surface_all.Labels = [Surface_all.Labels;Surf.Labels];
    
    if ~isfield(Surf,'faceLabels')
        if i == 1
            fL = 1;
        else
            fL = 1+max(Surface_all.faceLabels);
        end
        Surf.faceLabels = repmat(fL,[Nf,1]);
    end
    if update_fL
        Surf.faceLabels = ones(Nf,1)*i;
    end
    Surface_all.faceLabels = [Surface_all.faceLabels;Surf.faceLabels];
end
end

function output = surf_join(surfs_join,varargin)
options_default = structure('if_plot',1,'cmap','lines','alphas',0.5,'edge_smooth',10,'join_idx',{},'if_rematch',[],'if_recolor',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if length(surfs_join) == 1
    output = {surfs_join{1},[]};
    return
end

surf_join = Surf_func_ROI(surfs_join,'combine','cmap',cmap);

if_compute = isempty(join_idx);
if isempty(if_rematch)
    if_rematch = if_compute;
end
if if_compute
    join_idx = cell(2,1);
end

join_pos = cell(2,1);
for i = 1:2
    if if_compute
        surf_edge = Surf_func_bw(surfs_join{i},'zone_line','idx','edge','target','V');
        znlines{i} = surf_edge{1};
        id_valid = znlines{i}.idx_v>0;
        vIdx = znlines{i}.idx_v(id_valid);
        
        join_idx{i} = vIdx;
    end
    join_pos{i} = surfs_join{i}.Vertices(join_idx{i},:);
end

if if_rematch
    [join_Dist] = Net_nodes_distMap(join_pos{1},join_pos{2});
    [~,nearest_idx] = min(join_Dist,[],1);
    join_idx{2} = join_idx{2}(nearest_idx);
end

axes_func = @(x) [];

switch class(if_plot)
    case 'char'
        axes_func = if_plot;
        if_plot = 1;
    case 'function_handle'
        axes_func = if_plot;
        if_plot = 1;
end

if if_plot
    if ischar(cmap)
        eval(['cmap = ',cmap,'(6);'])
    end
end
for i = 1:2
    if if_plot
        next; plot_surf(surfs_join{i},'colors','L','alphas',alphas)
        
        zn_color = turbo(length(join_idx{i}));
        hold on; plot_scatter(join_pos{i},10,zn_color,'filled');
        axes_func(); caxis([0 size(cmap,1)]);
    end
end

reloc_idx = join_idx;
reloc_idx{2} = join_idx{2} + size(surfs_join{1}.Vertices,1);
remap_idx = 1:size(surf_join.Vertices,1);
temp = remap_idx; temp(reloc_idx{2}) = [];
remap_idx(temp) = 1:length(temp);
remap_idx(reloc_idx{2}) = reloc_idx{1};

surf_join = Surf_func_idx(surf_join,'remap_idx',{remap_idx});
surf_join = Surf_func_bw(surf_join,'zone_smooth','idx','all','iter',edge_smooth,'max_trans',2,'reg_bin',0);

if isempty(if_recolor)
    if_recolor = ~isfield(surf_join,'vColor');
end
if if_recolor
    surf_join = Surf_func_ROI(surf_join,'coloring','F',cmap);
end

if if_plot
    next; plot_surf(Surf_smooth(surf_join,'iter',10,'fix_id',reloc_idx{1}),'colors','F','alphas',alphas,'cmap_fix',cmap)
    axes_func(); caxis([0 size(cmap,1)]);
    
    colormap(cmap); lay_fig([12,4])
end

output = {surf_join,join_idx};
end

% 从surf减去子集
function [output] = surf_divide(Surface,varargin)
options_default = struct('idx',[],'target','','Labels','raw','faceLabels','raw','overlap',1,'remap',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

V = Surface.Vertices;
F = Surface.Faces;
Nv = size(V,1);
Nf = size(Surface.Faces,1);

switch target
    case 'V'
        if length(idx) < Nv
            temp = zeros(Nv,1);
            temp(idx) = 1;
            idx = temp>0;
        end
    case 'F'
        if length(idx) < Nf
            temp = zeros(Nf,1);
            temp(idx) = 1;
            idx = temp>0;
        end
end

unq_idx = unique(idx);
numSets = length(unq_idx);

if remap == 1 && numSets ~= 2
    error('numSets must be 2 when remap')
end

% divide_index最好为face的索引，如果不是(顶点的索引)，则需要推理
if length(idx) == Nv
    [idx_f] = Surf_func_idx(Surface,'V2F_idx',idx,'interp_type',1);
elseif length(idx) == Nf
    idx_f = idx;
else
    error('divide index error, it must be loggic index')
end

%% Create separate faces/vertices structures for each fSet
fvOut = repmat({struct('Faces',[],'Vertices',[])},numSets,1);

if ~isfield(Surface,'Labels')
    Surface.Labels = zeros(Nv,1);
end
if ~isfield(Surface,'vColor')
    Surface.vColor = [1,1,1]*0.7;
end

idx_v = cell(numSets,1);
for i = 1:numSets
    fvOut{i} = Surface;
    setF = F(idx_f == unq_idx(i),:);
    [unq_idx_v, ~, new_idx_v] = unique(setF);
    
    % 通常属性
    fields = fieldnames(Surface);
    for j = 1:length(fields)
        field = fields{j};
        if size(Surface.(field),1)>1
            if size(Surface.(field),1) == Nv
                fvOut{i}.(field)= Surface.(field)(unq_idx_v,:);
            end
            if size(Surface.(field),1) == Nf
                fvOut{i}.(field) = Surface.(field)(idx_f == unq_idx(i),:);
            end
        else
            fvOut{i}.(field) = Surface.(field);
        end
    end
    
    fvOut{i}.Faces = reshape(new_idx_v,size(setF));
    fvOut{i}.Vertices = V(unq_idx_v,:);
    idx_v{i} = unq_idx_v;
    
    if isfield(Surface,'faceLabels') && ~isempty(Surface.faceLabels)
        switch Labels
            case 'raw'
                fvOut{i}.faceLabels = Surface.faceLabels(idx_f == unq_idx(i),:);
            case 'id'
                fvOut{i}.faceLabels = i*ones(size(setF,1),1);
        end
    end
    
    switch Labels
        case 'raw'
            fvOut{i}.Labels = Surface.Labels(unq_idx_v,:);
        case 'id'
            fvOut{i}.Labels = i;
    end
    
    if remap% set the edge inex in front, so as to rejoin again
        znlines = Surf_func_bw(fvOut{i},'zone_line','idx','edge','target','F');
        zn_vIdx = nonzeros(znlines{1}.idx_v);
        
        remap_idx = 1:size(fvOut{i}.Vertices,1);
        remap_idx(zn_vIdx) = [];
        remap_idx= [zn_vIdx',remap_idx];
        fvOut{i} = Surf_func_idx(fvOut{i},'remap_idx',remap_idx);
        idx_v{i} = idx_v{i}(remap_idx);
        
        znlines = Surf_func_bw(fvOut{i},'zone_line','idx','edge','target','F');
        zn_vIdx = nonzeros(znlines{1}.idx_v);
    end
end

output = {fvOut,idx_f,idx_v};
end

function output = surf_split(Surface,varargin)
%% Organise faces into connected fSets that share nodes
V = Surface.Vertices;
F = Surface.Faces;

fSets = zeros(size(F,1),1,'uint32');
currentSet = 0;
while any(fSets==0)
    currentSet = currentSet + 1;
    %     fprintf('Connecting set #%d vertices...',currentSet);
    nextAvailFace = find(fSets==0,1,'first');
    openVertices = F(nextAvailFace,:);
    while ~isempty(openVertices)
        availFaceInds = find(fSets==0);
        [availFaceSub, ~] = find(ismember(F(availFaceInds,:), openVertices));
        fSets(availFaceInds(availFaceSub)) = currentSet;
        openVertices = F(availFaceInds(availFaceSub),:);
    end
    %     fprintf(' done! Set #%d has %d faces.\n',currentSet,nnz(fSets==currentSet));
end
numSets = currentSet;
%% Create separate faces/vertices structures for each fSet
% fvOut = repmat(struct('faces',[],'vertices',[]),numSets,1);
fvOut = cell(numSets,1);
BW = zeros(size(V,1),1);
for i = 1:numSets
    setF = F(fSets == i,:);
    [unqVertIds, ~, newVertIndices] = unique(setF);
    %     fvOut(currentSet).faces = reshape(newVertIndices,size(setF));
    %     fvOut(currentSet).vertices = V(unqVertIds,:);
    BW(unqVertIds) = i;
    
    Surface_div = Surf_func_ROI(Surface,'divide','idx',BW == i);
    fvOut{i} = Surface_div{end};
end

output = {fvOut,BW};
end

function output = surf_rm_isolated(Surface,varargin)
options_default = struct('area_min',50);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

[~,Lmap] = Surf_func_ROI(Surface,'split');
L = unique(Lmap);
L_area = hist(Lmap,L);
L_ismin = L(L_area<area_min);
is_small = ismember(Lmap,L_ismin);

[output] = surf_divide(Surface,is_small);
end

function output = surf_rm_nan(Surface,varargin)
options_default = struct('rm_v',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

V = Surface.Vertices;
idx_v = isnan(V(:,1));

if ~rm_v
    idx_f = Surf_func_idx(Surface,'V2F_idx',find(idx_v));
    Surface.Faces(find(idx_f),:) = [];
    %V(isnan(V)) = 0;
    Surface.Vertices = V;
else
    Surface = Surf_func_ROI(Surface,'divide',idx_v);
    Surface = Surface{1};
end

output = Surface;
end

function output = surf_rm_unuse(Surface,varargin)
V = Surface.Vertices; F = Surface.Faces;
v_used = V(:,1)*0;
v_used(unique(F(:))) = 1;
V(~v_used,:) = NaN;
Surface.Vertices = V;
output = surf_rm_nan(Surface,1);
end

function output = surf_reduce(Surface,varargin)
options_default = struct('bin',4,'limit',0,'with_fig',0,'use_fast',[],'if_remesh',0,'get_idx',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if limit>0
    bin = size(Surface.Vertices,1)/limit;
end

if with_fig
    figure
    vertices = Surface.Vertices;
    faces = Surface.Faces;
    p=trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),'EdgeColor','none','visible','off');
    p.FaceVertexCData = (1:length(p.FaceVertexCData))';
    p.Visible = 0;
else
    p = Surf_lower(Surface);
end

%%
if isempty(use_fast)
    use_fast = size(Surface.Faces,1) > size(Surface.Vertices,1); % 面比顶点还少说明是线条或者很分散
end

reduce_params = {};
switch use_fast
    case 1
        reduce_params = {'fast'};
end
%%
idx = [];
if bin>1
    if with_fig
        reducepatch(p,1/bin,reduce_params{:})
        Surface_reduce.Vertices = p.Vertices;
        Surface_reduce.Faces = p.Faces;
        
        ReducedVerticesNum = size(p.Vertices,1);
        p.FaceVertexCData=p.FaceVertexCData(1:ReducedVerticesNum);
        idx = p.FaceVertexCData;
    else
        [newF,newV] = reducepatch(p,1/bin,reduce_params{:});
        Surface_reduce.Vertices = newV;
        Surface_reduce.Faces = newF;
        
        if get_idx
            [idx] = knnsearch(p.Vertices,newV);
        end
    end
else
    output = {Surface,':'};
    return
end

% Surface_reduce.Labels = p.FaceVertexCData;
if isfield(Surface,'Labels')
    Surface_reduce.Labels = Surface.Labels(idx);
end
if isfield(Surface,'vColor') && size(Surface.vColor,1)>1
    Surface_reduce.vColor = Surface.vColor(idx,:);
end
% ReducedSurface.vNormals=p.VertexNormals;

if with_fig
    delete(gcf)
end

if if_remesh>0
    Surface_reduce = Surf_remesh(Surface_reduce,'iter',if_remesh);
end

output = {Surface_reduce,idx};
end

function output = surf_resize8c(Surface,varargin)
options_default = struct('sz',0.5,'if_close',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

reshape_Vs = @(V,n) set_val(reshape(V(:,1),n,[]),reshape(V(:,2),n,[]),reshape(V(:,3),n,[]));

Vs = Surface.Vertices;
Nv = size(Vs,1);
divN = sqrt(Nv);
[Xs,Ys,Zs] = reshape_Vs(Vs,divN);

% resize_func = @imresize;
resize_func = @(x,sz) imresize_grid(x,sz.*size(x));
Xs = resize_func(Xs,sz); Ys = resize_func(Ys,sz); Zs = resize_func(Zs,sz); 
temp = Surf_read_Data(surf2patch(Xs,Ys,Zs,'triangles'));

surf_rsz = Surface;
surf_rsz.Vertices = temp.Vertices;
surf_rsz.Faces = temp.Faces;

fields = fieldnames(Surface);
fields = setdiff(fields,{'Vertices'});
for i = 1:length(fields)
    field = fields{i};
    val = Surface.(field);
    if size(val,1) == Nv
        L = reshape(val,divN,divN,size(val,2)); L = imresize(L,sz);
        surf_rsz.(field) = reshape(L,[],size(val,2));
    elseif size(val,1) == 1
        surf_rsz.(field) = val;
    end
end

if isempty(if_close)
    if_close = 2*size(Vs,1) - size(Surface.Faces,1) == 4;
end

if if_close
    surf_rsz = surf_close8c(surf_rsz); surf_rsz = surf_rsz{1};
end

output = {surf_rsz};
end

function output = surf_close8c(Surface,varargin)
N = sqrt(size(Surface.Vertices,1));
d = 2/N;
[Xo,Yo,Zo] = Sph_map_base('mode','8c','pos',N,'if_plot',0,'lim',1-d);
sph_8c = Surf_read_Data(surf2patch(Xo,Yo,Zo,'triangles'));
sph_zl = Surf_func_bw(sph_8c,'zone_line'); sph_zl = sph_zl{1};

make_F = @(A,B) [A(1:end-1)',A(2:end)',B(1:end-1)';B(1:end-1)',A(2:end)',B(2:end)']; % 都顺时针
n = N/2-1;
F = [make_F(1:n+1,2*n+2:-1:n+2); make_F(2*n+2:3*n+2,4*n+3:-1:3*n+3); make_F(4*n+3:5*n+3,6*n+4:-1:5*n+4); make_F(6*n+4:7*n+4,[1,8*n+4:-1:7*n+5])];
F = [F; 1 2*n+2 4*n+3; 4*n+3 6*n+4 1];

F(:,1) = sph_zl.idx_v(F(:,1)); F(:,2) = sph_zl.idx_v(F(:,2)); F(:,3) = sph_zl.idx_v(F(:,3));
F = F(:,[2,1,3]); % 扭转法线方向
surf_close = Surface; surf_close.Faces = [surf_close.Faces; F];

output = {surf_close};
end

function Surface_refine = surf_refine(Surface,varargin)
options_default = struct('refine_rate',1,'blocks',1,'interp_type',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if prod(blocks)>1 % 启用分块再细化的方案，大幅加速
    v = Surface.Vertices;
    v_range = [min(v,[],1);max(v,[],1)];
    [block_IDX] = Mat_make_blocksIdx(v_range,blocks);
    disp_wait(0); process_over = 0; ii = 0; Surface_refine = cell(prod(blocks),1);
    while process_over < 1
        ii = ii+1;
        [sBlock] = Mat_make_blocksIdx(block_IDX,blocks,ii);
        [idx_x,idx_y,idx_z] = def_vars(sBlock,'x y z');
        [surf_roi,idx_v] = Surf_func_ROI(Surface,'inBox',[idx_x',idx_y',idx_z']);
        
        if sum(idx_v)>0
            surf_roi_refine = Surf_func_ROI(surf_roi,'refine','refine_rate',3,'blocks',1);
            Surface_refine{ii} = surf_roi_refine;
        end
        
        process_over = (ii == sBlock.n);
        disp_wait(sBlock.ii,sBlock.n)
    end
    Surface_refine = Surf_func_ROI(Surface_refine,'combine');
else
    if refine_rate>4
        refine_rate = 4;
    end
    
    [FV] = Surf_lower(Surface);
    for i = 1:refine_rate
        [FV2] = refinepatch(FV); %存储规则：前N个V不变，细化的新增在后面
        
        has_faceLabels = isfield(FV,'faceLabels') && ~isempty(FV.faceLabels);
        if has_faceLabels
            FV2.faceLabels = imresize(FV.faceLabels,[size(FV2.faces,1),1],'nearest');
        end
        
        FV = FV2;
    end
    Surface_refine = Surf_higher(FV);
end

Surface_refine.Vertices = real(Surface_refine.Vertices);
if has_faceLabels
    Surface_refine = Surf_func_idx(Surface_refine,'F2V_label','interp_type',interp_type);
end
end

function output = surf_rotate(Surface,varargin)
options_default = struct('azel',[0,0],'ang',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

% find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = pi*azel(1)/180;
    phi = pi*azel(2)/180;
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

alph = ang*pi/180;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
Trot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

surf_rot = Surface;
surf_rot.Vertices = Surface.Vertices*Trot;

output = surf_rot;
end

function output = surf_unify_normals(Surface,varargin)
options_default = struct('alighTo',1,'parts',[]); % 'in','out',or num
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

surf_raw = Surf_lower(Surface);
surf_part = surf_raw; surf_new = Surface;
F = surf_raw.Faces;

if isempty(parts)
    parts = round(size(F,1)/200);
end

part_idx = round(linspace(1,size(F,1),parts+1));
for i = 1:parts
    idx = part_idx(i):part_idx(i+1);
    surf_part.faces = F(idx,:);
    surf_part = unifyMeshNormals(surf_part,'alignTo',alighTo);
    surf_new.Faces(idx,:) = surf_part.faces;
    %fprintf('*')
end

output = surf_new;
end

function Surface = surf_boundary_smooth(Surface,varargin)
options_default = struct('nb_iterations',2,'ngb_degre',6,'max_trans',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

V = Surface.Vertices;
F = Surface.Faces;

% show_holes_and_boundary(V,T);
% shading interp;
boundaries = detect_mesh_holes_and_boundary(F);
V_out = smooth_mesh_boundaries(V,boundaries,nb_iterations,ngb_degre);

V_error_idx = mean(abs(V_out - V),2)>max_trans;
V_out(V_error_idx,:) = V(V_error_idx,:);

Surface.Vertices = V_out;
end

function Surface = surf_peel_normal(Surface,varargin)
options_default = struct('d',0,'iter',1,'sm_free',0,'fix_idx',[],'sm_vn',0,'sm_unif',5,'sm_all',2,...
    'given_vn',[],'reg_bin',0,'reg_res',5);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

Surface = Surf_func_ROI(Surface,'smooth','iter',sm_all);
for i = 1:iter
    vNormals = real(patchnormals(Surf_lower(Surface)));
    
    if ~isempty(given_vn)
        if_given = find(vecnorm(given_vn,2,2) > 0);
        vNormals(if_given,:) = given_vn(if_given,:);
    end
    
    if sm_vn>0
        surf_filt = Surface;
        surf_filt.Labels = vNormals;
        for j = 1:sm_vn
            evalc("[surf_filt] = Surf_func_filt(surf_filt,'V',@mean);");
        end
        vNormals = surf_filt.Labels;
    end
    
    V = Surface.Vertices;
    Vout = V + vNormals*d/iter;
    
    if reg_bin>0
        Vm = V(1:reg_bin:end,:);
        Vf = Vout(1:reg_bin:end,:);
        
        [~,T] = pts_register(Vm,Vf,'mode','bspline','res',reg_res);
        Vout = pts_warp(V,T);
    end
    Surface.Vertices = Vout;
    
    if sm_free>0
        Surface = Surf_func_ROI(Surface,'smooth','iter',sm_free,'fix_idx',fix_idx);
        Surface = Surf_func_ROI(Surface,'smooth','iter',sm_unif,'type',2,'fix_idx',fix_idx); % mesh vertices uniform distribute
        Surface = Surf_func_ROI(Surface,'smooth','iter',sm_all);
    end
end
end

function output = surf_peel_lap(Surface,varargin)
options_default = structure('pts_bin',2,'iter',4,'d','0.3','r_base',3,'res',6,'lap_opts',[2,1],'if_show',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

surf_bin = Surf_func_ROI(Surface,'resize8c',1/pts_bin);
surf_skel = surf_bin;
surf_skel.Vertices = Surf_lapSkeleton(surf_bin, structure('iter',iter,'WL',lap_opts(1),'WC',lap_opts(2)));

surf_skel = Surf_func_ROI(surf_skel,'resize8c',pts_bin);

mov_vec = surf_skel.Vertices-Surface.Vertices;
mov_d = vecnorm(mov_vec,2,2);
if ischar(d)
    dk = str2num(d);
    d = dk*max(mov_d);
else
    dk = d/max(mov_d);
end

for i = 1:length(d)
    rk = d(i)./mov_d;
    surf_peel{i} = Surface; surf_peel{i}.Vertices = (1-rk).*Surface.Vertices + rk.*surf_skel.Vertices;
end

if if_show
    plot_surf(Surface,'colors',[0.8,0.8,0.8],'alphas',0.25); hold on; plot_scatter(surf_skel.Vertices,2);
    hold on; plot_surf(surf_peel,'colors',(1:length(d))','alphas',0.5)
    hold on; plot_quiver(Surface.Vertices(1:66:end,:),mov_vec(1:66:end,:),0)
    caxis([1,length(d)]);
end

V = Surface.Vertices; v_std = mean(std(V,0,1)); divN = sqrt(size(V,1));
        
[Xo,Yo,Zo] = Sph_map_base('mode','8c','pos',divN,'if_plot',0,'lim',1);
sph_8c = Surf_read_Data(surf2patch(Xo,Yo,Zo,'triangles'));

sph_r = r_base*v_std;
vo = sph_r*sph_8c.Vertices; Vo = vo;

for i = 1:length(d)
    Vo = [Vo; vo*(1-dk(i))]; V = [V; surf_peel{i}.Vertices];
end

T_mulsph.sph_r = sph_r;
[~,T_mulsph.o2s] = pts_register(Vo,V,'bin',4,'mode','bspline','res',res);
[Vreg,T_mulsph.s2o] = pts_register(V,Vo,'bin',4,'mode','bspline','res',res);

if if_show
    next; plot_scatter(Vreg,3)
end

output = {surf_peel,T_mulsph,sph_8c};
end

function Surface = surf_peel_sphere(Surface,varargin)
options_default = structure('d',0.1,'v_sph',[],'T_sph',[],'if_offset',0,'offset_lim',[1/2,2]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if isempty(v_sph)
    v_sph = Surface.Vertices;
else
    if isstruct(v_sph) && isfield(v_sph,'T')
        v_sph = pts_warp(Surface.Vertices,v_sph,'bspline');
    end
end

if length(if_offset) == size(v_sph,1)
    if_offset(if_offset<offset_lim(1)) = offset_lim(1);
    if_offset(if_offset>offset_lim(2)) = offset_lim(2);
    d = d.*if_offset;
end

[Vs] = pts_warp(v_sph.*(1-d),T_sph,'bspline');

if isfield(T_sph,'rot')
    Vs = T_sph.cen + (Vs-T_sph.cen)/T_sph.rot;
end

Surface.Vertices = Vs;
end

function output = surf_peel_distgrad(Surface,varargin)
options_default = structure('map',[],'d',5,'iter',2,'smooth_dist',1,'ref',[],'smooth_surf',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if ~isempty(map)
    if isstruct(map)
        distmap = map.data;
    else
        distmap = map;
    end
else
    skind_fill = Surf_voxelise(Surface,'mode',1,'ref',ref);
    if mean(d)>0
        distmap = MRI_func_ROI(skind_fill.data==0,'dist_map');
    else
        distmap = MRI_func_ROI(skind_fill.data==0,'dist_map','if_inv',1,'ext',-1);
    end
    distmap = distmap.data;
    d = abs(d);
end

for i = 1:smooth_dist
    distmap = imgaussfilt3(distmap,1);
end
% [Gmag, Gaz, Gelev] = imgradient3(distmap);
% % Gmag = amp(Gmag,1e-3);
% Gmag = Gmag>0;
% [G{1},G{2},G{3}] = sph2cart(Gaz*pi/180, Gelev*pi/180,Gmag);
% G{2} = -G{2};

[G{1},G{2},G{3}] = imgradientxyz(distmap);
Gmag = 1e-5+sqrt(G{1}.^2+G{2}.^2+G{3}.^2);
G{1} = G{1}./Gmag; G{2} = G{2}./Gmag; G{3} = G{3}./Gmag;

v_raw = Surface.Vertices;
v_new = v_raw; G_new = 0*v_raw;

v_iters = repmat(v_raw,1,1,iter);
surf_sm = Surface;
for ii = 1:iter
    for i = 1:3
        G_new(:,i) = MRI_func_ROI(G{i},'sample',v_new);
    end

    if length(d) == 1
        len = d/iter; % 固定的总长度
    else
        len = d(:,ii);
    end

    v_new = v_new + G_new.*len;

    if smooth_surf>0
        surf_sm.Vertices = v_new;
        surf_sm = Surf_smooth(surf_sm,'iter',smooth_surf,'type',3);
        v_new = surf_sm.Vertices;
    end

    v_iters(:,:,ii) = v_new;
end

surf_inn = Surface; surf_inn.Vertices = v_new;
output = {surf_inn,G,v_iters,distmap};
end

function output = surf_attach_grad(Surface,varargin)
options_default = structure('map',[],'d',2,'iter',12,'d_decay',1.25,'sm_iter',1);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

[Gmag, Gaz, Gelev] = imgradient3(map);
[G{1},G{2},G{3}] = sph2cart(Gaz*pi/180, Gelev*pi/180,Gmag>0);
G{2} = -G{2};
Gmag_fuse = amp((2+amp(Gmag,0)).*exp(-amp(map,0)),0);

v_raw = Surface.Vertices;
v_new = v_raw; G_new = 0*v_raw;

iter = 12; d = 2;
Gmean = zeros(1,iter);
for ii = 1:iter
    for i = 1:3
        G_new(:,i) = MRI_func_ROI(G{i}.*Gmag_fuse,'sample',v_new);
    end
    Gmag_new = vecnorm(G_new,2,2);
    Gmean(ii) = mean(Gmag_new(:));
    v_new = v_new + G_new*d;
    d = d/1.25;
end

surf_new = Surface; surf_new.Vertices = v_new;
if sm_iter>0
    surf_new = Surf_smooth(surf_new,'iter',sm_iter,'type',3);
end

output = {surf_new, Gmean, d};
end

function output = surf_mapping_sphere(Surface,varargin)
options_default = structure('pre_rot',[],'pre_mov',[],'if_rot',1,'bin','2000','get_field',0,'colors',[],'fine',0,'r_rate',1,'post_func',[],'keep_area',0,'lock_zmax',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

Vs = double(Surface.Vertices);
Nv = size(Vs,1);

if ~isempty(pre_rot)
    v_cen = pre_rot{1};
    rotVecs = pre_rot{2};
    Vs = v_cen + (Vs-v_cen)*rotVecs;
end
if ~isempty(pre_mov)
    Vs = Vs + pre_mov;
end
Surface.Vertices = Vs;

%%
if fine>0 % 面数过少(不够平滑)时，共形变换不准确，例如对称脑曲面变换后两侧不对称
    surf_fine = Surf_func_ROI(Surface,'refine',1);
    Vk = surf_fine.Vertices; Fk = surf_fine.Faces;
else
    Vk = Vs; Fk = Surface.Faces;
end

Vc = mean(Vk,1);
Vd = mean(std(Vk,[],1));
Vks = (Vk-Vc)/Vd;

Vo = spherical_conformal_map(double(Vks),double(Fk));
if keep_area
    Vo = mobius_area_correction_spherical(double(Vks),double(Fk),Vo);
end

Vo = Vo(1:size(Vs,1),:);
Vo = 2*Vo*Vd+Vc;
sph_rot = Surface; sph_rot.Vertices = Vo;

if ischar(bin)
    pts_bin = max([1,round(Nv/str2double(bin))]);
else
    pts_bin = bin;
end

if if_rot % 尽可能旋转正
     surf_ref = Surface;
    if length(if_rot) > 1
        surf_ref.Vertices = if_rot;
    end
    
    % affine配准 ==> 用椭球去拟合
    [sph_fit,T_rot] = Surf_warp(sph_rot,surf_ref, 'mode','pointset_reg_affine','pts_bin',pts_bin);
    
    % 恢复为球
    Vf = sph_fit.Vertices; Vcf = mean(Vf,1);
    Rf = vecnorm(Vf-Vcf,2,2); Vf = (Vf-Vcf).*mean(Rf)./Rf + Vcf;
    sph_rot.Vertices = Vf;
    % figure; plot_surf_pair({sph_rot,surf_ref})
end

if lock_zmax
    [~,id] = max(Vks(:,3));
    p0 = Vo(id(1),:); % the selected vertex in the spherical parameterization
    p1 = [0,0,1]; % the target position
    % calculate cross and dot products
    C = cross(p0, p1) ;
    D = dot(p0, p1) ;
    NP0 = norm(p0) ; % used for scaling
    if ~all(C==0) % check for colinearity
        Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
        R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
    else
        R = sign(D) * (norm(p1) / NP0); % orientation and scaling
    end
    Vo = (R*Vo')';
end

Vo = r_rate*sph_rot.Vertices;
sph_r = mean(range(Vo,1))/2;
Vo = Vo - (max(Vo,[],1)+min(Vo,[],1))/2;

sph_rot.Vertices = Vo;
sph_rot.radius = sph_r;

if ~isempty(post_func)
    if ~iscell(post_func)
        post_func = {post_func};
    end
    for i = 1:length(post_func)
        cur_func = post_func{i};
        try
            sph_rot = cur_func(sph_rot);
        catch
            sph_rot.Vertices = cur_func(sph_rot.Vertices);
        end
    end
end
%%
T = [];
if get_field>0
    [~,T.o2s,reg_mse] = pts_register(Vo,Vs,'bin',pts_bin,'mode','bspline','res',get_field);
    [~,T.s2o,reg_inv_mse] = pts_register(Vs,Vo,'bin',pts_bin,'mode','bspline','res',get_field);
    
    if ~isempty(pre_rot)
        T.s2o.cen = v_cen; T.s2o.rot = rotVecs;
        T.o2s.cen = v_cen; T.o2s.rot = 1\rotVecs;
    end
end

output = {sph_rot,T};
end

function output = surf_mapping_rect(Surface,varargin)
V = Surface.Vertices; F = Surface.Faces;
Vr = rectangular_conformal_map(double(V),double(F));
surf_rect = Surface; surf_rect.Vertices(:,[1,2]) = Vr; surf_rect.Vertices(:,3) = 0;
output = surf_rect;
end

function output = surf_sph_alignAxis(Surface,varargin)
options_default = struct('axis_idx',[],'if_plot',0,'colors',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

Vo = Surface.Vertices;
[TH,PHI] = cart2sph(Vo(:,2),Vo(:,3),Vo(:,1)); % 轴在YZ平面

flat_bnd = Surface; Vop = [TH,PHI,0*TH]; flat_bnd.Vertices = Vop;

mov_idx = axis_idx; fix_idx = find(abs(Vop(:,2)) > 0.5); fix_idx = fix_idx(1:10:end);
mov_pts = Vop([mov_idx; fix_idx],[1,2]); ref_pts = mov_pts; ref_pts(1:size(mov_idx,1),2) = 0;
reg_c = [pi,pi/2];
T = []; [T.T,T.spacing,reg_pts] = point_registration([2*pi,pi],mov_pts+reg_c,ref_pts+reg_c,struct('MaxRef',4));
Vop_warp = bspline_trans_points_double(T.T,T.spacing,Vop(:,[1,2])+reg_c)-reg_c;
flat_bnd_warp = flat_bnd; flat_bnd_warp.Vertices(:,[1,2]) = Vop_warp;

[Vos(:,2),Vos(:,3),Vos(:,1)] = sph2cart(Vop_warp(:,1),Vop_warp(:,2),1+0*TH);
sphos_bnd = Surface; sphos_bnd.Vertices = Vos;

if if_plot
    color_val = colors; cmap_sign = jetw(255,'half');
    next; plot_surf(flat_bnd,'colors',color_val,'alphas',0.8); hold on; plot_scatter(Vop(mov_idx,:),'EdgeColor','none',10);
    colormapa(cmap_sign); caxis(1*[-1,1]); lighting none
    % next; plot_scatter({mov_pts,ref_pts,reg_pts-reg_c},'EdgeColor','none',10); axis tight; plot_legend
    next; plot_surf(flat_bnd_warp,'colors',color_val,'alphas',0.8); hold on; plot_scatter(Vop_warp(mov_idx,:),'EdgeColor','none',10);
    colormapa(cmap_sign); caxis(1*[-1,1]); lighting none
    next; plot_surf(sphos_bnd,'colors',color_val,'alphas',0.8); colormapa(cmap_sign); caxis(1*[-1,1]);
end

output = sphos_bnd;
end

function output = surf_inBox(Surface,varargin)
options_default = struct('pos',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

v = Surface.Vertices;
idx_v = (v(:,1)>=pos(1,1)) & (v(:,1)<pos(2,1));
idx_v = idx_v & (v(:,2)>=pos(1,2)) & (v(:,2)<pos(2,2));
idx_v = idx_v & (v(:,3)>=pos(1,3)) & (v(:,3)<pos(2,3));

if sum(idx_v)>0
    [Surface,idx_f,idx_v] = Surf_func_ROI(Surface,'divide',idx_v);
    Surface = Surface{2};
    idx_v = idx_v{2};
else
    Surface = [];
    idx_f = [];
    idx_v = [];
end

output = {Surface,idx_f,idx_v};
end

function output = surf_inside(Surface,varargin)
options_default = struct('pos',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

V_roi = Surface.Vertices;
box = [min(V_roi,[],1); max(V_roi,[],1)];

idx_v = (pos(:,1)>=box(1,1)) & (pos(:,1)<box(2,1));
idx_v = idx_v & (pos(:,2)>=box(1,2)) & (pos(:,2)<box(2,2));
idx_v = idx_v & (pos(:,3)>=box(1,3)) & (pos(:,3)<box(2,3));

in_idx = in_polyhedron(Surf_lower(Surface),pos(idx_v,:));
is_inside = idx_v; is_inside(is_inside>0) = in_idx;
output = {is_inside};
end

function output = surf_fill_hole(Surface,varargin)
if iscell(Surface)
    for i = 1:length(Surface)
        Surface{i} = Surf_func_ROI(Surface{i},'fill_hole',varargin{:});
    end
    Surface = {Surface};
    return
end

V = Surface.Vertices;
F = Surface.Faces;
Nv = size(V,1); Nf = size(F,1);
boundaries = detect_mesh_holes_and_boundary(F);

V_fill = zeros(length(boundaries),3);
F_fill = [];
for i = 1:length(boundaries)
    v_idx = (boundaries{i})';
    V_hole = V(v_idx,:);
    
    %     DT = delaunayTriangulation(V_hole); % 不是三角剖分，而是四面体剖分..
    %     figure; tetramesh(DT,'FaceAlpha',0.3);
    
    V_c = mean(V_hole,1); % 边界每条边都与中心相连即可
    V_fill(i,:) = V_c; % 新增一个顶点
    F_hole = [v_idx,v_idx([2:end,1]),ones(length(v_idx),1)*(Nv+i)];
    F_fill = [F_fill;F_hole];
end

Nv_fill = size(V_fill,1); Nf_fill = size(F_fill,1);
Surface.Vertices = [Surface.Vertices; V_fill]; Surface.Faces = [Surface.Faces; F_fill];

fields = fieldnames(Surface);
fields(ismember(fields,{'Vertices','Faces'})) = [];
for i = 1:length(fields)
    val = Surface.(fields{i});
    if isnumeric(val)
        if size(val,1) == Nv
            val = [val; nan(Nv_fill,size(val,2))];
        end
        if size(val,1) == Nf
            val = [val; nan(Nf_fill,size(val,2))];
        end
    end
    Surface.(fields{i}) = val;
end

output = {Surface,V_fill,F_fill};
end

function Surface = surf_fill_hole8c(Surface,varargin)
Vs = Surface.Vertices;
divN = sqrt(size(Vs,1));
[~,~,~,~,temp] = Sph_map_base('mode','8c_close','pos',divN,'if_plot',0); seal_faces = temp.seal_faces;
Surface.Faces = [Surface.Faces; seal_faces];
end

function output = surf_attach(Surface,varargin)
options_default = struct('pos',[],'bin',1,'roi',[],'extend',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

v = Surface.Vertices;
if extend>0
    if ~isfield(Surface,'vNormals')
        Surface.vNormals = real(patchnormals(Surf_lower(Surface)));
    end
    v = v*(1+extend*Surface.vNormals);
end

if isempty(roi)
    roi = [min(pos,[],1); max(pos,[],1)];
end

idx_ref = 1:bin:size(v,1); vbin = v(idx_ref,:);
is_out = vbin(:,1)<roi(1,1) | vbin(:,2)<roi(1,2) | vbin(:,3)<roi(1,3);
is_out = is_out | vbin(:,1)>roi(2,1) | vbin(:,2)>roi(2,2) | vbin(:,3)>roi(2,3);
idx_ref(is_out) = []; vbin = v(idx_ref,:);

% projMid_dist = Net_nodes_distMap(pos,'ref',vbin);
% [~,idx_proj] = min(projMid_dist,[],2);
[idx_proj] = knnsearch(vbin,pos);

idx_proj = idx_ref(idx_proj); idx_proj = idx_proj(:);
v_proj = v(idx_proj,:);
output = {v_proj,idx_proj};
end

function output = surf_attach_edge(Surface,varargin)
ref = []; type = [];
options_default = structure('ref',[],'type','all');
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

zl = Surf_func_bw(Surface,'zone_line','idx','edge'); zl = zl{1};
idx = find(zl.idx_v>0); pos = zl.pos(idx,:);
disp_map = Net_nodes_distMap(pos,ref);

switch type
    case 'all' % 所有贴到最靠近ref的地方
        [~,idx_match] = min(disp_map,[],2);
        pos_new = ref(idx_match,:);
    case 'each' % 离ref最近的几个点贴上去
        [~,idx_match] = min(disp_map,[],1);
        pos_new = pos;
        pos_new(idx_match,:) = ref;
end

Surface.Vertices(zl.idx_v(idx),:) = pos_new;
zl_new = zl; zl_new.pos(idx,:) = pos_new;
output = {Surface,zl,zl_new};
end

function P = surf_uniform_sampling(Surface, N)
% Generates N uniformly sampled points on a triangulated surface based on
% the area of each triangle.
vertices = Surface.Vertices';
faces = Surface.Faces';

face_num=size(faces,2); % number of faces
% Calculation of area of faces for sampling
area = zeros(face_num, 1);
for i=1:face_num
    t=faces(:,i);
    tcorr=vertices(:,t); % coordinates of all three vertices of i-th face
    u=tcorr(:,2)-tcorr(:,1); % vector joining vertex 1 to vertex 2
    v=tcorr(:,3)-tcorr(:,1); % vector joining vertex 1 to vertex 3
    area(i)=0.5*norm(cross(u,v)); % area of i-th face
end
ca=cumsum(area); % Cumulative sum of area in an array
total_area=ca(end);
f=size(faces,2);
P = zeros(3,N);
for i=1:N
    % Choosing face based on area of faces
    r=rand;
    ra=r*total_area; % random number between 0 and total_area
    % choosing a face based on random number ra
    m=1;n=f;
    while(abs(n-m)>2)
        j=m+round((n-m)/2);
        if ca(j)>=ra; n=j;
        else; m=j;
        end
    end
    if ra<=ca(n-1)
        if ra<=ca(m); n=m;
        else; n=n-1;
        end
    end
    
    t=faces(:,n);
    % A,B,C are the coordinates of three vertices of n-th the face
    A=vertices(:,t(1));     B=vertices(:,t(2));     C=vertices(:,t(3));
    r1=rand;r2=rand;
    % Based on paper "Shape Distributions", Osada et al.
    P(:,i)=(1-sqrt(r1))*A+sqrt(r1)*(1-r2)*B+sqrt(r1)*r2*C;
end

P = P';
end

function [output] = surf_space_sampling(Surface,varargin)
options_default = structure('d',10,'box',[],'dist_lim',20);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

[X,Y,Z] = meshgrid(box(1,1):d:box(2,1),box(1,2):d:box(2,2),box(1,3):d:box(2,3));
pts = [X(:),Y(:),Z(:)];
pts_dist = Surf_func_Val(Surface,'pts_dist','v',pts,'mode','p2p','pts_bin',10);
rm_idx = pts_dist<dist_lim; rm_pts = pts(rm_idx,:); pts(fix_idx,:) = [];

output = {pts,rm_pts};
end

function output = surf_sym(Surface,varargin)
mode = '';
options_default = structure('mode','icp','flip_cen',[],'flip_rot',[],'get_err',0,'if_plot',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

use_icp = 0; use_cpd = 0;
switch mode
    case 'icp'
        use_icp = 1;
    case 'cpd'
        use_cpd = 1;
end

V = Surface.Vertices;
V_std = std(V,0,1);

if isempty(flip_cen)
    flip_cen = mean(V,1);
end
if isempty(flip_rot)
    flip_rot = diag([-1,1,1]);
else
    if ischar(flip_rot)
        switch flip_rot
            case 'x'
                flip_rot = diag([-1,1,1]);
            case 'y'
                flip_rot = diag([1,-1,1]);
            case 'z'
                flip_rot = diag([1,1,-1]);
        end
    end
end

surf_flip = Surface; surf_sym = Surface;
surf_flip.Vertices = flip_cen + (Surface.Vertices-flip_cen)*flip_rot;
reg_err = 0.8;

%% rigid transform using icp registration between raw and flip surface
if use_icp
    [surf_reg,rotVecs] = Surf_warp(Surface,surf_flip,'pointcloud_reg_icp','pc_para',{'Metric','pointToPlane'}); % must use pointToPlane !
    
    % notice: rotVecs = (rotVecs + diag([1,1,1]))/2; % error, not rotation
    rotVecs = sqrtm(rotVecs); % T1*T1 = T2;
    surf_sym.Vertices = flip_cen + (Surface.Vertices-flip_cen)*rotVecs;
    
    T.rot = rotVecs; T.cen = flip_cen;
end

%% nonrigid transform using cpd registration between raw and flip surface
if use_cpd
    % build-in default paras:
    % 'MaxIterations',20,'OutlierRatio',0.1,'InteractionSigma',2,'SmoothingWeight',3 => 7s
    [surf_reg] = Surf_warp(surf_sym,surf_flip,'pointcloud_reg_cpd','norm_init',0, ...
        'param',{'MaxIterations',100,'OutlierRatio',0.1,'InteractionSigma',1e-4,'SmoothingWeight',0.5}); % 1e-4 for fast, 4s
    surf_reg = Surf_smooth(surf_reg,'iter',5,'type',1,'attach_surf',[1,3],'ref_surf',surf_flip); % => 3s
    
    surf_sym.Vertices = (surf_sym.Vertices + surf_reg.Vertices)/2;
    movVecs = surf_sym.Vertices - V;
    
    T.mov = movVecs; T.cen = flip_cen;
end

if get_err
    reg_err = Surf_func_Val(surf_reg,'pts_dist','v',surf_flip.Vertices);
end

if if_plot
    next; plot_surf(Surface,'colors',0.8,'alphas',0.25);
    hold on; plot_surf(surf_flip,'colors',[1,0.8,0.8],'alphas',0.5);
    hold on; plot_surf(surf_reg,'colors',[0.8,1,0.8],'alphas',0.8);
    next; plot_surf(surf_sym,'colors',reg_err,'alphas',1);
    if get_err
        colormapa(0.9*jetw(64,'half')); caxis([-1,1]*mean(V_std)/20); shading interp
    end
end

output = {surf_sym,T,reg_err};
end

function output = surf_partition(Surface,varargin)
options_default = structure('num',100,'D_gamma',-1,'if_plot',0);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

[f_s,f_cen] = Surf_func_Val(Surface,'area');

[A] = compute_dual_graph(double(Surface.Faces),double(Surface.Vertices));
D = Net_nodes_distMap(f_cen);
D = amp(exp(D_gamma*D).*A,0);

num_big = ceil(sqrt(num));
[L_mat,~,~,L_tree] = Net_func_ROI(D,'cluster',[num,num_big],2); % 分两级

labels = L_mat(:,1);
C = full(sparse(1:size(labels, 1), labels, 1));
clu_vSum = sum(C,1);
clu_sSum = sum(C.*f_s,1);

cmap_gap = 0.7*jet(num_big^2) + 0.3*repmat(get_cmaps(1,num_big),[num_big,1]); % 粗级别平滑，细级别变化
if if_plot
    next; colormapa(cmap_gap); plot_surf(Surface,'colors',labels);
    next; hist(clu_vSum); title('v num'); axis square
    next; hist(clu_sSum); title('area'); axis square
end

info = structure('vSum',clu_vSum,'sSum',clu_sSum,'L_tree',L_tree,'cmap',cmap_gap);
output = {L_mat,info};
end

function output = surf_slice(Surface,varargin)
options_default = structure('dim','z','idx',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

V = Surface.Vertices;
switch dim
    case 'x'
        V = V(:,[2,3,1]);
    case 'y'
        V = V(:,[1,3,2]);
    case 'z'
end

evalc("stlwrite('CurrentRegion.stl',Surface.Faces,V);");
triangles = stlread_triangles('CurrentRegion.stl');

slices_pointset = stl_slicing(triangles,idx);

output = {slices_pointset};
end
