function [varargout] = plot_surf(Surface,varargin)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
options_default = structure('colors','auto','alphas',1,'names',[],'cmap_fix',[],'label_start',1,...
    'AmbientStrength',0.6,'DiffuseStrength',0.5,'EdgeColor',[],'show_vn',[],'show_vl',[],'show_zl','','Parent',[],'show_legend',1,'add_light',1,'add_btn',1,'v_func',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout');

is_cell = iscell(Surface);
if ~is_cell
    Surface = {Surface};
end

ax = gca;
fig = gcf;

first_plot =  strcmp(ax.NextPlot,'replace'); % || isempty(ax.Children);
if isempty(Parent)
    if first_plot
        %     tile_layout = tiledlayout(1,2,'TileSpacing','none','Padding','tight');
        %     nexttile
        ax.Position([1,3]) = [0.1,0.7];
    end
else
    [new_figure] = setParent_forPlot(Parent);
end

surfs_h = [];
auto_EdgeColor = isempty(EdgeColor);
for i = 1:length(Surface)
    if iscell(Surface{i})
        surf_h = plot_surf(Surface{i},varargin{:});
    else
        try
            V = Surface{i}.Vertices;
            F = Surface{i}.Faces;
        catch
            continue
        end
        
        % 坐标变换
        if ~isempty(v_func)
            V = v_func(Surface{i},V);
            Surface{i}.Vertices = V;
        end
        
        if ~ischar(alphas)
            if length(alphas)>1
                if length(alphas) == length(Surface)
                    alpha = alphas(i);
                else
                    alpha = @(varargin) alphas;
                end
            else
                alpha = alphas;
            end
        else
            alpha = alphas;
        end
        if isnumeric(alpha)
            alpha(alpha<0) = 0; alpha(alpha>1) = 1;
        end
        
        if isempty(names)
            if ~isempty(inputname(1)) && ~is_cell
                name = [inputname(1),'.~'];
            else
                if isfield(Surface{i},'Comment')
                    name = [Surface{i}.Comment,'.~'];
                else
                    name = '';
                end
            end
        else
            if ~iscell(names)
                names = {names};
            end
            name = names{i};
        end
        
        if_blend = isa(alpha,'function_handle');
        if if_blend
            alpha_func = alpha;
            alpha = 'flat';
        end
        
        % 自动判断显示Edge：如果faces比vertices还少，说明显示的为一条线
        if auto_EdgeColor
            if size(F,1)<size(V,1) && size(F,2)==3
                EdgeColor = 'flat';
            else
                EdgeColor = 'none';
            end
        end
        
        hold on
        surf_h=trisurf(F,V(:,1),V(:,2),V(:,3),...
            'EdgeColor',EdgeColor,'AmbientStrength',AmbientStrength,'DiffuseStrength',DiffuseStrength,...
            'FaceAlpha',alpha,'FaceLighting','gouraud',... % 'BackFaceLighting','unlit'
            'Tag',name,'DisplayName',name);
        
        surf_h.UserData.SurfAlpha = alpha;
        
        if if_blend
            if isfield(Surface{i},'Labels')
                x = Surface{i}.Labels;
            elseif isfield(Surface{i},'faceLabels')
                x = Surface{i}.faceLabels;
            else
                x = 1;
            end
            
            try
                x = alpha_func(Surface{i},x);
                x = double(x);
            catch
                x = amp(x);
            end
           
            x(x<0) = 0; x(x>1) = 1;
            surf_h.FaceVertexAlphaData = x;
        end
        
        % --------------------- color ------------------%
        if iscell(colors)
            color = colors{i};
        else
            color = colors;
        end

        if isempty(color)
            color = 'auto';
        end
        
        if isnumeric(color)
            if size(color,1) < min([size(V,1),size(F,1)]) % 单色
                if size(color,1) == length(Surface)
                    color = colors(i,:);
                end
                if length(color)==1 && color<1
                    color = color*[1,1,1];
                end
                if length(color)==1 && color>=1
                    surf_h.FaceVertexCData = color*ones(size(V,1),1);
                else
                    surf_h.FaceColor = color;
                end
            else
                % 直接染色每个V或F
                surf_temp = Surface{i};
                if size(color,1) == size(V,1)
                    surf_temp.Labels = color;
                    color_temp = 'V';
                end
                if size(color,1) == size(F,1)
                    surf_temp.faceLabels = color;
                    color_temp = 'F';
                end
                try
                    surf_h = coloring_surf(surf_temp,surf_h,color_temp,cmap_fix,label_start);
                end
            end
        end
        
        if ischar(color) % L,F或者colormap映射
            surf_h = coloring_surf(Surface{i},surf_h,color,cmap_fix,label_start);
        end
        
        % --------------------- additional ------------------%
        %%
        if isequal(EdgeColor,'flat') && isfield(Surface{i},'Labels')
            surf_h.CData = Surface{i}.Labels;
        end
        
        if isequal(surf_h.EdgeColor,'flat') && size(surf_h.FaceVertexCData,1) ~= size(V,1)
            % surf_h.EdgeColor = [1,1,1]*0.5; % 不然会显示错误
        end
        
        
        if ~isempty(show_vn)
            % plot the normals
            if isfield(Surface{i},'vNormals') && ~isempty(Surface{i}.vNormals)
                vNormals = Surface{i}.vNormals;
            else
                vNormals = -patchnormals(Surf_lower(Surface{i}));
            end
            
            if ~iscell(show_vn)
                show_vn = {show_vn,':'}; % if_show_vn, show_vn_roi
            end
            
            vn_l = show_vn{1}; vn_roi = show_vn{2};
            if all(size(vn_l) == size(V))
                given_vn = find(vecnorm(vn_l,2,2) > 0);
                vNormals(given_vn,:) = vn_l(given_vn,:);
                vn_bin = ':';
            else
                idx = 1:size(V,1);
                idx = idx(vn_roi);
                vn_bin = idx(1:vn_l:end);
            end
            
            hold on;
            h = quiver3(V(vn_bin,1),V(vn_bin,2),V(vn_bin,3),vNormals(vn_bin,1),vNormals(vn_bin,2),vNormals(vn_bin,3),0.8);
            h.DisplayName = [name(1:end-2),'_vn'];
        end
        
        if ~isempty(show_vl)
            hold on;
            vs = 3;
            if size(show_vl,1) == size(V,1)
                vl = show_vl(:,1);
                if size(show_vl,2) == 2
                    vs = show_vl(:,2);
                end
            else
                if isfield(Surface{i},'Labels')
                    vl = Surface{i}.Labels;
                else
                    vl = show_vl(1);
                end
                if isnumeric(vl) && size(vl,1) == 1
                    vl = double(vl)*ones(size(V,1),1);
                end
                if isnumeric(vl) && size(show_vl,2) == 2
                    vs = vl*show_vl(2);
                end
            end
            vs(vs == 0) = NaN;
            h = plot_scatter(V,'EdgeColor','none',vs,vl);
            h.DisplayName = [name(1:end-2),'_vl'];
        end
        
        if ~isempty(show_zl)
            if ~iscell(show_zl)
                show_zl = {show_zl};
            end
            
            zl_target = show_zl{1};
            zl_color = 'auto';
            if length(show_zl)>1
                zl_color = show_zl{2};
            end
            
            if strcmp(show_zl,'edge')
                zl_idx = 'edge';
                zl_target = 'F';
            else
                zl_idx = 'all';
            end
            zls = Surf_func_bw(Surface{i},'zone_line','idx',zl_idx,'target',zl_target);
            hold on; plot_line(zls,'color',zl_color,show_zl{3:end})
        end
    end
    
    surfs_h = [surfs_h,surf_h];
end
axis equal
axis tight
xlabel('X');ylabel('Y');zlabel('Z')

if add_light
    h_light=findobj(gca,'Type','light');
    if isempty(h_light)
        CA=[125,10];
        LightPos=[cosd(90+CA(1)),sind(90+CA(1)),-tand(CA(2))];
        Light1=light('Parent',gca,'Position',LightPos,'Style','infinite','Tag','Light1');
        Light2=light('Parent',gca,'Position',-LightPos,'Style','infinite','Tag','Light2');
        %     light('Position',[-1000 0 0],'Style','local')
        %     light('Position',[1000 0 0],'Style','local')
        %     camlight left
        %     camlight right
        % camlight('headlight')
    end
end

if abs(mean(fig.Color)-0.94)<0.01
    fig.Color = [1,1,1];
end

% shading interp
% lighting gouraud
% material dull
ax = gca;
ax.Clipping = 'off';

if fig.Position(3)/fig.Position(4)<1.5 && first_plot
    fig.Position(3) = fig.Position(4)*1.5;
end

if show_legend
    lgd = plot_legend(first_plot);
end
% if first_plot
%     %     lgd.Position(1) = 0.8*fig.Position(3);
% end

if add_btn
    % IconFile = fullfile(matlabroot,'/toolbox/matlab/icons/help_ex.png');
    % try
    %     tb = ax.Toolbar;
    %     btn_glow = axtoolbarbtn(tb,'push','Icon',IconFile,'Tooltip','update glow');
    % catch
    %     tb = axtoolbar(ax,{'rotate','datacursor','pan','zoomin','zoomout','restoreview','export'});
    %     btn_glow = axtoolbarbtn(tb,'push','Icon',IconFile,'Tooltip','update glow');
    % end
    % btn_glow.ButtonPushedFcn = @update_glow;

    % dcm = datacursormode; dcm.UpdateFcn = @plot_surf_tips;
end

rotate3d;

if nargout > 0
    varargout = {surfs_h};
end
end

function surf_h = coloring_surf(Surface,surf_h,colors,cmap_fix,label_start)
% FaceVertexCData
h_coloring = Surf_func_ROI(Surface,'coloring','colors',colors,'cmap_fix',cmap_fix,'label_start',label_start,'mode','plot');
h_fields = fieldnames(h_coloring);

if isempty(h_fields)
    surf_h.FaceColor = [1,1,1]*0.7;
else
    for i = 1:length(h_fields)
        surf_h.(h_fields{i}) = h_coloring.(h_fields{i});
    end
end
end

function cdata = get_cdata(Surface,label_field,cmap_fix)
if isfield(Surface,'colormap')
    cmap = Surface.colormap;
else
    cmap = colormap;
end
labels = Surface.(label_field);
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
            labels = amp(labels,0)*length(cmap);
        case 'amp'
            cmap_fix = cmap;
            labels = amp(labels)*length(cmap);
    end
end

cdata = squeeze(ind2rgb(round(labels),cmap_fix));
cdata(isnan(labels),:) = NaN;
end

function update_glow(src,evt)
ax_c = src.Parent.Parent;
axes(ax_c)
plot_glow('ax',ax_c)

CVM_rotate3d
warning('off','MATLAB:Axes:UpVector')
end