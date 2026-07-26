function [varargout] = plot_mri(mri,varargin)
%UNTITLED6 此处显示有关此函数的摘要
%   position: xyz
cmap = []; info = []; zoom = [];
options_default = structure('act','on','pos',[],'ctype','index','coordinate','ACS','voxsize',[],'parent',[],'cmap','gray','cAxis','cmap','cLim',[],...
    'interp_method','auto','if_detrend',0,'if_contour',0,'tips_info',[],'cross_width',1,'cross_gap','0.05','cross_color',[],'cycle_width',0,...
    'moving_func',[],'enhance',1,'info',[],'titles',[],'verbose',0,'zoom',1,'if_proj',[],'thresh',[]);
[options, eval_str] = resolve_input(options_default,varargin);
eval(eval_str);

if iscell(mri)
    mri = mri'; mri = mri(:);
    switch act
        case 'win'
            plot_mri_win(mri,info)
            return
        otherwise
            switch act
                case 'on'
                    act = 'off';
            end
            for i = 1:length(mri)
                if iscell(cmap)
                    cur_cmap = cmap{i};
                else
                    cur_cmap = cmap;
                end
                if ~isempty(mri{i})
                    plot_mri(mri{i},varargin{:},'act',act,'cmap',cur_cmap);
                end
                if iscell(titles)
                    if i == 1
                        titles = titles'; titles = titles(:);
                    end
                    title(titles{i})
                end
                if i<length(mri)
                    nexttile;
                    axis off
                end
            end
            return
    end
end

sMRI = mri;
voxsize = [1,1,1];
if isstruct(mri)
    mri = mri.data;
    if isempty(mri)
        return
    end
    if isempty(voxsize)
        try
            voxsize = sMRI.voxsize;
        end
    end
end
% mri(isnan(mri)) = 0;

if ndims(mri)==4 && size(mri,4)>size(mri,3)
    mri = permute(mri,[1,2,4,3]);
end

if ~isempty(thresh)
    mri = mri-thresh;
    mri(mri<0) = 0;
end

if strcmp(act,'off')
    if isempty(parent)
        ax = gca;
        if isempty(ax.Children)
            parent = gca;
        end
    end
end

if if_detrend && size(mri,4)>1 % used for visualize deformfield
    mri_raw = mri;
    mri(isnan(mri)) = 0;
    
    [ny,nx,nz,nc] = size(mri);
    [X,Y,Z] = meshgrid(1:nx,1:ny,1:nz);
    mask = mean(mri,4)>0;
    for i = 1:nc
        mri_val = mri(:,:,:,i);
        regress_val = [X(mask),Y(mask),Z(mask),1+X(mask)*0];
        regress_b = regress(mri_val(mask),regress_val);
        mri_val = mri_val - (regress_b(1)*X + regress_b(2)*Y + regress_b(3)*Z + regress_b(4));

        mri_val = (mri_val - min(mri_val(mask))).*mask;
        mri(:,:,:,i) = mri_val;
    end
    mri_max = max(mri,[],[1,2,3]);
    % temp = mri(round(0.3*end):round(0.7*end),round(0.3*end):round(0.7*end),round(0.3*end):round(0.7*end),:); mri_min = min(temp,[],[1,2,3]);
    temp = reshape(mri,[],size(mri,4));
    temp(mean(temp,2)==0,:) = [];
    temp = prctile(temp,10,1);
    mri_min = permute(temp,[1,4,3,2]);
    mri = (mri-mri_min)./(mri_max-mri_min)*255;
    mri(mri<0) = 0; mri(mri>255) = 255;
    mri = mri/255;
end

char_cmap =  ischar(cmap);
if char_cmap
    cmap = str2func(cmap);
    cmap = cmap();
    cmap(1,:) = [0,0,0];
end

if strcmp(interp_method,'auto')
    color_diff = mean2(abs(cmap(1:end-1,:) - cmap(2:end,:)));
    if color_diff<0.005
        interp_method = 'bilinear';
    else
        interp_method = 'nearest';
    end
end

plot_dim = '';
only_1dim = ismember(act,{'x','y','z'});
if only_1dim
    plot_dim = act;
    act = 'off';
end

switch act
    case 'win'
        plot_mri_win(mri,info)
    case 'on'
        if isempty(parent)
            fig = figure;
        else
            fig = parent;
        end

        mri(isnan(mri)) = 0;
        mri(isinf(mri)) = 0;
        os_h = orthosliceViewer(mri,'Parent',fig,'ScaleFactors',[1,1,1],'DisplayRangeInteraction','off');
        fig.UserData.orthoslice = os_h;
        fig.UserData.ortho = [];
        fig.UserData.ortho.vol = sMRI;
        fig.UserData.ortho.locs = [];
        fig.UserData.ortho.size = size(sMRI);
        fig.UserData.ortho.voxsize = voxsize;

        if isempty(parent)
            fig.Position(2) = 10;
        end

        scale_rate = mean(fig.Position(3:4))/600;
        fig.Position(3:4) = fig.Position(3:4)/scale_rate;
        panel_h = os_h.Parent;

        if strcmp(coordinate,'ACS')
            [hXY hYZ hXZ] = getAxesHandles(os_h);
            panel_h.SizeChangedFcn = '';

            hXY.Units = 'normalized'; hYZ.Units = 'normalized'; hXZ.Units = 'normalized';
            hXZ.YDir = 'normal'; hYZ.YDir = 'normal'; hXY.YDir = 'normal';
            hXZ.Tag = 'XZ'; hYZ.Tag = 'YZ'; hXY.Tag = 'XY';

            axes(hYZ);
            view([90 -90])

            % hXZ.DataAspectRatioMode = 'Auto'; hYZ.DataAspectRatioMode = 'Auto'; hXY.DataAspectRatioMode = 'Auto';
            py_sz = size(mri,[1:3]).*voxsize;
            % hXZ.DataAspectRatio = py_sz([3,2,1]);
            % hYZ.DataAspectRatio = py_sz([3,1,1]);
            % hXY.DataAspectRatio = py_sz([2,1,1]);

            hXZ.Position([3,4]) = hXZ.Position([3,4]);
            hYZ.Position([3,4]) = hYZ.Position([3,4]);
            hXY.Position([3,4]) = hXY.Position([3,4]);

            hXY.Position(2) = 0.05;
            hYZ.Position(2) = hXY.Position(2)+hXY.Position(4)+0.05;
            hXZ.Position(2) = hYZ.Position(2);
            hYZ.Position(4) = hXZ.Position(4);
            hYZ.Position(1) = hXZ.Position(1)+hXZ.Position(3)+0.05;

            fig.WindowScrollWheelFcn = @figScroll;
        end

        if ~isempty(if_proj)
            projZ_imgs = 0*mri;
            LH = os_h.DisplayRange;
            
            mri_max = max(mri(1:2:end));
            projImgs = plot_projImgs(mri,'if_plot',0,'dim','xyz','if_grad',0,'kb',if_proj,'LH',LH,'cmax',mri_max,'samp',300,'fuse',0.7);

            hXY.UserData.proj_imgs = projImgs.zImgs;
            hYZ.UserData.proj_imgs =  projImgs.xImgs;
            hXZ.UserData.proj_imgs =  projImgs.yImgs;
            moving_func = @moving_proj;

            img_interp = interp_method;
            hXY.Children(end).Interpolation = img_interp;
            hYZ.Children(end).Interpolation = img_interp;
            hXZ.Children(end).Interpolation = img_interp;

            fig.UserData.ortho.vol_max = gather(mri_max);

            clear('a_all')
            % mri = gather(mri);
        end

        if ~isempty(tips_info)
            os_h.DisplayRangeInteraction = 'off';

            region_id_h = uicontrol('Style','edit','String','id','Parent',panel_h,'Units','normalized','Tag','region_id');
            region_id_h.Position = [hYZ.Position(1),hXY.Position(2)+0.07,hYZ.Position(3)*0.35,0.05];

            region_tag_h = uicontrol('Style','edit','String','region','Parent',panel_h,'Units','normalized','Tag','region_tag');
            region_tag_h.Position = [hYZ.Position(1)+hYZ.Position(3)*0.4,hXY.Position(2)+0.07,hYZ.Position(3)*0.6,0.05];

            region_fullname_h = uicontrol('Style','edit','String','region fullname','Parent',panel_h,'Units','normalized','Tag','region_fullname');
            region_fullname_h.Position = [hYZ.Position(1),hXY.Position(2),hYZ.Position(3),0.05];

            if isempty(moving_func)
                moving_func = @moving_log;
            end
            addlistener(os_h,'CrosshairMoved',@(src,evt,moving_func) disp_info(src,evt,moving_func));
        end
        if ~isempty(moving_func)
            addlistener(os_h,'CrosshairMoving',moving_func);
            addlistener(os_h,'CrosshairMoved',moving_func);
            fig.UserData.moving_func = moving_func;
        end

        colormap(cmap);
        %         os_h.DisplayRange = [0,size(cMap,1)];

        m = uimenu(fig,'Text','&Plot_MRI');
        mitem = uimenu(m,'Text','Clear location log'); mitem.MenuSelectedFcn = @(varargin) clear_locLog(varargin{:});
        mitem = uimenu(m,'Text','Plot location log'); mitem.MenuSelectedFcn = @(varargin) plot_locLog(varargin{:});
    case 'off'
        mri_4d = mri;
        if ndims(mri) == 4
            chs_max = 3;
        else
            chs_max = 1;
        end

        imgsize = size(mri,1:3);
        if isempty(pos)
            pos = round([0.5,0.46,0.5].*imgsize);
            plot_type = 'slice';
        else
            if ischar(pos) && strcmp(pos,'max_proj')
                pos = [1,1,1];
                plot_type = 'max_projection';
            else
                pos = pos([2,1,3]);
                % relative pos
                idx = pos<1; pos(idx) = pos(idx).*imgsize(idx);
                plot_type = 'slice';
            end
        end

        pos = round(pos);

        img_mtg = [];
        for chs = 1:chs_max
            mri = mri_4d(:,:,:,chs);

            if strcmp(coordinate,'ACS')
                imgXZ = []; imgYZ = []; imgXY = [];
                switch plot_type
                    case 'slice'
                        imgXZ = rot90(squeeze(mri(pos(1),:,:)),1);
                        imgYZ = rot90(squeeze(mri(:,pos(2),:)),1);
                        imgXY = flipud(mri(:,:,pos(3)));
                    case 'max_projection'
                        [imgYZ,imgXZ,imgXY] = make_projection_map(mri);
                        imgYZ = rot90(imgYZ,2);
                        imgXZ = flipud(imgXZ);
                        imgXY = flipud(imgXY);
                end
            end

            if ~isempty(voxsize)
                array2str = @(x) strjoin(split(num2str(x)),', ');
                voxsize_norm = voxsize/min(voxsize);
                imgsize_re = round(imgsize.*voxsize_norm) ;

                if chs == 1 && any(voxsize ~= 1) && verbose
                    disp(['volume size = [',array2str(size(mri)), '], voxsize = [',array2str(voxsize),']',';  show size = [',array2str(imgsize_re),']', ', pixsize = ',array2str(min(voxsize))]);
                end

                imgXZ = imresize(imgXZ,imgsize_re([3,2]),interp_method);
                imgYZ = imresize(imgYZ,imgsize_re([3,1]),interp_method);
                imgXY = imresize(imgXY,imgsize_re([1,2]),interp_method);
                pos_re = pos.*voxsize_norm;
            else
                imgsize_re = imgsize;
            end

            img_cat = [imgXZ,imgYZ;imgXY, zeros(imgsize_re(1))];
            img_mtg = cat(3,img_mtg,img_cat);
        end


        if isnumeric(enhance)
            img_mtg = enhance*img_mtg;
        else
            img_mtg = enhance(img_mtg);
        end

        img_labels = flipud([imgXZ*0+1,imgYZ*0+2;imgXY*0+3,zeros(imgsize_re(1))+4]);
        %         if ischar(cMap)
        %             eval(['cMap = ',cMap,';']);
        %             cMap(1,:) = [0,0,0];
        %         end

        if isempty(parent)
            ax = gca;
            need_new_figure = norm(ax.Position - [0.1300 0.1100 0.7750 0.8150])<0.01;
        else
            need_new_figure = 0;
            axis(parent);
        end

        % cmap = colormap(cmap);
        if mean(cmap(1,:))<0.5
            cmap(1,:) = [0,0,0];
        end

        if isempty(cross_color)
            linecolor = 1 - cmap(1,:);
        else
            linecolor = cross_color;
        end

        if need_new_figure
            close(gcf);

            if ndims(img_mtg) == 3
                img_mtg = double(img_mtg)/255; %rgb
            end

            showc(flipud(img_mtg))
            axis xy
            colormapa(cmap);
        else
            img_mtg = flipud(img_mtg);
            switch ctype
                case 'index'
                    img_h = showc(img_mtg);
                case 'rgb'
                    mtg_max = max(img_mtg(:));
                    if mtg_max>255
                        img_mtg = uint16(img_mtg);
                    elseif mtg_max<=1
                        img_mtg = im2uint8(img_mtg);
                    else
                        img_mtg = uint8(img_mtg);
                    end

                    if char_cmap
                        try
                            cmap = imresize_grid(cmap,double([1+max(img_mtg(:),[],'omitnan'),3]));
                        catch
                            cmap = gray;
                        end
                    end

                    if ndims(img_mtg) == 2
                        img_rgb = ind2rgb(img_mtg,cmap);
                    else
                        img_rgb = double(img_mtg)/255; %rgb
                    end

                    img_h = showc(img_rgb);
            end

            img_h.Interpolation = 'bilinear';
            axis xy
            axis off
            ax.Colormap = cmap;
        end
        %%
        %pos_re = pos.*voxsize_norm([2,1,3]);

        % plot box
        hold on
        plot([0 1+imgsize_re(2)+imgsize_re(1)],imgsize_re(1)*[1 1],'Color',linecolor,'LineWidth',2);
        hold on
        plot(imgsize_re(2)*[1 1],[0 1+imgsize_re(1)+imgsize_re(3)],'Color',linecolor,'LineWidth',2);

        % plot cross line
        if ischar(cross_gap)
            cross_gap = str2double(cross_gap).*min(size(img_mtg));
        end
        add_gap = @(x,gap) [x-gap,NaN,x+gap];
        if cross_width>0
            hold on
            xx = [1 add_gap(pos_re(2),cross_gap) imgsize_re(2)];
            plot(xx,pos_re(1)*(xx*0+1),'Color',linecolor,'LineWidth',cross_width);
            hold on
            xx = [1 add_gap(pos_re(2),cross_gap), imgsize_re(2)+[add_gap(pos_re(1),cross_gap), imgsize_re(1)]];
            plot(xx,imgsize_re(1)+pos_re(3)*(xx*0+1),'Color',linecolor,'LineWidth',cross_width);
            hold on
            yy = [1 add_gap(pos_re(1),cross_gap) imgsize_re(1)+[add_gap(pos_re(3),cross_gap), imgsize_re(3)]];
            plot(pos_re(2)*(yy*0+1),yy,'Color',linecolor,'LineWidth',cross_width);
            hold on
            yy = [imgsize_re(1) imgsize_re(1)+[add_gap(pos_re(3),cross_gap), imgsize_re(3)]];
            plot((imgsize_re(2)+pos_re(1)*(yy*0+1)),yy,'Color',linecolor,'LineWidth',cross_width);
        end
        if cycle_width>0
            hold on
            viscircles([pos_re(2),pos_re(1); pos_re(2),pos_re(3)+imgsize_re(1); pos_re(1)+imgsize_re(3),pos_re(3)+imgsize_re(1)],...
                cycle_width/2*[1,1,1],'Color',linecolor);
        end

        % plot contur / deform grid
        if if_contour
            if 1 %(if_contour == 1)
                switch chs_max
                    case 1
                        label_list = {[1,2,3]};
                    case 3
                        label_list = {[1,3],[2,3],[1,2]};
                end
            end
            for i = 1:chs_max
                img = img_mtg(:,:,i);
                img(~ismember(img_labels, label_list{i})) = NaN;
                hold on; contour(img,if_contour,'Color',linecolor)
            end
        end

        if only_1dim
            switch plot_dim
                case 'x'
                    xLim = imgsize_re(2)+[1,imgsize_re(1)]; yLim = imgsize_re(1)+[1,imgsize_re(3)];
                    cen_cur = [imgsize_re(2)+pos_re(1),imgsize_re(1)+pos_re(3)];
                case 'y'
                    xLim = [1,imgsize_re(2)]; yLim = imgsize_re(1)+[1,imgsize_re(3)];
                    cen_cur = [pos_re(2),imgsize_re(1)+pos_re(3)];
                case 'z'
                    xLim = [1,imgsize_re(2)]; yLim = [1,imgsize_re(1)];
                    cen_cur = [pos_re(1),pos_re(2)];
            end

            if zoom>1
                xLim = cen_cur(1) + diff(xLim)/(2*zoom)*[-1,1];
                yLim = cen_cur(2) + diff(yLim)/(2*zoom)*[-1,1];
            end

            xlim(xLim); ylim(yLim)
        end

        os_h = flipud(img_mtg);
end
if strcmp(cAxis,'cmap') && ~char_cmap
    caxis([0,size(cmap,1)-1]);
end

fig = gcf;
if isempty(tips_info)==0
    dcm = datacursormode;
    %     dcm.Enable = 'on';
    dcm.UpdateFcn = @display_tips;
    fig.UserData.tips_info = tips_info;
    dcm.Enable = 'off';
end

if nargout>0
    varargout = {os_h};
end
end

function disp_info(src,evt,moving_func)
evname = evt.EventName;
switch(evname)
    case{'CrosshairMoved'}
        os_h = src;
        all_ui = os_h.Parent.Children;

        region_id_h = findobj(all_ui,'Tag','region_id');
        region_tag_h = findobj(all_ui,'Tag','region_tag');
        region_fullname_h = findobj(all_ui,'Tag','region_fullname');

        [hXY hYZ hXZ] = getAxesHandles(os_h);
        img_h = findobj(hXY.Children,'Type','Image');
        img = img_h.CData;
        xyz = evt.CurrentPosition;
        id = img(xyz(2),xyz(1));
        region_id_h.String = num2str(id);

        fig = gcf;
        tips_info = fig.UserData.tips_info;
        curosr_info = tips_info.(['id_',num2str(id)]);
        region_tag_h.String = curosr_info.datatip;
        region_fullname_h.String = curosr_info.datainfo;

        fig.UserData.ortho.locs(end+1,:) = xyz;

        moving_func(src,evt);
    case{'CrosshairMoving'}
        xyz = evt.CurrentPosition;
        fig = gcf;
        fig.UserData.ortho.locs(end+1,:) = xyz;
end
end

function moving_log(src,evt)
xyz = evt.CurrentPosition;
fig = gcf;
fig.UserData.ortho.locs(end+1,:) = xyz;
end

function txt = display_tips(~,info)
x = info.Position(1);
y = info.Position(2);
fig = gcf;
tips_info = fig.UserData.tips_info;
img = info.Target.CData;
cursor_idx = img(y,x);
curosr_info = tips_info.(['id_',num2str(cursor_idx)]);
txt = {['(' num2str(x) ', ' num2str(y) '): id = ',num2str(cursor_idx)];[curosr_info.datatip,': ',curosr_info.datainfo]};
end

function figScroll(src,evt)
mov = sign(evt.VerticalScrollCount);
ax = gca;
mov_idx = [0,0,0];
switch ax.Tag
    case 'XZ'
        mov_idx(2) = mov;
    case 'YZ'
        mov_idx(1) = mov;
    case 'XY'
        mov_idx(3) = mov;
    otherwise
        return
end

fig = gcf;
try
    os_h = fig.UserData.orthoslice;
    os_h.SliceNumbers = os_h.SliceNumbers + mov_idx;
    % fig.UserData.ortho.locs(end+1,:) = os_h.SliceNumbers;

    if isfield(fig.UserData,'moving_func')
        fig.UserData.moving_func(src,os_h.SliceNumbers);
    end
end
end

function clear_locLog(varargin)
fig = gcf;
fig.UserData.ortho.locs = [];
end

function plot_locLog(varargin)
fig = gcf;
pos = fig.UserData.ortho.locs;
mri = fig.UserData.ortho.vol;
voxsize = fig.UserData.ortho.voxsize;
iso_thresh = prctile(mri.data(1:3:end),30);

pos = pos.*voxsize([2,1,3]);
figure; next; plot_isosurf(mri,iso_thresh);
hold on; plot_line(pos)
end

function plot_mri_win(mri,info)
temp_dir = [pwd,'\temp_mat\'];
mkdir_silent(temp_dir)

if ~iscell(mri)
    mri = {mri};
end

for i = 1:length(mri)
    nii_file = temp_dir+"temp"+i+".nii.gz";
    nii = mri{i};
    if isstruct(nii)
        if isempty(info)
            info = nii;
        end
        nii = nii.data;
    end
    nii(round(end/2),round(end/2),round(end/2)) = 0;
    MRI_save_File(single(nii),nii_file,'info',info);
    winopen(nii_file);
    % val_cen = nii(round(end/2),round(end/2),round(end/2));
end

% pause(2)
% val_cen = "0.00000";
% win_title = "0x0x0= "+val_cen+"  ";
% [coords, titles, HWNDs] = listwindows(char(win_title));
% for i = 1:length(MRIs)
%     setwinpos(HWNDs(i), int32([500*i,200,500, 500]))
% end
end

function moving_proj(src,evt)
fig = gcf;
os_h = fig.UserData.orthoslice;
[hXY hYZ hXZ] = getAxesHandles(os_h);
mri_max = fig.UserData.ortho.vol_max;

if ~isnumeric(evt)
    xyz = evt.CurrentPosition;
else
    xyz = evt;
end

c_base = 0.8;

hXY.Children(end).CData = hXY.UserData.proj_imgs(:,:,xyz(3));
hXY.Children(end).AlphaData = c_base + 10*hXY.Children(end).CData/mri_max;

hYZ.Children(end).CData = squeeze(hYZ.UserData.proj_imgs(:,xyz(1),:));
hYZ.Children(end).AlphaData = c_base + 10*hYZ.Children(end).CData/mri_max;

hXZ.Children(end).CData = squeeze(hXZ.UserData.proj_imgs(xyz(2),:,:))';
hXZ.Children(end).AlphaData = c_base + 10*hXZ.Children(end).CData/mri_max;
end