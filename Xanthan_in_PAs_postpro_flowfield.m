% % % % % % % % % % Plot velocity fields % % % % % % % % % %

clear; close all; clc

magni = 0.067; % the magnification of the objective (um/pixel)
PIV_step_size = 16; % the step length of uPIV calculation (in pixel)

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

exp2proc = uigetdir('Z:\Processing & Results\PIV & PTV\', 'Choose the mother folder of experiments:');

Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.

layersNUM = length(Selected_Names);
Z = zeros(1, layersNUM);  
vx_z_phy_mean = zeros(1, layersNUM); vy_z_phy_mean = zeros(1, layersNUM); 
need_new_mask = 1;
for ii = 1:layersNUM

    %%%%%% About the mask loading (start) %%%%%%
    if need_new_mask == 1
        disp(extractBetween(Selected_Names{ii},'Z','um'));

        [maskfilename, maskpathname] = uigetfile([exp2proc, '\.mat'], 'Choose the mask:');

        load(fullfile(maskpathname, maskfilename));
        % viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

        answer = questdlg('Would you like to input a different mask?', ...
            'About mask', ...
            'Yes','No','No');
        % Handle response
        switch answer
            case 'Yes'
                need_new_mask = 1;
            case 'No'
                need_new_mask = 0;
        end
    end
    %%%%%% About the mask loading (end) %%%%%%

    theFOLDER = Selected_Names{ii}(length(exp2proc)+2:end);
    Z(ii) = str2double(extractBetween(theFOLDER,'Z','um'));
    theFILE = dir(fullfile(exp2proc, theFOLDER, '*.mat'));
    theFILE = theFILE.name;
    deltaT = str2double(extractBetween(theFILE,'dt','us'));

    load(fullfile(exp2proc,theFOLDER,theFILE)); % load the ensemblePIV results

    vx_phy = u_filt{1,1} * magni / deltaT; % velocity in m/s.
    vy_phy = v_filt{1,1} * magni / deltaT; % velocity in m/s.
    x_phy = x{1,1} * magni; % velocity in um.
    y_phy = y{1,1} * magni; % velocity in um.

    % to mask out the pillars
    velocityMask = ~circleMask(PIV_step_size:PIV_step_size:end-PIV_step_size, ...
        PIV_step_size:PIV_step_size:end-PIV_step_size);
    velocityMask_double = double(velocityMask);
    velocityMask_double(~velocityMask) = nan;
    vx_phy = vx_phy .* velocityMask_double'; vy_phy = vy_phy .* velocityMask_double';

    figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
    contourf(x_phy, y_phy, sqrt(vx_phy.^2 + vy_phy.^2),100,'LineStyle','none'); hold on 
    shading interp; axis equal;
    quiver(x_phy, y_phy, vx_phy, vy_phy, 'Color','k');
%     xlabel('$x\ (\mathrm{\mu m})$');
%     ylabel('$y\ (\mathrm{\mu m})$');
    axis off
    cmocean('speed');

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\PIV & PTV\Figures\', exp2proc(end-30:end-22), theFOLDER(end-8:end), 'flowfield.png'],'Resolution',100)



    %%%%% calculate flow type indicator %%%%%
    dx = (x_phy(1,2) - x_phy(1,1))*1e-6; dy = (y_phy(2,1) - y_phy(1,1))*1e-6; 

    [du_dx,du_dy] = gradient(vx_phy); 
    [dv_dx,dv_dy] = gradient(vy_phy); % NOTICE the direction!!

    du_dx = du_dx / dx; dv_dx = dv_dx / dx;
    du_dy = du_dy / dy; dv_dy = dv_dy / dy;

    VGT_sym = [du_dx(:), 0.5*(du_dy(:)+dv_dx(:)), dv_dy(:), 0.5*(du_dy(:)+dv_dx(:))];
    VGT_anti = [zeros(numel(du_dx),1), 0.5*(dv_dx(:)-du_dy(:)), zeros(numel(du_dx),1), 0.5*(du_dy(:)-dv_dx(:))];

    % flow type indicator (lambda)
    flow_type = nan(numel(du_dx),1);
    for jj = 1: numel(du_dx)

        norm_E = norm(reshape(VGT_sym(jj, :), 2, 2),'fro');
        norm_O = norm(reshape(VGT_anti(jj, :), 2, 2),'fro');

        flow_type(jj) = (norm_E-norm_O) / (norm_E+norm_O);

    end
    flow_type = reshape(flow_type, size(x_phy, 1), size(y_phy, 2));

    figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
    contourf(x_phy, y_phy, flow_type,100,'LineStyle','none');
    shading interp; axis equal; axis off
    cmocean('balance'); caxis([-1 1])
    c = colorbar;
    c.Label.String = 'FlowType';
    c.TickLabelInterpreter = 'LaTeX';
    c.FontSize = 18;

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\PIV & PTV\Figures\', exp2proc(end-30:end-22), theFOLDER(end-8:end), 'flowType.png'],'Resolution',100)

end



