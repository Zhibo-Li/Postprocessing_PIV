% Post processing for uPIV experiments in the pillar arrays
%
% First part: plot velocity field and flowType parameter in the whole
% channel at a given height (Z-scanning). Need to choose the masks of
% pîllar arrays in *.mat format.
%
% Second part: plot in a lattice (or unit cell). 
%
% NOTICE 1: need to be very careful about the flow direction in the
% experiments, especially when plot in a lattice!
%
% *NOTICE 2: particle images are transposed if the PIVlab calculations are 
% based on *.im7 format!!!



%% Velocity fields in whole channel in pillar array

clear; close all; clc

magni = 0.067; % the magnification of the objective (um/pixel)
PIV_step_size = 8; % the step length of uPIV calculation (in pixel)

mother_save_path = 'Z:\Processing & Results\PIV & PTV\Figures\';

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

exp2proc = uigetdir('Z:\Processing & Results\PIV & PTV\', 'Choose the mother folder of experiments:');
[~, exp2proc_folderName, ~] = fileparts(exp2proc);

if ~exist(fullfile(mother_save_path, exp2proc_folderName), 'dir')
    mkdir(fullfile(mother_save_path, exp2proc_folderName));
end

Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.

layersNUM = length(Selected_Names);
Z = nan(1, layersNUM);  
need_new_mask = 1;
for ii = 1:layersNUM

    %%%%%% About the mask loading (start) %%%%%%
    if need_new_mask == 1
        
        [~, theFOLDER, ~] = fileparts(Selected_Names{ii}); % NB: should not have '.' in the folder name.
        disp(theFOLDER); 

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


    try
        Z(ii) = str2double(extractBetween(theFOLDER,'Z','um'));
    catch
        Z(ii) = NaN; 
    end
    theFILE = dir(fullfile(exp2proc, theFOLDER, '*.mat'));
    theFILE = theFILE.name;
    deltaT = str2double(extractBetween(theFILE,'dt','us'));

    load(fullfile(exp2proc,theFOLDER,theFILE)); % load the ensemblePIV results

    vx_phy = v_filt{1,1}' * magni / deltaT; % velocity in m/s.
    %   vx_phy = -(-v_filt{1,1})' * magni / deltaT;        % Here is the reason.
    vy_phy = -u_filt{1,1}' * magni / deltaT; % velocity in m/s.
    x_phy = y{1,1}' * magni; % velocity in um.    
    %   x_phy = (size(circleMask, 1) - (size(circleMask, 1) - y{1,1}) )' * magni;         % Here is the reason.
    y_phy = (size(circleMask, 1) - x{1,1})' * magni; % velocity in um.
    % Consider NOTICE 2 and the fact that y{1,1} and v_filt{1,1} are in
    % image coordinates when they are calculated. 

    % to mask out the pillars
    velocityMask = ~circleMask(PIV_step_size:PIV_step_size:end-PIV_step_size, ...
        PIV_step_size:PIV_step_size:end-PIV_step_size);
    velocityMask_double = double(velocityMask);
    velocityMask_double(~velocityMask) = nan;
    vx_phy = vx_phy .* velocityMask_double; vy_phy = vy_phy .* velocityMask_double;

    figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
    contourf(x_phy, y_phy, sqrt(vx_phy.^2 + vy_phy.^2),100,'LineStyle','none'); hold on 
    shading interp; axis equal;
    quiver(x_phy, y_phy, vx_phy, vy_phy, 'Color','k');
%     xlabel('$x\ (\mathrm{\mu m})$');
%     ylabel('$y\ (\mathrm{\mu m})$');
    axis off
    cmocean('speed');

%     f=gcf;
%     exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, [theFOLDER, '_flowfield.png']),'Resolution',100);


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
%     c = colorbar;
%     c.Label.String = 'FlowType';
%     c.TickLabelInterpreter = 'LaTeX';
%     c.FontSize = 18;

    f=gcf;
    exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, [theFOLDER, '_flowType.png']),'Resolution',100);

end





%% Velocity fields within one lattice

clear; close all; clc
set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

magni = 0.067; % the magnification of the objective (um/pixel)
PIV_step_size = 8; % the step length of uPIV calculation (in pixel)
Ctr2Ctr = 30; % the pillar center-to-center distance (um)
Ctr2Ctr_pixel = Ctr2Ctr/magni; % the pillar center-to-center distance (pixel)

Plot_resol = 0.02; % plotting resolution in the lattice
% time_interval = 1/15*(500/20); % time interval for sampling (s) [NB: needed only for time fluctuation analysis]

mother_save_path = 'Z:\Processing & Results\PIV & PTV\Figures\';

prompt = {'Enter tilt angle:'};
dlgtitle = 'Input';
fieldsize = [1 45];
definput = {'15'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
Array_angle = str2double(answer{1, 1});
RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);

answer_plot = questdlg('Would you like to center the pillar in lattice?', ...
    'About plotting', ...
    'Yes','No','No');

exp2proc = uigetdir('Z:\Processing & Results\PIV & PTV\', 'Choose the mother folder of experiments:');
[~, exp2proc_folderName, ~] = fileparts(exp2proc);

if ~exist(fullfile(mother_save_path, exp2proc_folderName), 'dir')
    mkdir(fullfile(mother_save_path, exp2proc_folderName));
end

Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.

layersNUM = length(Selected_Names);
Z = nan(1, layersNUM); 
vx_z_phy_mean = zeros(1, layersNUM); vy_z_phy_mean = zeros(1, layersNUM); 
need_new_mask = 1;
for ii = 1:layersNUM

    %%%%%% About the mask loading (start) %%%%%%
    if need_new_mask == 1

        [~, theFOLDER, ~] = fileparts(Selected_Names{ii}); % NB: should NOT have '.' in the folder name.
        disp(theFOLDER);

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
    r_pillar = mean(radii)/Ctr2Ctr_pixel; 
    %%%%%% About the mask loading (end) %%%%%%



    %%%%%% About the coordinates rotation (start) %%%%%%
    centers_flip = centers;
    centers_flip(:, 2) = size(circleMask, 1) - centers_flip(:, 2);  % flip to x-y coordinate (NOT image coordinate)
    centers_new = (RotMatrix_init * centers_flip')'; % rotate based on the design
    % viscircles(centers_new, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

    % correct the rotation degree
    sorted_centers_new = sortrows(centers_new,2);
    sorty = sort(centers_new(:, 2));
    sorty_ind = find(diff(sorty) > 100); % Notice: change here!
    corrected_angle = zeros(length(sorty_ind)-1, 1);
    for jj = 1: length(sorty_ind)-1
        to_be_fitted = sorted_centers_new(sorty_ind(jj)+1:sorty_ind(jj+1), :);
        fit_linear = fit(to_be_fitted(:, 1), to_be_fitted(:, 2), 'poly1');
        k = fit_linear.p1;
        corrected_angle(jj) = atand(k);
    end
    corrected_angle = mean(corrected_angle);
    RotMatrix_correct = rotz(-Array_angle-corrected_angle);
    RotMatrix_correct = RotMatrix_correct(1:2, 1:2);
    centers_corrected = (RotMatrix_correct * centers_flip')';
    % viscircles(centers_corrected, radii,'LineStyle','--','LineWidth', 0.5, 'Color', 'r'); axis equal; hold on
    %%%%%% About the coordinates rotation (end) %%%%%%



    %%%%%% About the data loading (start) %%%%%%
    try
        Z(ii) = str2double(extractBetween(theFOLDER,'Z','um'));
    catch
        Z(ii) = NaN; 
    end
    theFILE = dir(fullfile(exp2proc, theFOLDER, '*.mat'));
    theFILE = theFILE.name;
    deltaT = str2double(extractBetween(theFILE,'dt','us'));

    load(fullfile(exp2proc,theFOLDER,theFILE)); % load the ensemblePIV results

    theFOLDER_no = str2double(extract(theFOLDER(end-2:end), digitsPattern));

    vx_phy = v_filt{1,1}' * magni / deltaT; % velocity in m/s.
    % vx_phy = -(-v_filt{1,1})' * magni / deltaT;        % Here is the reason.
    vy_phy = -u_filt{1,1}' * magni / deltaT; % velocity in m/s.
    x_phy = y{1,1}' * magni; % velocity in um.    
    %  x_phy = (size(circleMask, 1) - (size(circleMask, 1) - y{1,1}) )' * magni;         % Here is the reason.
    y_phy = (size(circleMask, 1) - x{1,1})' * magni; % velocity in um.
    % Consider NOTICE 2 and the fact that y{1,1} and v_filt{1,1} are in
    % image coordinates when they are calculated. 
    x_pixel = y{1,1}'; % x-coordinates in pixel
    y_pixel = (size(circleMask, 1) - x{1,1})'; % y-coordinates in pixel

    % to mask out the pillars
    velocityMask = ~circleMask(PIV_step_size:PIV_step_size:end-PIV_step_size, ...
        PIV_step_size:PIV_step_size:end-PIV_step_size);
    velocityMask_double = double(velocityMask);
    velocityMask_double(~velocityMask) = nan;
    vx_phy = vx_phy .* velocityMask_double; vy_phy = vy_phy .* velocityMask_double;
    v_phy_mag = sqrt(vx_phy.^2 + vy_phy.^2);
    %%%%%% About the data loading (end) %%%%%%



    %%%%%% calculate flow type indicator (start) %%%%%%
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
    %%%%%% calculate flow type indicator (end) %%%%%%


    %%%%%% Rotate into new coordinate and grid the data (start) %%%%%%
    XY_rotated = (RotMatrix_correct * [x_pixel(:) y_pixel(:)]')';
%     figure; plot(XY_rotated(:,1), XY_rotated(:,2)); axis equal; hold on
%     plot(centers_corrected(:,1), centers_corrected(:,2),'ro'); hold on

    UV_rotated = (RotMatrix_correct * [vx_phy(:) vy_phy(:)]')';

    % NOTICE: vx and vy should be also converted!
    vx_phy_rotated = UV_rotated(:,1); vy_phy_rotated = UV_rotated(:,2); 
    v_phy_mag_in_col = v_phy_mag(:);
    flow_type_in_col = flow_type(:);

    U_rotated_grid = cell(1/Plot_resol+1, 1/Plot_resol+1);
    V_rotated_grid = cell(1/Plot_resol+1, 1/Plot_resol+1);
    v_phy_mag_grid = cell(1/Plot_resol+1, 1/Plot_resol+1);
    flow_type_grid = cell(1/Plot_resol+1, 1/Plot_resol+1);

    for jj = 1:size(v_phy_mag_grid, 1)
        for kk = 1:size(v_phy_mag_grid, 2)

            switch answer_plot
                case 'No'
                    lower_lim = (jj-3/2)*Plot_resol*Ctr2Ctr_pixel;
                    upper_lim = (jj-1/2)*Plot_resol*Ctr2Ctr_pixel;
                    left_lim = (kk-3/2)*Plot_resol*Ctr2Ctr_pixel;
                    right_lim = (kk-1/2)*Plot_resol*Ctr2Ctr_pixel;
                case 'Yes'
                    lower_lim = (jj-3/2)*Plot_resol*Ctr2Ctr_pixel-Ctr2Ctr_pixel/2;
                    upper_lim = (jj-1/2)*Plot_resol*Ctr2Ctr_pixel-Ctr2Ctr_pixel/2;
                    left_lim = (kk-3/2)*Plot_resol*Ctr2Ctr_pixel-Ctr2Ctr_pixel/2;
                    right_lim = (kk-1/2)*Plot_resol*Ctr2Ctr_pixel-Ctr2Ctr_pixel/2;
            end
            % 'relative' edges of each grid around the pillar
            % remove '+Ctr2Ctr_pixel/2' to get the lattice in between four pillars

%             plot(right_lim+centers_corrected(9,1), lower_lim+centers_corrected(9,2),'m*'); hold on
            % to show position of the lattice close to the 9th pillar!

            jj_star = size(v_phy_mag_grid, 1)-jj+1; % for convenience

            for pp = 1:size(centers_corrected,1)

                select_ind = find(XY_rotated(:,2)>lower_lim+centers_corrected(pp,2) & ...
                    XY_rotated(:,2)<upper_lim+centers_corrected(pp,2) & ...
                    XY_rotated(:,1)>left_lim+centers_corrected(pp,1) & ...
                    XY_rotated(:,1)<right_lim+centers_corrected(pp,1));

                U_rotated_grid{jj_star,kk} = [U_rotated_grid{jj_star,kk}; vx_phy_rotated(select_ind)];
                V_rotated_grid{jj_star,kk} = [V_rotated_grid{jj_star,kk}; vy_phy_rotated(select_ind)];
                v_phy_mag_grid{jj_star,kk} = [v_phy_mag_grid{jj_star,kk}; v_phy_mag_in_col(select_ind)];
                flow_type_grid{jj_star,kk} = [flow_type_grid{jj_star,kk}; flow_type_in_col(select_ind)];

            end

            [N,edges] = histcounts(U_rotated_grid{jj_star,kk},'Normalization','probability');
            U_rotated_grid{jj_star,kk} = N*movmean(edges',2,"Endpoints","discard"); % weighted average
            [N,edges] = histcounts(V_rotated_grid{jj_star,kk},'Normalization','probability');
            V_rotated_grid{jj_star,kk} = N*movmean(edges',2,"Endpoints","discard");
            [N,edges] = histcounts(v_phy_mag_grid{jj_star,kk},'Normalization','probability');
            v_phy_mag_grid{jj_star,kk} = N*movmean(edges',2,"Endpoints","discard");
            [N,edges] = histcounts(flow_type_grid{jj_star,kk},'Normalization','probability');
            flow_type_grid{jj_star,kk} = N*movmean(edges',2,"Endpoints","discard");

        end
    end

    [lattice_X, lattice_Y] = meshgrid(0:Plot_resol:1, 1:-Plot_resol:0); % NOTICE: the origin is at left-bottom corner.
    U_rotated_grid = cell2mat(U_rotated_grid);
    V_rotated_grid = cell2mat(V_rotated_grid);
    v_phy_mag_grid = cell2mat(v_phy_mag_grid);
    flow_type_grid = cell2mat(flow_type_grid);
    
    %%%%%% Rotate into new coordinate and grid the data (end) %%%%%%



    %%%%%%%%%%%%%%%%%% Figures (start) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch answer_plot
        case 'Yes'
            %%% Figures: velocity field in lattice (pillar centered)
            figure('color', 'w','units','normalized','outerposition',[0 0 1 1]);
            contourf(lattice_X, lattice_Y, v_phy_mag_grid/max(max(v_phy_mag_grid)), 100,'LineStyle','none'); hold on
            shading interp; axis equal; axis off
            quiver(lattice_X, lattice_Y, U_rotated_grid, V_rotated_grid, 'Color','k');
            cmocean('speed'); caxis([0 5]);
            viscircles([0.5, 0.5], r_pillar,'Color','b'); 
            camroll(180)

            f=gcf;
            exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, ...
                [theFOLDER, '_flowfield_in_lattice_ave.png']),'Resolution',100)

            %%% Figures: flowtype in lattice (pillar centered)
            figure('color', 'w','units','normalized','outerposition',[0 0 1 1]);
            contourf(lattice_X, lattice_Y, flow_type_grid, 100,'LineStyle','none');
            shading interp; axis equal; axis off
            cmocean('balance'); caxis([-1 1])
%             c = colorbar;
%             c.Label.String = 'FlowType';
%             c.TickLabelInterpreter = 'LaTeX';
%             c.FontSize = 18;
            viscircles([0.5, 0.5], r_pillar,'Color','b'); 
            camroll(180)

            f=gcf;
            exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, ...
                [theFOLDER, '_flowType_in_lattice_ave.png']),'Resolution',100)

        case 'No'
            %%% Figures: velocity field in lattice (in between four pillars)
            figure('color', 'w','units','normalized','outerposition',[0 0 1 1]);
            contourf(lattice_X, lattice_Y, v_phy_mag_grid/max(max(v_phy_mag_grid)), 100,'LineStyle','none'); hold on
            shading interp; axis equal; axis off
            quiver(lattice_X, lattice_Y, U_rotated_grid, V_rotated_grid, 'Color','k');
            cmocean('speed'); caxis([0 1]);
            viscircles([0, 0; 0, 1; 1, 0; 1, 1], r_pillar*ones(4,1),'Color','b');
            xlim([0 1]); ylim([0 1]); 
            camroll(180)

            f=gcf;
            exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, ...
                [theFOLDER, '_flowfield_in_lattice_ave.png']),'Resolution',100)
            

            %%% Figures: flowtype in lattice (in between four pillars)
            figure('color', 'w','units','normalized','outerposition',[0 0 1 1]);
            contourf(lattice_X, lattice_Y, flow_type_grid, 100,'LineStyle','none');
            shading interp; axis equal; axis off; 
            cmocean('balance'); caxis([-1 1])
%             c = colorbar;
%             c.Label.String = 'FlowType';
%             c.TickLabelInterpreter = 'LaTeX';
%             c.FontSize = 18;
            viscircles([0, 0; 0, 1; 1, 0; 1, 1], r_pillar*ones(4,1),'Color','b');
            xlim([0 1]); ylim([0 1]); 
            camroll(180)

            f=gcf;
            exportgraphics(f,fullfile(mother_save_path, exp2proc_folderName, ...
                [theFOLDER, '_flowType_in_lattice_ave.png']),'Resolution',100)

            close all
    end


% % %     U_rotated_grid_all(:,:,theFOLDER_no) = U_rotated_grid;
% % %     V_rotated_grid_all(:,:,theFOLDER_no) = V_rotated_grid;
% % %     v_phy_mag_grid_all(:,:,theFOLDER_no) = v_phy_mag_grid/max(max(v_phy_mag_grid));
% % %     flow_type_grid_all(:,:,theFOLDER_no) = flow_type_grid;

end

% % % v_mean_phy_mag_grid = mean(v_phy_mag_grid_all, 3);
% % % v_prime_phy_mag_grid_all = v_phy_mag_grid_all - repmat(v_mean_phy_mag_grid, [1 1 size(v_phy_mag_grid_all, 3)]);
% % % normalized_v_prime_phy_mag_grid_all = v_prime_phy_mag_grid_all ./ repmat(v_mean_phy_mag_grid, [1 1 size(v_phy_mag_grid_all, 3)]);
% % % 
% % % normalized_v_prime_time = squeeze(normalized_v_prime_phy_mag_grid_all(11,:,:)); 
% % % % Notice that the flow field in the lattice has been rotated!!!
% % % % So, the rows are along the vertical direction in the flow field figures.
% % % [lattice_time, lattice_YY] = meshgrid(0:time_interval:(size(normalized_v_prime_time,2)-1)*time_interval, 0:Plot_resol:1);
% % % 
% % % % Figures: velocity time fluctuations sampled on a line.
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 400]); % plot the v-profiles
% % % imagesc(lattice_time(1,:), lattice_YY(:,1)', normalized_v_prime_time); 
% % % set(gca, 'YTick', [], 'FontSize', 24,'TickLabelInterpreter','latex');  
% % % xlabel('$\rm{Time}\ (s)$');
% % % cmocean('balance'); caxis([-0.4 0.4]);
% % % c = colorbar;
% % % c.FontName = 'Times New Roman';
% % % c.Label.String = '$u^{\prime}/<u>_t$';
% % % c.Label.Interpreter = 'LaTeX';
% % % c.FontSize = 24;
% % % 
% % % set(gcf,'renderer','Painters');
% % % print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\PIV & PTV' ...
% % %     '\Figures\PEO\100nL_normalized_v_prime_time.eps']);
% % % 
% % % 
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 300]); % plot the v-profiles
% % % plot(lattice_time(1,:), normalized_v_prime_time(11,:),':k','LineWidth',2); hold on
% % % plot(lattice_time(1,:), normalized_v_prime_time(5,:),':m','LineWidth',2); hold on
% % % plot(lattice_time(1,:), normalized_v_prime_time(17,:),':r','LineWidth',2); hold on
% % % set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.3, 'FontSize', 24,'TickLabelInterpreter','latex')
% % % ylim([-0.3 0.3])
% % % ylabel('$u^{\prime}/<u>_t$','FontSize', 24,'Interpreter', 'latex');
% % % xlabel('$\rm{Time}\ (s)$','FontSize', 24,'Interpreter', 'latex');
% % % 
% % % set(gcf,'renderer','Painters');
% % % print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\PIV & PTV' ...
% % %     '\Figures\PEO\100nL_normalized_v_prime_time_on_one_point.eps']);