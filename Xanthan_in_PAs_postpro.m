% % % % % % % % % % Plot velocity profiles % % % % % % % % % %

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%% Blank area (along z-direction) %%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

magni = 0.067; % the magnification of the objective (um/pixel)
PIV_step_size = 16; % the step length of uPIV calculation (in pixel)

% EXP 1
exp2proc = 'Z:\Processing & Results\PIV & PTV\20240208_Xanthan100ppm_filtered';
n = VicFc_Xanthan_uPIV_fitting_velocity_profile_blank_channel(magni, exp2proc, 'r', 'o', 1, 1);
text(20, 110, ['n (100 ppm)=', num2str(n)],'FontName','times new roman','FontSize', 16)

% EXP 2
exp2proc = 'Z:\Processing & Results\PIV & PTV\20240206_Xanthan500ppm_filtered';
n = VicFc_Xanthan_uPIV_fitting_velocity_profile_blank_channel(magni, exp2proc, 'm', '*', 1, 0);
text(20, 85, ['n (500 ppm)=', num2str(n)],'FontName','times new roman','FontSize', 16)

% EXP 3
exp2proc = 'Z:\Processing & Results\PIV & PTV\20240220_Xanthan1000ppm_filtered';
n = VicFc_Xanthan_uPIV_fitting_velocity_profile_blank_channel(magni, exp2proc, 'g', '^', 1, 0);
text(20, 60, ['n (1000 ppm)=', num2str(n)],'FontName','times new roman','FontSize', 16)

% % % EXP 4 (water) 
% exp2proc = 'Z:\Processing & Results\PIV & PTV\20240220_Xanthan1000ppm_filtered';
% n = VicFc_Xanthan_uPIV_fitting_velocity_profile_blank_channel(magni, exp2proc, 'k', 's', 1, 0);
% text(20, 135, ['n (water)=', num2str(n)],'FontName','times new roman','FontSize', 16)

xlabel('$z\ (\mathrm{\mu m})$');
ylabel('$u_x\ (\mathrm{\mu m/s})$');
ylim([0 500])
legend('$\rm{100\,ppm}$','','$\rm{500\,ppm}$','','$\rm{1000\,ppm}$','fitting','location','northeast')

f=gcf;
exportgraphics(f,'Z:\Processing & Results\PIV & PTV\Figures\Xanthan_Ux-Z_20240206_2024208_20240220.png','Resolution',100)

% f=gcf;
% savefig(f,[exp2proc,'\Xanthan_Ux-Z.fig'])
% set(f,'renderer','Painters');
% print('-depsc2','-tiff','-r100','-painters',[exp2proc,'\Xanthan_Ux-Z.eps'])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%% Blank area (along y-direction) %%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

magni = 0.067; % the magnification of the objective (um/pixel)

exp2proc = uigetdir('Z:\Processing & Results\PIV & PTV\');

Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.

layersNUM = length(Selected_Names);
Z = zeros(1, layersNUM);  UUU_alongZ = zeros(1, layersNUM); VVV_alongZ = zeros(1, layersNUM);
UUU_alongZ_err = zeros(1, layersNUM); VVV_alongZ_err = zeros(1, layersNUM); 
for ii = 1:layersNUM

    theFOLDER = Selected_Names{ii}(length(exp2proc)+2:end);
    Z(ii) = str2double(extractBetween(theFOLDER,'Z','um'));
    theFILE = dir(fullfile(exp2proc, theFOLDER, '*.mat'));
    theFILE = theFILE.name;
    deltaT = str2double(extractBetween(theFILE,'dt','us'));

    load(fullfile(exp2proc,theFOLDER,theFILE)); % load the ensemblePIV results

    vx_phy = u_filt{1,1} * magni / deltaT; % velocity in m/s.
    vy_phy = v_filt{1,1} * magni / deltaT; % velocity in m/s.
    x_phy = x{1,1} * magni; % x-coordinates in um.
    y_phy = y{1,1} * magni; % y-coordinates in um.

    vx_y_phy_mean(:, ii) = mean(vx_phy, 2); % V_x (along y) of no.ii layer (real height is Z(ii))
    y_phy_mean(:, ii) = mean(y_phy, 2);

%     vx_z_phy_mean(ii) = mean(mean(vx_phy(round(size(vx_phy,1)/4),:))); % V_x (along z) of no.ii layer (real height is Z(ii))

end

counter = 1;
figure('color', 'w');  set_plot(gcf, gca)
for jj = 1:length(Z)
    plot(y_phy_mean(:, jj),  vx_y_phy_mean(:, jj)*1e6, 'LineWidth', 2, ...
        'LineStyle','none', 'Marker','o','MarkerSize', 6); hold on
    legend_txt{counter} = ['$Z=',num2str(Z(jj)),'\mu{m}$'];
    counter = counter + 1;
end
xlabel('$y\ (\mathrm{\mu m})$');
ylabel('$u_x\ (\mathrm{\mu m/s})$');
legend(legend_txt,'location','southeast','FontSize',12);

f=gcf;
exportgraphics(f,'Z:\Processing & Results\PIV & PTV\Figures\Xanthan_Ux-Y.png','Resolution',100)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%% Pillar array area %%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

prompt = {'Enter tilt angle:'};
dlgtitle = 'Input';
fieldsize = [1 45];
definput = {'30'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
Array_angle = str2double(answer{1, 1});
RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);

layersNUM = length(Selected_Names);
Z = zeros(1, layersNUM);  
vx_z_phy_mean = zeros(1, layersNUM); vy_z_phy_mean = zeros(1, layersNUM); 
need_new_mask = 1;
for ii = 1:layersNUM

    %%%%%% About the masks (start) %%%%%%
    if need_new_mask == 1
        disp(extractBetween(Selected_Names{ii},'Z','um'));

        [maskfilename, maskpathname] = uigetfile([exp2proc, '\.mat'], 'Choose the mask:');

        load(fullfile(maskpathname, maskfilename));
        % viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

        centers_flip = centers;
        centers_flip(:, 2) = size(circleMask, 1) - centers_flip(:, 2);  % flip to x-y coordinate (NOT image coordinate)
        centers_new = (RotMatrix_init * centers_flip')'; % rotate based on the design
        % viscircles(centers_new, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

        % correct the rotation degree
        sorted_centers_new = sortrows(centers_new,2);
        sorty = sort(centers_new(:, 2));
        sorty_ind = find(diff(sorty) > 100);
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

        % get the approximate (after-gridlization) pillar center positions
        sortxx = sort(centers_corrected(:, 1)); sortyy = sort(centers_corrected(:, 2));
        sortxx_ind = find(diff(sortxx) > 100); sortyy_ind = find(diff(sortyy) > 100);
        PAs_X = zeros(length(sortxx_ind)+3, 1); PAs_Y = zeros(length(sortyy_ind)+3, 1);
        for jj = 1: length(sortxx_ind)-1
            PAs_X(jj+2) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
        end
        for jj = 1: length(sortyy_ind)-1
            PAs_Y(jj+2) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
        end
        % make a grid of pillar CoMs (enlarge the PAs area)
        Ctr2Ctr_x = mean(diff(PAs_X(3:end-2))); Ctr2Ctr_y = mean(diff(PAs_Y(3:end-2)));
        PAs_X(2) = PAs_X(3) - Ctr2Ctr_x; PAs_X(1) = PAs_X(2) - Ctr2Ctr_x;
        PAs_X(end-1) = PAs_X(end-2) + Ctr2Ctr_x; PAs_X(end) = PAs_X(end-1) + Ctr2Ctr_x;
        PAs_Y(2) = PAs_Y(3) - Ctr2Ctr_y; PAs_Y(1) = PAs_Y(2) - Ctr2Ctr_y;
        PAs_Y(end-1) = PAs_Y(end-2) + Ctr2Ctr_y; PAs_Y(end) = PAs_Y(end-1) + Ctr2Ctr_y;
        [PAs_X_mesh,PAs_Y_mesh] = meshgrid(PAs_X,PAs_Y);
        % PAs_xy = complex(PAs_X_mesh(:), PAs_Y_mesh(:));
        PAs_xy = [PAs_X_mesh(:), PAs_Y_mesh(:)];
        % viscircles(PAs_xy, 100*ones(size(PAs_xy,1),1),'LineStyle','--','LineWidth', 0.5, 'Color', 'b'); axis equal; hold on

        % centers of a unit cell
        PAs_gap_X_mesh = movmean(PAs_X_mesh, 2, 2); PAs_gap_X_mesh = PAs_gap_X_mesh(2:end, 2:end);
        PAs_gap_Y_mesh = movmean(PAs_Y_mesh, 2, 1); PAs_gap_Y_mesh = PAs_gap_Y_mesh(2:end, 2:end);
        PAs_gap_xy = [PAs_gap_X_mesh(:), PAs_gap_Y_mesh(:)];
        % viscircles(PAs_gap_xy, 50*ones(size(PAs_gap_xy,1),1),'LineStyle','--','LineWidth', 0.5, 'Color', 'm'); axis equal; hold on

        % Rotate back
        PAs_gap_xy_origin = (RotMatrix_correct\PAs_gap_xy')'; % (inv(RotMatrix_correct) * PAs_gap_xy')'
        
        % NOTICE: need to adjust the coordinates again (to velocity field plotting coordinates)!!!
        PAs_gap_xy_origin_tmp = PAs_gap_xy_origin;
        PAs_gap_xy_origin(:, 2) = PAs_gap_xy_origin_tmp(:, 1); % Y = X
        PAs_gap_xy_origin(:, 1) = size(circleMask, 1) - PAs_gap_xy_origin_tmp(:, 2); % X = image_height - Y

        PAs_gap_xy_origin(PAs_gap_xy_origin(:,1)<0, :) = [];
        PAs_gap_xy_origin(PAs_gap_xy_origin(:,1)>size(circleMask, 1), :) = [];
        PAs_gap_xy_origin(PAs_gap_xy_origin(:,2)<0, :) = [];
        PAs_gap_xy_origin(PAs_gap_xy_origin(:,2)>size(circleMask, 2), :) = [];
        % viscircles(PAs_gap_xy_origin, 10*ones(size(PAs_gap_xy_origin,1),1),'LineStyle','--','LineWidth', 3, 'Color', 'c'); axis equal; hold on

        % normalize the crosspoints by the calculation gap length (for V-profile along Z).
        LocSelectX = round(PAs_gap_xy_origin(:,1)/PIV_step_size); 
        LocSelectY = round(PAs_gap_xy_origin(:,2)/PIV_step_size);
        crossnum = size(PAs_gap_xy_origin, 1);

        % normalize the X by the calculation gap length (for V-profile along flow direction).
        theNormX = round(PAs_gap_xy_origin(:,2)/PIV_step_size);

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
    %%%%%% About the masks (end) %%%%%%


    theFOLDER = Selected_Names{ii}(length(exp2proc)+2:end);
    Z(ii) = str2double(extractBetween(theFOLDER,'Z','um'));
    theFILE = dir(fullfile(exp2proc, theFOLDER, '*.mat'));
    theFILE = theFILE.name;
    deltaT = str2double(extractBetween(theFILE,'dt','us'));

    load(fullfile(exp2proc,theFOLDER,theFILE)); % load the ensemblePIV results

    vx_phy = u_filt{1,1} * magni / deltaT; % velocity in m/s.
    vy_phy = v_filt{1,1} * magni / deltaT; % velocity in m/s.
    x_phy = x{1,1} * magni; % x-coordinates in um.
    y_phy = y{1,1} * magni; % y-coordinates in um.

    % the U and V along the z-direction.
    UU = []; VV = []; % calculate the speed near the cross-points
    for pp = 1:crossnum
%         imshow(vx_phy, []); hold on; axis equal;
%         viscircles([LocSelectX(pp),LocSelectY(pp)], 1,'LineStyle','--','LineWidth', 3, 'Color', 'g');
        try
            UU = [UU, mean(mean(vx_phy(LocSelectY(pp)-2:LocSelectY(pp)+2,...
                LocSelectX(pp)-2:LocSelectX(pp)+2),'omitnan'))];
            VV = [VV, mean(mean(vy_phy(LocSelectY(pp)-2:LocSelectY(pp)+2,...
                LocSelectX(pp)-2:LocSelectX(pp)+2),'omitnan'))];
        catch
            warning('The intersect point is out of the field.');
            UU = UU; VV = VV;
        end
    end
    vx_z_phy_mean(ii) = mean(UU,'omitnan');  vy_z_phy_mean(ii) = mean(VV,'omitnan'); % the averaged the U, V over cross-points.


% %     % the U and V along the flow-direction.
% %     UU = []; VV = [];
% %     if 24<Z(ii) && Z(ii)<25   % Z = 25.07: middle plane.
% %         for pp = 1:length(theNormX)
% %             UU = vx_phy(theNormX(pp), :);
% %             VV = vy_phy(theNormX(pp), :);
% %             XX = x_phy(theNormX(pp), :);
% %             YY = y_phy(theNormX(pp), :);
% %         end
% %         vx_x_phy_mean = UU;  vy_x_phy_mean = VV;
% %         x_phy_mean = XX;  y_phy_mean = YY;
% %     end

end

norm_v = sqrt(vx_z_phy_mean.^2 + vy_z_phy_mean.^2);

figure('color', 'w');  set_plot(gcf, gca)
plot(Z, norm_v*1e6,'Color', [77,175,74]/255, 'LineWidth', 2, ...
    'LineStyle','none', 'Marker','o','MarkerSize', 10); hold on

% % % fitting: v(r) = v_max*(1-(|r|/halfChannel_H)^(1+1/n)) --> 'v_max', 'halfChannel_H', 'n'
% % ft = fittype( 'v_max * (1-(abs(x-halfChannel_H)/halfChannel_H)^(1+1/n))', 'independent', 'x', 'dependent', 'y', 'coefficients', {'v_max', 'halfChannel_H', 'n'});
% % ftopts = fitoptions( 'Method', 'NonlinearLeastSquares','StartPoint', [max(norm_v)*1e6 max(Z)/2 1]);
% % ftopts.Display = 'Off';
% % [fitobject, GOF, fitinfo]  = fit(Z', norm_v'*1e6, ft, ftopts);
% % X = 0:0.1:max(Z);
% % plot(X, feval(fitobject,X),'--k','linewidth',1.5)

xlabel('$z\ (\mathrm{\mu m})$');
ylabel('$|v|\ (\mathrm{\mu m/s})$');
f=gcf;
exportgraphics(f,'Z:\Processing & Results\PIV & PTV\Figures\Xanthan_PAs_normV-Z_20240208.png','Resolution',100)
% 
% figure('color', 'w');  set_plot(gcf, gca)
% plot(x_phy_mean, vx_x_phy_mean*1e6,'Color', [55,126,184]/255, 'LineWidth', 2, ...
%     'LineStyle','none', 'Marker','o','MarkerSize', 10); hold on
% xlabel('$x\ (\mathrm{\mu m})$');
% ylabel('$u_x\ (\mathrm{\mu m/s})$');
% f=gcf;
% exportgraphics(f,'Z:\Processing & Results\PIV & PTV\Figures\Xanthan_PAs_Ux-X.png','Resolution',100)
% 
% figure('color', 'w');  set_plot(gcf, gca)
% plot(x_phy_mean, vy_x_phy_mean*1e6,'Color', [228,26,28]/255, 'LineWidth', 2, ...
%     'LineStyle','none', 'Marker','o','MarkerSize', 10); hold on
% xlabel('$x\ (\mathrm{\mu m})$');
% ylabel('$u_y\ (\mathrm{\mu m/s})$');
% f=gcf;
% exportgraphics(f,'Z:\Processing & Results\PIV & PTV\Figures\Xanthan_PAs_Uy-X.png','Resolution',100)
