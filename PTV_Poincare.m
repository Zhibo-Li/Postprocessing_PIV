%% Poincaré plot or recurrence plot (RP) for tracers from PTV experiments.

clear; close all; clc
set(0, 'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

magni = 0.1; % the magnification of the objective (um/pixel)
Ctr2Ctr = 30; % the pillar center-to-center distance (um)
Ctr2Ctr_pixel = Ctr2Ctr/magni; % the pillar center-to-center distance (pixel)

prompt = {'Enter tilt angle:'};
dlgtitle = 'Input';
fieldsize = [1 45];
definput = {'10'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
Array_angle = str2double(answer{1, 1});
RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);

mother_save_path = 'Z:\Processing & Results\PIV & PTV\Figures\';

exp2proc = uigetdir('Z:\Processing & Results\PIV & PTV\20240906-PTV\', 'Choose the folder of experiment:');
[exp2proc_path, exp2proc_case, ~] = fileparts(exp2proc);
[~, exp2proc_folderName, ~] = fileparts(exp2proc_path);

% % if ~exist(fullfile(mother_save_path, exp2proc_folderName), 'dir')
% %     mkdir(fullfile(mother_save_path, exp2proc_folderName));
% % end

[maskfilename, maskpathname] = uigetfile([exp2proc_path, '\.mat'], 'Choose the mask:');

load(fullfile(maskpathname, maskfilename));
% viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal
r_pillar = mean(radii)/Ctr2Ctr_pixel;



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
% % % viscircles(centers_corrected, radii,'LineStyle','--','LineWidth', 0.5, 'Color', 'r'); axis equal; hold on
%%%%%% About the coordinates rotation (end) %%%%%%



%%%%%% get the approximate (after-gridlization) pillar center positions (start) %%%%%%
sortxx = sort(centers_corrected(:, 1)); sortyy = sort(centers_corrected(:, 2));
sortxx_ind = find(diff(sortxx) > 100); sortyy_ind = find(diff(sortyy) > 100);
PAs_X = zeros(length(sortxx_ind)+1, 1); PAs_Y = zeros(length(sortyy_ind)+1, 1);
for jj = 1: length(sortxx_ind)-1
    PAs_X(jj+1) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
end
for jj = 1: length(sortyy_ind)-1
    PAs_Y(jj+1) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
end
PAs_X(1) = PAs_X(2) - mean(diff(PAs_X(2:end-1)));
PAs_X(end) = PAs_X(end-1) + mean(diff(PAs_X(2:end-1)));
PAs_Y(1) = PAs_Y(2) - mean(diff(PAs_Y(2:end-1)));
PAs_Y(end) = PAs_Y(end-1) + mean(diff(PAs_Y(2:end-1)));
% [XX, YY] = meshgrid(PAs_X, PAs_Y); plot(XX(:), YY(:))
%%%%%% get the approximate (after-gridlization) pillar center positions (end) %%%%%%

% average gap along y-direction
ave_y_gap = mean(diff(PAs_Y));



% % % % % %%
% % % % % % load *.mat file of streamlines from uPIV results
% % % % % load('Z:\Processing & Results\PIV & PTV\20240822-uPIV-streamlines\63x_Tra0o3um_Chip20211007_FlAng00d_10nL_Z=mid.mat');
% % % % % 
% % % % % n = 1; % counting for different trajectories
% % % % % for ii = 10: length(xy_stream)-10 % There are some empty trajectories
% % % % %     % read the data
% % % % %     traj_xy_i = xy_stream{ii}/magni; % already in x-y coordinate (NOT image coordinate)
% % % % % 
% % % % %     % calculate the confinement ratio: net distance / total distance travelled
% % % % %     persistence = sqrt(sum((traj_xy_i(1, :) - traj_xy_i(end, :)).^2)) / sum(sqrt(sum(diff(traj_xy_i).^2, 2)));
% % % % %     % skip the trajectory if it is too short or its confinement ratio is small
% % % % %     if size(traj_xy_i, 1) < 10 || persistence < 0.85
% % % % %         continue
% % % % %     end
% % % % % 
% % % % %     % remove NaNs in the traj_xy_i
% % % % %     traj_xy_i(isnan(traj_xy_i(:, 1)), :) = [];
% % % % %     traj_xy_i = (RotMatrix_correct * traj_xy_i')';
% % % % %     plot(traj_xy_i(:, 1), traj_xy_i(:, 2)); hold on
% % % % % % % %     viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal
% % % % % 
% % % % %     % divide the trajectory into pieces according to their y-position
% % % % %     fiber_Y_indicator = traj_xy_i(:, 2);
% % % % %     for kk = 1: length(PAs_Y)-1
% % % % %         fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
% % % % %     end
% % % % % 
% % % % %     % the 'entering lattice' positions
% % % % %     ind_ToBeMoved = min(abs(repmat(traj_xy_i(:, 1), 1, length(PAs_X)) - PAs_X'), [], 2) > 60;
% % % % %     traj_xy_i(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
% % % % %     fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
% % % % %     fiber_X = traj_xy_i(:, 1);
% % % % %     For_Poincare = nan(length(PAs_X), 3);
% % % % %     for kk = 1: length(PAs_X)
% % % % %         to_be_fitted2 = traj_xy_i(abs(fiber_X-PAs_X(kk))<=60, :);
% % % % %         if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
% % % % %             fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
% % % % %             For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;
% % % % % 
% % % % %             For_Poincare(kk, 2) = kk;
% % % % %             For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==traj_xy_i(:,1));
% % % % %         end
% % % % %     end
% % % % %     For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;
% % % % % 
% % % % %     Info.map{n} = For_Poincare;
% % % % %     Info.FlowAngle(n) = Array_angle;
% % % % %     n = n + 1;
% % % % % 
% % % % % end
%%



% load all the *. csv files in folder 'exp2proc' and process them them one by one in a for loop
csv_files = dir(fullfile(exp2proc, '*.csv'));
n = 1; % counting for different trajectories
for file = csv_files'
    % read the data
    data = readmatrix(fullfile(file.folder, file.name)); 
    traj_ID_xy = data(2:end, [3 5 6 9]); % ID, x, y, Frame_no

    % find the unique trajectory IDs to have different trajectories
    traj_ID = unique(traj_ID_xy(:, 1));
    for i = 1: length(traj_ID)
        traj_ID_i = traj_ID(i);
        traj_xy_i = traj_ID_xy(traj_ID_xy(:, 1) == traj_ID_i, 2:4);
        % sort traj_xy_i based on the first column (x)
        traj_xy_i = sortrows(traj_xy_i, 3); traj_xy_i(:,3) = [];
        % flip to x-y coordinate (NOT image coordinate)
        traj_xy_i(:, 2) = size(circleMask, 1) - traj_xy_i(:, 2);
%         plot(traj_xy_i(:,1), traj_xy_i(:,2)); axis equal; hold on

        % calculate the confinement ratio: net distance / total distance travelled
        persistence = sqrt(sum((traj_xy_i(1, :) - traj_xy_i(end, :)).^2)) / sum(sqrt(sum(diff(traj_xy_i).^2, 2)));  
        % skip the trajectory if it is too short or its confinement ratio is small
        if size(traj_xy_i, 1) < 10 || persistence < 0.85
            continue
        end

        
        traj_xy_i = (RotMatrix_correct * traj_xy_i')';
% % %         plot(traj_xy_i(:, 1), traj_xy_i(:, 2)); hold on
% % %         viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal

        % divide the trajectory into pieces according to their y-position
        fiber_Y_indicator = traj_xy_i(:, 2);
        for kk = 1: length(PAs_Y)-1
            fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
        end

        % the 'entering lattice' positions
        ind_ToBeMoved = min(abs(repmat(traj_xy_i(:, 1), 1, length(PAs_X)) - PAs_X'), [], 2) > 60;
        traj_xy_i(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
        fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
        fiber_X = traj_xy_i(:, 1);
        For_Poincare = nan(length(PAs_X), 3);
        for kk = 1: length(PAs_X)
            to_be_fitted2 = traj_xy_i(abs(fiber_X-PAs_X(kk))<=60, :);
            if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
                [fit_linear2, gof] = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
                % remove bad fitting
% % %                 if gof.rsquare < 0.9
% % %                     continue
% % %                 end
                For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;

                For_Poincare(kk, 2) = kk;
                For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==traj_xy_i(:,1));
            end
        end
        For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;

        Info.map{n} = For_Poincare;
        Info.FlowAngle(n) = Array_angle;
        n = n + 1;
        
% % %         % Temporary figure
% % %         Lattice_in = For_Poincare;
% % % 
% % %         Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
% % %         out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);
% % % 
% % %         Lattice_in = Lattice_in(out_in_diff==0)';
% % %         Lattice_out = Lattice_out(out_in_diff==0)';
% % % 
% % %         figure
% % %         plot(Lattice_in, Lattice_out, 'k.', 'LineStyle', 'none','MarkerSize', 8); 
% % %         hold on; axis equal
% % %         xlim([0 1]); ylim([0 1]);
    end
end



%% Figures

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];

%         figure('color', 'w'); set(gcf, 'Position', [100 100 500 500]);
%         plot(Lattice_in, Lattice_out, '.', 'LineStyle', 'none','MarkerSize', 13);
%         xlim([0 1]); ylim([0 1]); ax=gca; ax.FontSize = 15;
%         Info.name{jj}

    end

    f = figuresetting('centimeters',10,10,'times new roman',20,'off',2,'off','off');
    f.figure('');
    plot(Lattice_in_all, Lattice_out_all, 'kx', 'LineStyle', 'none','MarkerSize', 5);
    cmocean('amp'); caxis([0 50]);
%     hcb=colorbar; 
%     set(hcb,'TickLabelInterpreter','latex','Fontsize',20);
%     title(hcb,'$L(\mu{m})$','FontSize', 20,'Interpreter', 'latex');
    title_txt = strcat('$\alpha=', num2str(C(ii)), '^{\circ}$');

    f.interp_font('latex')
    f.axes('linear',[0 1],'linear',[0 1],'$y_{\mathrm{in}}/\lambda$','$y_{\mathrm{out}}/\lambda$',20);
    f.axes_ticks([0:0.2:1], [0:0.2:1]);
    grid on
    text(0.05, 0.9, title_txt,'FontSize', 20, ...
        'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])

    hold on
    % plot for tracer
    Simu_deg = C(ii);
    Simu_data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_' ...
        '0811exp_45deg\',num2str(Simu_deg),'deg\Data\Streamline_forPoincare_moreLines.csv']);
    XXYY_Simu = Simu_data(1:end, 13:14);  % (x, y) of the streamline
    IntegrationTime = Simu_data(1:end, 4);  % use to separate the different streamlines

    locs = find(IntegrationTime==0);
    ave_ele = size(XXYY_Simu, 1) / (length(locs)-1);  % average element number of each trajectory

    RotMatrix = rotz(-Simu_deg); RotMatrix = RotMatrix(1:2, 1:2);
    XXYY_Simu = (RotMatrix * XXYY_Simu')';  % after rotation

    PAs_X = (-30:3:30)*1e-4; PAs_Y = (-30:3:30)*1e-4;  % pillar positions
    ave_y_gap = 3e-4;

    for streamline_i = round((length(locs)-150)/2):3:length(locs)-round((length(locs)-150)/2)

        XXYY_Simu_i = XXYY_Simu(locs(streamline_i):locs(streamline_i+1)-1, :);

        XXYY_Simu_i(XXYY_Simu_i(:, 1) > 0.0015, :) = [];
        XXYY_Simu_i(XXYY_Simu_i(:, 1) < -0.0015, :) = [];
        fiber_Y_indicator = XXYY_Simu_i(:, 2);
        for kk = 1: length(PAs_Y)-1
            fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
        end

        % the 'entering lattice' positions
        ind_ToBeMoved = min(abs(repmat(XXYY_Simu_i(:, 1), 1, length(PAs_X)) - PAs_X), [], 2) > 6e-6;
        XXYY_Simu_i(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
        fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
        fiber_X = XXYY_Simu_i(:, 1);
        For_Poincare = nan(length(PAs_X), 3);
        for kk = 1: length(PAs_X)
            to_be_fitted2 = XXYY_Simu_i(abs(fiber_X-PAs_X(kk))<=6e-6, :);
            if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
                fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
                For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;
                % Correct: add '-PAs_Y(1)' @ 20230217
                For_Poincare(kk, 2) = kk;
                For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1));
                % Correct: change 'fiber_Y_indicator(kk)' to fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1))
            end
        end

        For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;
        Lattice_in = For_Poincare;
        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        plot(Lattice_in, Lattice_out, 'cx', 'LineStyle', 'none','MarkerSize', 3); hold on

    end

    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',fullfile(mother_save_path, ...
        exp2proc_folderName, [exp2proc_case, 'Poincare.eps']));

    exportgraphics(gcf,fullfile(mother_save_path, exp2proc_folderName, ...
        [exp2proc_case, 'Poincare.png']),'Resolution',100)

end



