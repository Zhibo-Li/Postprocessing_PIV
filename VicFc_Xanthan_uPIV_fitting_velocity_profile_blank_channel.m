function fluid_index = VicFc_Xanthan_uPIV_fitting_velocity_profile_blank_channel(magni, exp2proc, linecolor, marker, if_scan_whole_channel, fig_new)
% Plot viscosity vs. shear rate curve of Xanthan suspension
%
% INPUT:
% magni: the magnification of the objective (um/pixel)
% exp2proc: dir which stores the data (uPIV result folders).
% linecolor: curve color (1*3 vector or str like 'r').
% marker: the type of the marker. e.g. 'o'.
% if_scan_whole_channel: if the data scans the whole channel height.
% fig_new: if plot a new figure. 0: no.
% 
% OUTPUT:
% fluid_index: fitted fluid index n based on power-law model

Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.

layersNUM = length(Selected_Names);
Z = zeros(1, layersNUM);  vx_z_phy_mean = zeros(1, layersNUM); 
for ii = 1:layersNUM

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

    vx_z_phy_mean(ii) = mean(mean(vx_phy(round(size(vx_phy,1)/4),:))); % V_x (along z) of no.ii layer (real height is Z(ii))

end

if nargin > 5 && fig_new ~= 0
    figure('color', 'w');  set_plot(gcf, gca);
end
plot(Z, abs(vx_z_phy_mean)*1e6,'Color', linecolor, 'LineWidth', 2, ...
    'LineStyle','none', 'Marker', marker,'MarkerSize', 10); hold on

% % % fitting: v(r) = v_max*(1-(|r|/halfChannel_H)^(1+1/n)) --> 'v_max', 'n'
% % ft = fittype( 'v_max * (1-x^(1+1/n))', 'independent', 'x', 'dependent', 'y', 'coefficients', {'v_max', 'n'});
% % ftopts = fitoptions( 'Method', 'NonlinearLeastSquares','StartPoint', [max(abs(vx_z_phy_mean))*1e6 1]);
% % ftopts.Display = 'Off'; 
% % [fitobject, GOF, fitinfo]  = fit(abs(2*Z/max(Z)-1)', abs(vx_z_phy_mean)'*1e6, ft, ftopts);
% % X = 0:0.1:max(Z);
% % plot(X, feval(fitobject,abs(2*X/max(X)-1)),'--k','linewidth',1.5)

if if_scan_whole_channel
    % fitting: v(r) = v_max*(1-(|r|/halfChannel_H)^(1+1/n)) --> 'v_max', 'halfChannel_H', 'n'
    ft = fittype( 'v_max * (1-(abs(x-x_shift-halfChannel_H)/halfChannel_H)^(1+1/n))', ...
        'independent', 'x', 'dependent', 'y', 'coefficients', {'v_max', 'halfChannel_H', 'n', 'x_shift'});
    ftopts = fitoptions( 'Method', 'NonlinearLeastSquares','StartPoint', [max(abs(vx_z_phy_mean))*1e6 max(Z)/2 1 0]);
    ftopts.Display = 'Off';
    [fitobject, GOF, fitinfo]  = fit(Z', abs(vx_z_phy_mean)'*1e6, ft, ftopts);
    X = 0:0.1:max(Z);
    plot(X, feval(fitobject,X),'--k','linewidth',1.5)
else
    % fitting: v(r) = v_max*(1-(|r|/halfChannel_H)^(1+1/n)) --> 'v_max', 'halfChannel_H', 'n' [For the cases that don't scan whole channel]
    ft = fittype( 'v_max * (1-(abs(x-x_shift-halfChannel_H)/halfChannel_H)^(1+1/n))', ...
        'independent', 'x', 'dependent', 'y', 'coefficients', {'v_max', 'halfChannel_H', 'n', 'x_shift'});
    ftopts = fitoptions( 'Method', 'NonlinearLeastSquares','StartPoint', [max(abs(vx_z_phy_mean))*1e6 max(Z) 1 0]);
    ftopts.Display = 'Off';
    [fitobject, GOF, fitinfo]  = fit(Z', abs(vx_z_phy_mean)'*1e6, ft, ftopts);
    X = 0:0.1:2*max(Z);
    plot(X, feval(fitobject,X),'--k','linewidth',1.5)
end

fluid_index = fitobject.n;

