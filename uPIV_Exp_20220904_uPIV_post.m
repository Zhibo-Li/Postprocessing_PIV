%% tmp code for checking the influence of the fiber on the flow.

%%% convert *.vc7 to *.mat and convert to physical unit (um/s) and plotting
% direction: x -- spanwise
%            y -- streamwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

magni = 0.42; % the magnification of the objective (um/pixel)

pathname = uigetdir(['Z:\Experimental Data (RAW)\FSI - Rigid Fiber &  Individual Obstacle' ...
    '\20220904_Chip20220624_DualTri_SU8_Tracer\Lateral=5nL_Main=1nL_Lens=10_' ...
    'Tracer=1.1um_Atte=100_Z=mid_Phi100um_03\AddGeometricMask_02\ImgPreproc\']);
listing = dir(fullfile(pathname, '*.vc7')); suffiex = 'vc7'; % Choose this if you need to convert.
% listing = dir(fullfile(pathname, '*.mat')); suffiex = 'mat'; % Choose this if there are *.mat files.

for ii = 1:length(listing)

    theONE = listing(ii).name;
%     Z(ii) = str2double(extractBetween(theONE,'Z=','um'));
    deltaT = 10000;
    if strcmp(suffiex,'vc7')

        v = loadvec(fullfile(pathname,theONE));  % to calculate the physical field

        ave_field_phy = v; % in physical units.
        ave_field_phy.vx = v.vx * magni / deltaT * 1e6; % velocity in um/s.
        ave_field_phy.vy = v.vy * magni / deltaT * 1e6; % velocity in um/s.
        ave_field_phy.x = v.x * magni; % velocity in um.
        ave_field_phy.y = v.y * magni; % velocity in um.
        ave_field_phy.unitvx = '\mu{m}/s'; ave_field_phy.unitvy = '\mu{m}/s';
        ave_field_phy.unitx = '\mu{m}'; ave_field_phy.unity = '\mu{m}';

%         save([pathname,filesep,theONE(1:end-4),'.mat'], "ave_field_phy");

    else

        load(fullfile(pathname,theONE))

    end

    %%%%% the U and V along the z-direction (channel height) %%%%%
    if ii == 1
        imshow(imbinarize(ave_field_phy.vy));   
        selrect = getrect; close; % select the area to average.
    end
    select_U = ave_field_phy.vx(selrect(2):selrect(2)+selrect(4), selrect(1):selrect(1)+selrect(3));
    select_V = ave_field_phy.vy(selrect(2):selrect(2)+selrect(4), selrect(1):selrect(1)+selrect(3));
    UUU_alongZ(ii) = mean(mean(select_U));
    VVV_alongZ(ii) = mean(mean(select_V)); 

end

t = 1000/15 * [0:length(UUU_alongZ)-1];
%%% plot V along the z-direction
figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
plot(t, sqrt(VVV_alongZ.^2 + UUU_alongZ.^2), 'LineStyle','none', 'Marker','.','MarkerSize', 30); hold on
xlabel('Time\ (ms)','FontSize',22,'Interpreter', 'latex');
ylabel('$|U|\ \rm{({\mu}m/s)}$','FontSize',22,'Interpreter', 'latex');
set(gca,'FontSize',20)

set_plot(gcf, gca)

% f=gcf;
% exportgraphics(f,'F:\Processing & Results\Characterization - Particles & Fibers\20220704-SU8_Fibers_Size\lengths_distribution.eps')

function set_plot( current_fig, current_axes )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% To avoid figures with random size (matlab bug)
pause(0.1);    
drawnow; 

lineStyle = {'-','--',':','-.'};
assignin('caller','lineStyle',lineStyle);

set(current_fig, ...
    'Position',[100 100 800 600], ...
    'Color','white', ...
    'DefaultTextInterpreter', 'latex')

set(current_axes, ...
    'Box', 'On', ...
    'XGrid', 'On', ...
    'YGrid', 'On', ...
    'GridAlpha', 0.5, ...
    'FontSize', 24, ...
    'NextPlot','replacechildren', ...
    'TickLabelInterpreter','latex')

current_axes.XLabel.Color = [0 0 0];
current_axes.YLabel.Color = [0 0 0];
end



