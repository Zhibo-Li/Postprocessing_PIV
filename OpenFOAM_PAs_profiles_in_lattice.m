%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

for Array_angle = [0 45]
    
    % read simulation data
    Data = readcell(['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_MidPlane_Inlattice_components.csv']);
    Data_Titles = Data(1, :);
    Data_withoutTitles = cell2mat(Data(2:end, :));
    
    % find the velocity components
    velocity_column_ind = find(contains(Data_Titles, 'U:'));
    velocity = Data_withoutTitles(:, velocity_column_ind);
    
    % find the velocity gradient components
    velocity_gradient_column_ind = find(contains(Data_Titles, 'grad(U):'));
    velocity_gradient = Data_withoutTitles(:, velocity_gradient_column_ind);
    % reshape the velocity gradient matrix
    flowType_Zhibo = zeros(size(velocity_gradient, 1), 1);
    for ii = 1: size(velocity_gradient, 1)

        data = reshape(velocity_gradient(ii, :), [3, 3]);
        data_sym = (data+data.')/2;
        data_anti = (data-data.')/2;
        
        norm_E = norm(data_sym,'fro');
        norm_O = norm(data_anti,'fro');
        flowType_Zhibo(ii) = (norm_E-norm_O) / (norm_E+norm_O); % To check the flow-type parameter.

    end
    
    % find the flowType
    flowType_column_ind = find(contains(Data_Titles, 'flowType'));
    flowType = Data_withoutTitles(:, flowType_column_ind);
    
    % find the flowStrength
    flowStrength_column_ind = find(contains(Data_Titles, 'flowStrength'));
    flowStrength = Data_withoutTitles(:, flowStrength_column_ind);
    
    % find the coordinates
    position_column_ind = find(contains(Data_Titles, 'Points:'));
    positiion = Data_withoutTitles(:, position_column_ind);
    % shift the position by minus the start point and only choose x and y coordinates
    positiion = positiion - positiion(1, :); positiion = positiion(:, 1:2);
    normalized_position = sqrt(sum(positiion.^2, 2)) / max(sqrt(sum(positiion.^2, 2)));
    
    % plot the velocity gradient profiles in the lattice
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
    
    plot(normalized_position, velocity_gradient(:,1), 'k*'); hold on;
    plot(normalized_position, velocity_gradient(:,2), 'ro'); hold on;
    plot(normalized_position, velocity_gradient(:,4), 'bd', 'MarkerSize', 1); hold on;
    plot(normalized_position, velocity_gradient(:,5), 'c^'); hold on;
    
    xlim([0 1]);
    xlabel('Normalized Positions','FontSize', 24,'FontName', 'Times New Roman');
    ylabel('Values','FontSize', 24,'FontName', 'Times New Roman');
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24,'TickLabelInterpreter','latex')
    legend('dU/dx', 'dU/dy', 'dV/dx', 'dV/dy', 'Location', 'best', 'FontSize', 24,'FontName', 'Times New Roman');
    
    hhh = gcf;
    set(hhh,'Units','Inches');
    pos = get(hhh,'Position');
    set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hhh, '-dpdf',['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\figures\VelocityGradient_in_lattice.pdf']);

    % plot the flozType and flowStrength profiles in the lattice
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);

    plot(normalized_position, flowType, 'k*'); hold on;
    plot(normalized_position, flowStrength, 'ro'); hold on;

    xlim([0 1]);
    xlabel('Normalized Positions','FontSize', 24,'FontName', 'Times New Roman');
    ylabel('Values','FontSize', 24,'FontName', 'Times New Roman');
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24,'TickLabelInterpreter','latex')
    legend('flowType', 'flowStrength', 'Location', 'best', 'FontSize', 24,'FontName', 'Times New Roman');

    hhh = gcf;
    set(hhh,'Units','Inches');
    pos = get(hhh,'Position');
    set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hhh, '-dpdf',['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\figures\flowType_and_flowStrength_in_lattice.pdf']);
    
    close all
    
end