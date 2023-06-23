%% Umag
figure('color', 'w'); set(gcf, 'Position', [100 100 650 200]);

cmocean('speed');
caxis([0 5])
c = colorbar;
c.Label.String = '$|U|/|U_0|$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 14;
c.Location = 'north';
axis off

f=gcf;
savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
    '\PAs_Umag_Colorbar.fig'])
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\5-flexible_fiber_array\PAs_Umag_Colorbar.eps']);



%% flowType
figure('color', 'w'); set(gcf, 'Position', [100 100 650 200]);

cmocean('balance');
caxis([-1 1])
c = colorbar;
c.Label.String = 'Flow-type parameter: $\xi$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 14;
c.Location = 'north';
axis off

f=gcf;
savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
    '\PAs_flowType_Colorbar.fig'])
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\5-flexible_fiber_array\PAs_flowType_Colorbar.eps']);



%% flowStrength
figure('color', 'w'); set(gcf, 'Position', [100 100 650 200]);

cmocean('dense');
caxis([0 3])
c = colorbar;
c.Label.String = 'Flow-strength parameter: $\sigma$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 14;
c.Location = 'north';
axis off

f=gcf;
savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
    '\PAs_flowStrength_Colorbar.fig'])
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\5-flexible_fiber_array\PAs_flowStrength_Colorbar.eps']);

   