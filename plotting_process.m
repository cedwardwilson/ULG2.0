%% Plotting generation method
f1 = openfig('/media/chris/Project1/uDALES_veg/ULG2.0/process/5g_step1.fig')
ax1 = gca; % get handle to axes of figure
pause(0.1)
f2 = openfig('/media/chris/Project1/uDALES_veg/ULG2.0/process/5g_step2.fig')
ax2 = gca; % get handle to axes of figure
pause(0.1)
f3 = openfig('/media/chris/Project1/uDALES_veg/ULG2.0/process/5g_step3.fig')
ax3 = gca; % get handle to axes of figure
pause(0.1)
f4 = openfig('/media/chris/Project1/uDALES_veg/ULG2.0/process/5g_step4.fig')
ax4 = gca; % get handle to axes of figure
pause(0.1)

f = figure;
tcl=tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
ax1.Parent=tcl;
ax1.Layout.Tile=1;
ax2.Parent=tcl;
ax2.Layout.Tile=2;
ax3.Parent=tcl;
ax3.Layout.Tile=3;
ax4.Parent=tcl;
ax4.Layout.Tile=4;

% fig1 = get(ax1,'children'); %get handle to all the children in the figure
% fig2 = get(ax2,'children');

