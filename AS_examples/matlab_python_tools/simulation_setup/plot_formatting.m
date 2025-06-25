%matlab plot formatting

linewidth = 3;
fontsize = 17;
legendfontsize = 15;
markersize = 10;

set(groot,'defaultlinelinewidth',linewidth)
set(groot,'defaultaxesfontsize',fontsize)
set(groot,'defaultaxestickdir','in')
set(groot,'defaultaxesticklength',[.02 .02]);
set(groot,'defaultaxesXMinorTick','on');
set(groot,'defaultaxesYMinorTick','on');
set(groot,'defaultaxesTickLabelInterpreter', 'latex');
set(groot,'defaultTextInterpreter', 'latex');
set(groot,'DefaultAxesXGrid','on');
set(groot,'DefaultAxesYGrid','on');
set(groot,'DefaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual')
set(groot,'DefaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'DefaultAxesBox','on');
set(groot,'DefaultLegendFontSize',legendfontsize);
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'defaultlinemarkersize',markersize)
set(groot,'defaultfigureunits','pixels');
set(groot,'defaultfigureposition', [100 100 800 700]);
