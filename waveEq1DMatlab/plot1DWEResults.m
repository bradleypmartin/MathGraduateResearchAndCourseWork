% Here, you can plot snapshots of wave problem results obtained by running
% the runDriver1DWE.m script.

viewTime = 0.75;
% time for viewing snapshot (in the unitless time of the
% simulation; i.e. a wave front in a region with c = 1 will
% have traversed 1 spatial unit in 1 unit of time).
              
% computing index of data to view at viewTime
timeIndex = round(viewTime/tstep)+1;

% plotting pressure/stress data
figure(1);
plot(xgrid,fxt(:,timeIndex),'-k','linewidth',1.5);
hold on;

% plotting interfaces that mark boundary of the heterogeneous region
plot(zeros(1,10),linspace(-1,1,10),'--k');
plot(interfaceWidth*ones(1,10),linspace(-1,1,10),'--k');
hold off;

% enlarging fonts and labeling axes
set(gca,'fontSize',14);
title('Plot of $f$ (stress) at requested time; 1D wave eq.','interpreter','latex');
xlabel('$x$ (unitless)','interpreter','latex');
ylabel('$f$ (stress; unitless)','interpreter','latex');

% scaling axes appropriately for current snapshot
axis([-1.2 1.2 min(-1.2,min(min(fxt(:,timeIndex)))-0.1) ...
    max(1.2,max(max(fxt(:,timeIndex)))+0.1)])