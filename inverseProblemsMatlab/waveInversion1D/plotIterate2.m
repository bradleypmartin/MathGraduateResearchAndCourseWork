endTimePlot = endtime;

uPlot = [];
uPlotLog = [];
xPlot = [];
tPlot = [];

tTemp = tstep;
timeIndex = 1;

while tTemp <= endTimePlot
    uPlot = [uPlot ufxtK(1+(timeIndex-1)*2*N:N+(timeIndex-1)*2*N,1)];
    uPlotLog = [uPlotLog log(abs(ufxtK(1+(timeIndex-1)*2*N:N+(timeIndex-1)*2*N,1)))];
    xPlot = [xPlot xgrid];
    tPlot = [tPlot ones(size(xgrid))*tTemp];
    timeIndex = timeIndex+1;
    tTemp = tTemp+tstep;
end

figure(1)
pcolor(xPlot,tPlot,uPlot); shading('interp')
%axis equal
set(gca,'fontSize',12)
colorbar
xlabel('x (L)')
ylabel('time (t)')
zlabel('u(x,t)')
title('S-T colormap of "u" (6th guess; inv. problem; 1st ord. reg.*)')

figure(2)
pcolor(xPlot,tPlot,uPlotLog); shading('interp')
%axis equal
set(gca,'fontSize',12)
colorbar
xlabel('x (L)')
ylabel('time (t)')
zlabel('u(x,t)')
title('S-T colormap of log(u) (sixth guess; inv. problem; 1st ord. reg.)')
