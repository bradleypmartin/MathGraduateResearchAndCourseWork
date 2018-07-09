endTimePlot = endtime;

uPlot = [];
xPlot = [];
tPlot = [];

tTemp = tstep;
timeIndex = 1;

while tTemp <= endTimePlot
    uPlot = [uPlot ufxt(1+(timeIndex-1)*2*N:N+(timeIndex-1)*2*N,1)];
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
title('Space-time colormap of "u" (first guess; inv. problem)')
%{
% interface plotting
xIntPlot = 0.5*ones(5,1);
tIntPlot = linspace(0,1.5,5);
hold on;
plot(xIntPlot,tIntPlot,'--k')
hold off;
%}