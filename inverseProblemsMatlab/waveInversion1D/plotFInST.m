endTimePlot = endtime;
plotType = 0; % (0 = u; 1 = f)

FPlot = [];
xPlot = [];
tPlot = [];

tTemp = tstep;
timeIndex = 1;

while tTemp <= endTimePlot
    FPlot = [FPlot F(1+plotType*N+(timeIndex-1)*2*N:N+plotType*N+(timeIndex-1)*2*N,1)];
    xPlot = [xPlot xgrid];
    tPlot = [tPlot ones(size(xgrid))*tTemp];
    timeIndex = timeIndex+1;
    tTemp = tTemp+tstep;
end

figure(1)
surf(xPlot,tPlot,FPlot); shading('interp')
%axis([-0.2 1.2 -1.2 1.2])
set(gca,'fontSize',12)
colorbar
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
