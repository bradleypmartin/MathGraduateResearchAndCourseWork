% initializing figure
figure(1), set(gcf, 'Color','white')

% building interface illustration
xInterface = 0.5*ones(5,1); uInterface = linspace(-1.5,1.5,5);

k = 50;
uPlot = ufxt(1+(k-1)*2*N:N+(k-1)*2*N,1);

tSlider = k*tstep*ones(5,1); uSlider = linspace(-1.5,1.5,5);

subplot(1,2,1)
set(gca,'fontSize',12);

plot(xgrid,uPlot,'-k');
title('"u" (true model) over entire domain - sampling location circled')
xlabel('x (L)'); ylabel('local particle velocity (L/t)');
hold on;
plot(xInterface,uInterface,'--k')
axis([-0.1 1.1 -1.1 1.1])
plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',5)
plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',10)
plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',20)
hold off;

subplot(1,2,2)
plot(timeVec,knownu);
hold on;
%plot(tSlider,uSlider,'--b')
plot(tSlider(1,1),uPlot(20,1),'ob','MarkerSize',5)
plot(tSlider(1,1),uPlot(20,1),'ob','MarkerSize',10)
plot(tSlider(1,1),uPlot(20,1),'ob','MarkerSize',20)
hold off;
axis([0 1.5 -1.5 1.5])
set(gca,'fontSize',12)
title('"u" at sampling location - simulated seismogram (data)')
xlabel('time (t)')
ylabel('u (data) at sampling location')