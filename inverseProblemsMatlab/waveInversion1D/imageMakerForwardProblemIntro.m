% initializing figure
figure(1), set(gcf, 'Color','white')

% building interface illustration
xInterface = 0.5*ones(5,1); uInterface = linspace(-1.5,1.5,5);

tSlider = 0*ones(5,1); uSlider = linspace(-1.5,1.5,5);

k = 100;
uPlot = ufxt(1+(k-1)*2*N:N+(k-1)*2*N,1);

set(gca,'fontSize',12);

plot(xgrid,uPlot,'-k');
title('"u" (true model) over entire domain')
xlabel('x (L)'); ylabel('local particle velocity (L/t)');
hold on;
plot(xInterface,uInterface,'--k')
axis([-0.1 1.1 -1.1 1.1])
hold off;