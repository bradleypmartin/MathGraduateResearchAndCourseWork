% initializing figure
figure(1), set(gcf, 'Color','white')

% building interface illustration
xInterface = 0.5*ones(5,1); uInterface = linspace(-1.5,1.5,5);

tSlider = 0*ones(5,1); uSlider = linspace(-1.5,1.5,5);

subplot(1,2,1)
set(gca,'fontSize',12);
title('"u" (true model) over entire domain - sampling location circled')
plot(xgrid,zeros(size(xgrid)),'-k');
xlabel('x (L)'); ylabel('local particle velocity (L/t)');
hold on;
plot(xInterface,uInterface,'--k')
axis([-0.1 1.1 -1.1 1.1])
plot(xgrid(20,1),0,'ob','MarkerSize',5)
plot(xgrid(20,1),0,'ob','MarkerSize',10)
plot(xgrid(20,1),0,'ob','MarkerSize',20)
hold off;

subplot(1,2,2)
plot(timeVec,knownu);
hold on;
%plot(tSlider,uSlider,'--b')
plot(tSlider(1,1),0,'ob','MarkerSize',5)
plot(tSlider(1,1),0,'ob','MarkerSize',10)
plot(tSlider(1,1),0,'ob','MarkerSize',20)
hold off;
axis([0 1.5 -1.5 1.5])
set(gca,'fontSize',12)
title('"u" at sampling location - simulated seismogram (data)')
xlabel('time (t)')
ylabel('u (data) at sampling location')

set(gca, 'nextplot','replacechildren', 'Visible','off');

%# create AVI object
nFrames = 300;
vidObj = VideoWriter('InvProblemIntro.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);

%# create movie
for k=1:nFrames

    uPlot = ufxt(1+(k-1)*2*N:N+(k-1)*2*N,1);

    subplot(1,2,1)
    plot(xgrid,uPlot,'-k');
    set(gca,'fontSize',12)
    title('"u" (true model) over entire domain - sampling location circled')
    xlabel('x (L)'); ylabel('local particle velocity "u" (L/t)');
    hold on;
    plot(xInterface,uInterface,'--k')
    plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',5)
    plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',10)
    plot(xgrid(20,1),uPlot(20,1),'ob','MarkerSize',20)
    axis([-0.1 1.1 -1.5 1.5])
    hold off;

    subplot(1,2,2)
    
    tSlider = tstep*k*ones(5,1);
    
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
    ylabel('local particle velocity "u" (data) at sampling location')
    
   writeVideo(vidObj, getframe(figure(1)));
       
end
%close(gcf)

%# save as AVI file, and open it using system video player
close(vidObj);
winopen('InvProblemIntro.avi')