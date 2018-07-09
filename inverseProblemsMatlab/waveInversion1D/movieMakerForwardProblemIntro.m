% initializing figure
figure(1), set(gcf, 'Color','white')

% building interface illustration
xInterface = 0.5*ones(5,1); uInterface = linspace(-1.5,1.5,5);

tSlider = 0*ones(5,1); uSlider = linspace(-1.5,1.5,5);

set(gca,'fontSize',12);
title('"u" (true model) over entire domain')
plot(xgrid,zeros(size(xgrid)),'-k');
xlabel('x (L)'); ylabel('local particle velocity (L/t)');
hold on;
plot(xInterface,uInterface,'--k')
axis([-0.1 1.1 -1.1 1.1])
hold off;

set(gca, 'nextplot','replacechildren', 'Visible','off');

%# create AVI object
nFrames = 300;
vidObj = VideoWriter('ForwardProblemIntro.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);

%# create movie
for k=1:nFrames

    uPlot = ufxt(1+(k-1)*2*N:N+(k-1)*2*N,1);

    plot(xgrid,uPlot,'-k');
    set(gca,'fontSize',12)
    title('"u" (true model) over entire domain')
    xlabel('x (L)'); ylabel('local particle velocity "u" (L/t)');
    hold on;
    plot(xInterface,uInterface,'--k')
    axis([-0.1 1.1 -1.5 1.5])
    hold off;
   
    tSlider = tstep*k*ones(5,1);
    
   writeVideo(vidObj, getframe(figure(1)));
       
end
%close(gcf)

%# save as AVI file, and open it using system video player
close(vidObj);
winopen('ForwardProblemIntro.avi')