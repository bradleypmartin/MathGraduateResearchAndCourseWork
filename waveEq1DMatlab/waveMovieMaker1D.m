% Here, we're going to make a high-quality .avi video file of our wave
% simulation results.  This is adapted code provided by Mike Nicholas
% (Colorado School of Mines) and used in the Fall 2017 Intro to Scientific
% Computing course (MATH 307).  It was adapted by Brad Martin
% (bmartin@mines.edu) on 180119 for use with wave and heat simulation 
% problems.

% Make sure to run the runDriver1DWE.m script first so this video maker has
% data to work with.

M = VideoWriter('wave1D.avi');  % opening a new .avi video file
M.FrameRate = 60;               % declaring a 60 FPS framerate
M.Quality   = 100;              % choosing high-quality recording (100 MAX)
open(M);

figure(1);  % We'll use snapshots in figure(1) to write the video
            % frame by frame.
for i = 1:(60*endtime*slowMult+1)
    cla;                                   % clearing the frame
    viewTime = (i-1)/(60*slowMult);        % time within sim. scale
    timeIndex = round(viewTime/tstep)+1;   % index of data to write
    figure(1);
    % plotting pressure/stress data
    plot(xgrid,fxt(:,timeIndex),'-k','linewidth',1.5);
    hold on;
    % plotting interfaces
    plot(zeros(1,10),linspace(-1,1,10),'--k');
    plot(interfaceWidth*ones(1,10),linspace(-1,1,10),'--k');
    hold off;
    set(gca,'fontSize',14) % enlarging fonts; labeling plot/axes
    title('Plot of $f$ (stress) over time; 1D wave eq.','interpreter','latex');
    xlabel('$x$ (unitless)','interpreter','latex');
    ylabel('$f$ (stress; unitless)','interpreter','latex');
    % setting axis limits
    axis([-1.2 1.2 min(-1.2,min(min(fxt))-0.1) max(1.2,max(max(fxt))+0.1)])
    
    % Capturing current frame and store into movie matrix
    frame = getframe(gcf);  % get current figure (includes axes info)
    writeVideo(M,frame);    % store it in the video M
end

close(M);  % needed so you can open your avi with other applications