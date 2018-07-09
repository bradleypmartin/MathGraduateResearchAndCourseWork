% Script RESULTSDISPLAY helps display data from RBF-FD modeling of wave
% transport (2D elastic wave equation).
xlb = 0; xub = 1; ylb = 0; yub = 1;

c2Plot = 0;

timeindex = 1;
% Index of the time point at which you want to examine
% data (index in the data storage vector)

N = size(xnodes,1);

plottype = 2;
% 1 = u; 2 = v; 3 = sigma-xx;
% 4 = sigma-xy; 5 = sigma-yy
if plottype < 6
    %
    if c2Plot == 0
        figure(1)
        resultspcolor(xnodes,ynodes,...
            real(datastore(1+N*(plottype-1):plottype*N,timeindex)),xlb,xub,xlb,xub,0.2,500);
        %colormap('gray')
        colormap('jet')
        colorbar
        xlabel('x'); ylabel('y')
    else
        figure(1)
        resultspcolor(xnodes,ynodes,lambdaMuVec,0,1,0,1,0.2,500);
        colormap('jet')
        colorbar
        xlabel('x'); ylabel('y')
    end
    %}
    
else
    %{
    figure(1)
resultspcolor(xnodes,ynodes,...
    log((datastore(1+N*(1-1):1*N,timeindex).^2+...
    datastore(1+N*(2-1):2*N,timeindex).^2).^0.5)/log(10),0,1,0,1,0.2,200);
colormap('jet'); caxis([-7 -2]); colorbar
xlabel('x'); ylabel('y')
    %}
    figure(1)
    resultspcolor(xnodes,ynodes,...
        (datastore(1+N*(1-1):1*N,timeindex).^2+...
        datastore(1+N*(2-1):2*N,timeindex).^2).^0.5,0,1,0,1,0.2,500);
    colormap('jet'); colorbar
    xlabel('x'); ylabel('y')
end
%{
axis([xlb xub xlb xub -0.5 1.0])
caxis([-0.5 1.0])
view([-85 15])
%}