% Script RESULTSDISPLAY helps display data from RBF-FD modeling of wave
% transport (2D elastic wave equation).

c2Plot = 0;

timeindex = 11;
% Index of the time point at which you want to examine
% data (index in the data storage vector)

N = size(xnodes,1);

plottype = 1;        
% 1 = u; 2 = v; 3 = sigma-xx; 
% 4 = sigma-xy; 5 = sigma-yy
if plottype < 4
%
if c2Plot == 0
    figure(1)
    resultssurf(xnodes,ynodes,...
        real(datastore(1+N*(plottype-1):plottype*N,timeindex)),0,1,0,1,0.2,500);
    colormap('jet')
    colorbar
    xlabel('x'); ylabel('y')
else
    figure(1)
    resultspcolor(xnodes,ynodes,tripleVec,0,1,0,1,0.2,500);
    colormap('jet')
    colorbar
    xlabel('x'); ylabel('y')
end
%}

else
    %{
    figure(1)
resultspcolor(xnodes,ynodes,...
    log((datastoreRBFnaive(1+N*(1-1):1*N,timeindex).^2+...
    datastoreRBFnaive(1+N*(2-1):2*N,timeindex).^2).^0.5)/log(10),0,1,0,1,0.2,200);
colormap('jet'); caxis([-7 -3]); colorbar
xlabel('x'); ylabel('y')
%}
figure(2)
resultspcolor(xnodes,ynodes,...
    (datastoreRBFnaive(1+N*(1-1):1*N,timeindex).^2+...
    datastoreRBFnaive(1+N*(2-1):2*N,timeindex).^2).^0.5,0,1,0,1,0.2,200);
colormap('jet'); colorbar
xlabel('x'); ylabel('y')
end