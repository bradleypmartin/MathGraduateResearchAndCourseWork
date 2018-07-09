figure(4); plot(xgrid,ufxtK(2*bigSize+1:2*bigSize+N));
set(gca,'fontSize',12)
xlabel('x')
ylabel('1 / rho')
%axis([0 1 0 1.1])
%{
figure(4); plot(xgrid,log(abs(ufxtK(2*bigSize+1:2*bigSize+N))));
set(gca,'fontSize',12)
xlabel('x')
ylabel('log ( 1 / rho )')
%}
figure(6); plot(xgrid,ufxtK(2*bigSize+N+1:2*bigSize+2*N));
set(gca,'fontSize',12)
xlabel('x')
ylabel('rho * c^2')
%axis([0 1 0 7])
%{
figure(6); plot(xgrid,log(abs(ufxtK(2*bigSize+N+1:2*bigSize+2*N))));
set(gca,'fontSize',12)
xlabel('x')
ylabel('log ( rho * c^2 )')
%}