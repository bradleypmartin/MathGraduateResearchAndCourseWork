knownsigma = 0;

knownuInd = [20];
knownfInd = [20];
ZOswitchu = 1;
ZOswitchf = 1;

knownuAlpha = 1;
knownfAlpha = 1;

knownu = zeros(ceil(endtime/tstep),size(knownuInd,1));
knownf = zeros(ceil(endtime/tstep),size(knownfInd,1));

timeVec = (linspace(0,endtime,ceil(endtime/tstep)))';

for m = 1:ceil(endtime/tstep)
   for n = 1:size(knownuInd,1)
       knownu(m,n) = ufxt(knownuInd(n,1)+(m-1)*2*N,1);
   end
   for n = 1:size(knownfInd,1)
       knownf(m,n) = ufxt(knownfInd(n,1)+N+(m-1)*2*N,1);
   end
end

knownu = knownu+randn(size(knownu))*knownsigma;
knownf = knownf+randn(size(knownu))*knownsigma;

figure(1); plot(timeVec,knownu);
%figure(2); plot(timeVec,knownf);

set(gca,'fontSize',12)
xlabel('t')
ylabel('u (data) at sampling location')