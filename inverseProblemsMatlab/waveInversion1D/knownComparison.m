knownuComp = zeros(size(knownu));
knownfComp = zeros(size(knownu));

knownIndices = [20];

for m = 1:size(knownu,1)
   for n = 1:size(knownIndices,1)
   knownuComp(m,n) = ufxtK(knownIndices(n,1)+(m-1)*2*N,1);
   knownfComp(m,n) = ufxtK(knownIndices(n,1)+N+(m-1)*2*N,1);
   end
end

figure(3); plot(timeVec,knownu);
set(gca,'fontSize',12)
xlabel('time point')
ylabel('u readout at sample location (x = 0.10); known')
axis([0 1.5 -1.5 1.5])
figure(4); plot(timeVec,knownuComp); axis([0 1.5 -1.5 1.5])
set(gca,'fontSize',12)
xlabel('time point')
ylabel('u readout at sample location (x = 0.10); INV solution')

%figure(5); plot(knownf);
%figure(6); plot(knownfComp);
figure(7); plot(timeVec,knownu-knownuComp); axis([0 1.5 -1.5 1.5])
set(gca,'fontSize',12)
xlabel('time point')
ylabel('known u - inverted u, sample location (x = 0.10)')