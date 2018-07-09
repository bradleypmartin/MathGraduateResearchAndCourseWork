%%% Brad Martin, Inv. Meth. HW6, prob. 2 (revisiting HW3, prob 1b)

load('G.txt'); load('m0.txt'); load('d.txt')  % loading data, m0
load('TikSolution.mat')  % loading 1st order Tikhonov solution

rows = size(G,1); cols = size(G,2);  % sizing up G

%%% PART B - Kalman filter (Bayes solution); uncorrelated error
yVec = (0.0250:0.05:0.9750)';
CPartB = diag(ones(cols,1)*2000^2);

m0 = m0+0.3*sin(1*pi*yVec);
GAug = [G; CPartB^-0.5]; dAug = [d; CPartB^-0.5*m0];
mPartB = GAug\dAug;

%%% PART D - plotting the solutions together
%

figure(3); plot(yVec,m0,'--k',yVec,mPartA,'sq-k',...
    yVec,mPartB,'-*k')
legend('m0','m (Part A)','m (Part B)');
xlabel('y value'); ylabel('model value m(y)')
axis([0 1 0.1 1])
%}
%%% END PROBLEM 1

(G*mPartB-d)'*(G*mPartB-d)
(G*m0-d)'*(G*m0-d)