%%% Brad Martin, Inv. Meth. HW6, prob. 2 (revisiting HW3, prob 1b)

load('G.txt'); load('m0.txt'); load('d.txt')  % loading data, m0
load('TikSolution.mat')  % loading 1st order Tikhonov solution
                         % (saved from HW3)

rows = size(G,1); cols = size(G,2);  % sizing up G

% Kalman filter (Bayes solution); uncorrelated error -
% reproduced from HW3, part B
yVec = (0.0250:0.05:0.9750)'; Cm0 = diag(ones(cols,1)*2000^2);
GAug = [G; CPartB^-0.5]; dAug = [d; Cm0^-0.5*m0];
mPartB = GAug\dAug;

figure(1); plot(yVec,m0,'--k',yVec,mPartA,'sq-k',...
    yVec,mPartB,'-*k')
legend('m0, HW3, P. 1b','m (Tik. 1st O.)','m (Part B)');
xlabel('y value'); ylabel('model value m(y)'); axis([0 1 0.1 1])

% Bayes solution - m01 determined by smooth perturbation to m0
m01 = m0+0.3*sin(1*pi*yVec); dAug = [d; Cm0^-0.5*m01]; m01post = GAug\dAug;

figure(2); plot(yVec,m01,'--k',yVec,mPartA,'sq-k',...
    yVec,m01post,'-*k')
legend('m01 (error = 2000)','m (Tik. 1st O.)','Bayes solution from m01');
xlabel('y value'); ylabel('model value m(y)'); axis([0 1 0.1 1])

% Bayes solution - m02 determined by oscillatory perturbation to m0
m02 = m0+0.3*sin(4*pi*yVec); dAug = [d; Cm0^-0.5*m02]; m02post = GAug\dAug;

figure(3); plot(yVec,m02,'--k',yVec,mPartA,'sq-k',...
    yVec,m02post,'-*k')
legend('m02 (error = 2000)','m (Tik. 1st O.)','Bayes solution from m02');
xlabel('y value'); ylabel('model value m(y)'); axis([0 1 0.1 1])

% Bayes solution - m01 smooth perturbation; std. error = 500
CmSurer = diag(ones(cols,1)*500^2); GAug = [G; CmSurer^-0.5];
dAug = [d; CmSurer^-0.5*m01]; m01post = GAug\dAug;

figure(4); plot(yVec,m01,'--k',yVec,mPartA,'sq-k',...
    yVec,m01post,'-*k')
legend('m01 (error = 500)','m (Tik. 1st O.)','Bayes solution, m01');
xlabel('y value'); ylabel('model value m(y)'); axis([0 1 0.1 1])

% Bayes solution - m02 oscillatory perturbation; std. error = 500
dAug = [d; CmSurer^-0.5*m02]; m02post = GAug\dAug;

figure(5); plot(yVec,m02,'--k',yVec,mPartA,'sq-k',...
    yVec,m02post,'-*k')
legend('m02 (error = 500)','m (Tik. 1st O.)','Bayes solution, m02');
xlabel('y value'); ylabel('model value m(y)'); axis([0 1 0.1 1])

% END PROBLEM 2