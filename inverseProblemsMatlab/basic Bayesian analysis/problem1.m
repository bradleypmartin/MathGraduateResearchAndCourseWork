% Brad Martin - Inverse Methods - HW6, Problem 1

% We'll start by creating a lognormal distribution for mu1.
% +20 from 10 (linear) = *3 from 10 in log space.
mu1centerlog = log(10); mu1errorlog = log(3);   

mu1evalslog = linspace(mu1centerlog-4*mu1errorlog,...
    mu1centerlog+2*mu1errorlog,200);

% normal (Gaussian) distribution in log space
mu1PDF = exp(-(mu1evalslog-mu1centerlog).^2/(2*mu1errorlog^2));
mu1evals = exp(mu1evalslog); % changing log scale back to linear scale

% plotting prior distribution for mu1
figure(1); plot(mu1evals,mu1PDF,'-k');
xlabel('mu1 value'); ylabel('relative probability density')

% We continue by creating a normal distribution for mu2.
mu2center = 50; mu2error = 5;

mu2evals = linspace(mu2center-3*mu2error,...
    mu2center+3*mu2error,200);
mu2PDF = exp(-(mu2evals-mu2center).^2/(2*mu2error^2));

% plotting prior distribution for mu2
figure(2); plot(mu2evals,mu2PDF,'-k');
xlabel('mu2 value'); ylabel('relative probability density')

% Evaluating the bivariate distribution of mu1 and mu2
[MU1 MU2] = meshgrid(mu1evals,mu2evals);
[MU1PDF MU2PDF] = meshgrid(mu1PDF,mu2PDF);

MUbivariatePDF = MU1PDF.*MU2PDF; % assuming independence
% plotting prior bivariate distribution of mu1 and mu2
figure(3); C1 = contour(MU1,MU2,MUbivariatePDF,'-k'); clabel(C1);
xlabel('mu1'); ylabel('mu2')

% Creating the likelihood function for the sum d = mu1 + mu2
dcenter = 80; derror = 10;

devals = MU1+MU2;
likelihoodF = exp(-(devals-dcenter).^2/(2*derror^2));

figure(4); C2 = contour(MU1,MU2,likelihoodF,'-k'); clabel(C2,'manual');
xlabel('mu1'); ylabel('mu2')

% Evaluating the posterior joint distribution 
% for mu1 and mu2 through Bayes

posteriorD = MUbivariatePDF.*likelihoodF;

figure(5); C3 = contour(MU1,MU2,posteriorD,'-k'); clabel(C3);
xlabel('mu1'); ylabel('mu2')

% Approximating posterior mu1 and mu2 distributions by Riemann sums

posteriormu1 = sum(posteriorD,1);
figure(6); plot(mu1evals,posteriormu1/max(posteriormu1),'-k');
xlabel('mu1 value'); ylabel('relative probability density')

posteriormu2 = (sum(posteriorD,2))';
figure(7); plot(mu2evals,posteriormu2/max(posteriormu2),'-k');
xlabel('mu2 value'); ylabel('relative probability density');

% END PROBLEM 1