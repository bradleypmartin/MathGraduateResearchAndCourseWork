% Brad Martin - Inverse Methods, HW6, problem 3 (Aster 11.4)

% PART B
bevals = linspace(0,1,200); % I'll plot over b = [0,1].

% in part B, P(d1 = 0|b) is just the probability that we get one roll
% of heads, given a certain bias. Given the convention in the problem
% statement, this is proportional to b.
dprobsB = bevals;

% since we've got an uninformative prior, P(d1 = 0|b) is directly
% proportional to P(b|d1 = 0), so these are qualitatively equivalent.
figure(1); plot(bevals,dprobsB,'-k')
xlabel('bias (b)'); ylabel('relative probability density (P(b|d1=0))')

% PART C
% Part C is a little trickier - P(d1 = 0, ... d5 = 1|b) is proportional to
% the probability of getting a head then 4 tails with a certain bias, or
% (1-b)^4*b^1.  As in part B, with an uninformative prior, our a
% posteriori distribution just looks like our likelihood function.
dprobsC = (1-bevals).^4.*bevals.^1;
figure(2); plot(bevals,dprobsC,'-k')
xlabel('bias (b)'); ylabel('relative probability density P(b|1h,4t)')

% PART E
% Here we just follow Bayes' theorem as done in Problem 1: we get
% P(b|(data from part C)) ~ P(b)*P((data from part C)|b).
bPDF = exp(-10*(0.5-bevals).^2);
bprobsE = dprobsC.*bPDF;
figure(3); plot(bevals,bprobsE,'-k')
xlabel('bias (b)'); ylabel('relative probability density (post., part E)')

% END PROBLEM 3