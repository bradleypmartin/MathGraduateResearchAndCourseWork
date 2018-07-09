%%% Brad Martin, Inv. Meth. HW3, prob. 1

load('G.txt'); load('m0.txt'); load('d.txt')  % loading data, m0

%%% PART A - 1st-order Tikhonov regularization

rows = size(G,1); cols = size(G,2);  % sizing up G

L = zeros(cols-1,cols);       % initializing 1st-order L block
zeroBlock = zeros(cols-1,1);  % augmented RHS

for counter1 = 1:cols-1         % filling out the L block
   L(counter1,counter1) = -1;   % (proportional to FD approx of 
   L(counter1,counter1+1) = 1;  % 1st deriv.)
end

% below: setting bounds for and testing a range of alpha values

alphaStart = 1E-9; alphaEnd = 1E0;
numTrials = 30;
resNormStack = zeros(1,numTrials); LNormStack = zeros(1,numTrials);

logAlphaStart = log(alphaStart); logAlphaEnd = log(alphaEnd);
logAlphaVec = linspace(logAlphaStart,logAlphaEnd,numTrials);
alphaVec = exp(logAlphaVec);

index = 16;   % circled point on L-curve corresponding to 16th
              % alpha value (4.5 E-05)

for counter1 = 1:numTrials
    GAug = [G; alphaVec(1,counter1)*L]; dAug = [d; zeroBlock];
    mTemp = GAug\dAug;
    resNormStack(1,counter1) = norm(d-G*mTemp);
    LNormStack(1,counter1) = norm(L*mTemp);
end

% plotting the L-curve

figure(1);
plot(log(resNormStack),log(LNormStack),'.k',...
    log(resNormStack(1,index)),log(LNormStack(1,index)),'ok',...
    'MarkerSize',10);
xlabel('log ( || Gm-d || )'); ylabel('log ( || Lm || )')

% producing a model from the corner of the L-curve (alpha = 4.5 E-05)

GAug = [G; alphaVec(1,16)*L]; dAug = [d; zeroBlock];
mPartA = GAug\dAug;
    
    
%%% PART B - Kalman filter; uncorrelated error

CPartB = diag(ones(cols,1)*2000^2);

GAug = [G; CPartB^-0.5]; dAug = [d; CPartB^-0.5*m0];
mPartB = GAug\dAug;


%%% PART C - Kalman filter; block-correlated error

CBlock = 0.1*eye(4)+0.9*ones(4);
CPartC = eye(20);

% Below: creating the block-correlation matrix

for blockCounter = 1:5
   CPartC(1+(blockCounter-1)*4:blockCounter*4,...
       1+(blockCounter-1)*4:blockCounter*4) = 2000^2*CBlock;
end

GAug = [G; CPartC^-0.5]; dAug = [d; CPartC^-0.5*m0];
mPartC = GAug\dAug;


%%% PART D - plotting the solutions together

yVec = (0.0250:0.05:0.9750)';
figure(2); plot(yVec,m0,'--k',yVec,mPartA,'sq-k',...
    yVec,mPartB,'-*k',yVec,mPartC,'o-k')
legend('m0','m (Part A)','m (Part B)','m (Part C)');
xlabel('x value'); ylabel('model value m(x)')

%%% END PROBLEM 1