N = 200; tstep = 0.005; endtime = 1.5; dampConst = 300; naiveFlag = 0;
spaceSigma = 0.015; timeSigma = 0.015; exLoc = 0.2; exTime = 0.1;

% orig (monochrome) 1st guess
intStack = [0.50001 1 1];

% one int test
%intStack = [0.50001 2 1.5];

% close ints test
%intStack = [0.4001 2 1; 0.4301 1 1];

% 3 interfaces base
%intStack = [0.31001 2 1; 0.56001 3 1; 0.83001 2 1];

% Naive 200 3 interfaces seq
%intStack = [0.31001 3.2^0.5 1];
%intStack = [0.31001 3.15^0.5 1; 0.53001 7.25^0.5 1];

% AC 200 3 interfaces seq
%intStack = [0.31001 2 1];
%intStack = [0.31001 2.9^0.5 1; 0.52001 7.5^0.5 1];

% Naive 400 3 interfaces seq
%intStack = [0.31001 3.6^0.5 1];
%intStack = [0.31001 3.55^0.5 1; 0.545001 7.75^0.5 1];

% AC 400 3 interfaces seq
%intStack = [0.31001 2 1];
%intStack = [0.31001 4^0.5 1; 0.56001 8.92^0.5 1];

[xgrid ufxt bigOperator bigSize] = ...
    FD4ACST7(N,tstep,endtime,dampConst,naiveFlag,...
    spaceSigma,timeSigma,exLoc,exTime,intStack);

% 400 x 400 - t = 24.4 sec
% 200 x 200 - t = 3.2 sec
% 100 x 100 - t = 0.43 sec
% 50 x 50 - t = 0.077 sec