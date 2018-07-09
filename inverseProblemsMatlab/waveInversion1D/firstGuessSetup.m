% Zero initial guess
%ufxtK = [zeros(size(ufxt)); ones(2*N,1)];
%ufxtK = [zeros(size(ufxt)); ones(N,1); 4*ones(N,1)];

% "monochrome" Initial guess
ufxtK = [ufxt; ones(2*N,1)];

% Softball guess
%ufxtK = [ufxt; ones(N,1); 4*ones(N,1)];

% Naive - 3 ints sequential 200
% ufxtK = [ufxt; ones(200,1); ones(63,1); 3.2*ones(137,1)];
%ufxtK = [ufxt; ones(200,1); ones(63,1); 3.15*ones(44,1); 7.25*ones(93,1)];

% AC - 3 ints sequential 200
%ufxtK = [ufxt; ones(200,1); ones(63,1); 4*ones(137,1)];
%ufxtK = [ufxt; ones(200,1); ones(63,1); 2.9*ones(42,1); 7.5*ones(95,1)];

% Naive - 3 ints sequential 400
%ufxtK = [ufxt; ones(400,1); ones(125,1); 3.6*ones(275,1)];
%ufxtK = [ufxt; ones(400,1); ones(125,1); 3.55*ones(94,1); 7.75*ones(181,1)];

% AC - 3 ints sequential 400
%ufxtK = [ufxt; ones(400,1); ones(125,1); 4*ones(275,1)];
%ufxtK = [ufxt; ones(400,1); ones(125,1); 4*ones(100,1); 8.92*ones(175,1)];