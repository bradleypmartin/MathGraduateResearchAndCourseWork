%%% Brad Martin - Inverse Methods, HW 4, problem 1

load('blur.mat')  % Loading G and d from Aster's materials (prob. 6.5)

%%% PART A - Viewing sparsity structure of G...

% ...in a full 10,000 x 10,000 frame...
figure(1); spy(G); xlabel('column index'); ylabel('row index')

% ...zoomed in to 500 x 500 around the main diagonal...
figure(2); spy(G); axis([0 500 0 500]);
xlabel('column index'); ylabel('row index')

% ...and finally, zoomed to 30 x 30 around the main diagonal.
figure(3); spy(G); axis([0 30 0 30])
xlabel('column index'); ylabel('row index')


%%% PART B - memory requirements assessment

gNonZeros = nnz(G);   % counting non-zero values in G matrix

memReqSpG = gNonZeros*32/8  % sparse G storage requirement
memReqFullG = 10000^2*32/8  % full G storage requirement
memReqSVDG = memReqFullG*2+10000*32/8  % SVD storage requirement


%%% PART C - plotting the raw, blurred image in vector d.

img = reshape(d,100,100);
figure(4); imagesc(img); colormap('gray')
xlabel('x location in image (out of 100)');
ylabel('y location in image (out of 100)');


%%% PART D/E - using CGLS algorithm to solve for m

numIterations = 100;   % we're requested to run 100 iterations of CGLS.

% Below: setting initial algorithm elements as directed
% for the CGLS routine (Aster p155)

mPrev = zeros(size(G,2),1); pPrev = zeros(size(G,1),1);
betaPrev = 0; sPrev = -d; rPrev = G'*sPrev;

resNormStackCGLS = zeros(numIterations+1,1);   % stacks to store res.
modelNormStackCGLS = zeros(numIterations+1,1); % and model norms
resNormStackCGLSnormEq = zeros(numIterations+1,1);

resNormStackCGLS(1,1) = norm(d);
resNormStackCGLSnormEq(1,1) = norm(G'*d);
modelNormStackCGLS(1,1) = 0;

% Below: we'll store models at selected iterations in mStorageMatrix.

mStorageIndices = [10 25 40 55]';
mStorageMatrix = zeros(size(G,2),size(mStorageIndices,1));

for k = 1:numIterations
    pNext = -rPrev+betaPrev*pPrev;                  % steps to left
    alpha = norm(rPrev)^2/((G*pNext)'*(G*pNext));   % correspond w/
    mNext = mPrev+alpha*pNext;                      % steps 1-6 in
    sNext = sPrev+alpha*G*pNext;                    % CGLS algorithm
    rNext = G'*sNext;                               % (Aster, p155)
    betaNext = (rNext'*rNext)/(rPrev'*rPrev);
    
    resNormStackCGLS(k+1,1) = norm(G*mNext-d);  % storing res. and
    modelNormStackCGLS(k+1,1) = norm(mNext);    % model norms each iter.
    resNormStackCGLSnormEq(k+1,1) = norm(G'*G*mNext-G'*d);
    
    % Below: storing new algorithm elements as previous elements for
    % next iteration
    
    mPrev = mNext; pPrev = pNext; betaPrev = betaNext;
    sPrev = sNext; rPrev = rNext;
    
    for counter = 1:size(mStorageIndices,1)
       if k == mStorageIndices(counter,1)      % storing desired
           mStorageMatrix(:,counter) = mNext;  % models (PART E)
       end
    end
end

% Plotting residual norm vs. iteration number...

figure(5); plot((0:100)',log(resNormStackCGLS)/log(10),'.k')
xlabel('iteration number');
ylabel('log10 of residual norm ( log10(|| Gm - d ||) ) ')

% ...and model norm vs. iteration number.

figure(6); plot((0:100)',modelNormStackCGLS,'.k')
xlabel('iteration number'); ylabel('model norm (|| m ||) ')

% For PART E, we plot the model at desired iterations (10,25,40,55).

for k = 1:4
img = reshape(mStorageMatrix(:,k),100,100);
figure(k+6); imagesc(img); colormap('gray')
xlabel('x location in image (out of 100)');
ylabel('y location in image (out of 100)');
end


%%% PART F - using the Steepest Descent algorithm to solve for m

% Below: setting initial algorithm elements as directed
% for the CGLS routine (Aster p155)

mPrev = zeros(size(G,2),1); sPrev = -d; rPrev = G'*sPrev;

resNormStackSD = zeros(numIterations+1,1);   % stacks to store res.
modelNormStackSD = zeros(numIterations+1,1); % and model norms

resNormStackSD(1,1) = norm(d);
modelNormStackSD(1,1) = 0;

for k = 1:numIterations
    
    % CGLS method modified into SD method (pk = -rk; beta not needed)
    
    alpha = (rPrev'*rPrev)/((G*rPrev)'*(G*rPrev));
    mNext = mPrev-alpha*rPrev;    
    sNext = sPrev-alpha*G*rPrev;
    rNext = G'*sNext;
    
    resNormStackSD(k+1,1) = norm(G*mNext-d);  % storing res. and
    modelNormStackSD(k+1,1) = norm(mNext);    % model norms each iter.
    
    % Below: storing new algorithm elements as previous elements for
    % next iteration
    
    mPrev = mNext; sPrev = sNext; rPrev = rNext;
    
end

% Plotting residual norm vs. iteration number (incl. CGLS)...

figure(11); plot((0:100)',log(resNormStackCGLS)/log(10),'.k',...
    (0:100)',log(resNormStackSD)/log(10),'--k')
xlabel('iteration number');
ylabel('log10 of residual norm ( log10(|| Gm - d ||) ) ')
legend('CGLS method','SD method')

% ...and model norm vs. iteration number (incl. CGLS).

figure(12); plot((0:100)',modelNormStackCGLS,'.k',...
    (0:100)',modelNormStackSD,'--k')
xlabel('iteration number'); ylabel('model norm (|| m ||) ')
legend('CGLS method','SD method')


%%% PART G - using MATLAB's pcg() function

GtG = G'*G;  % setting up normal equations for solution with pcg()
Gtd = G'*d;

% using pcg() without preconditioner
[x0,fl0,rr0,it0,rv0] = pcg(GtG,Gtd,1e-4,100);

% using pcg() with incomplete Cholesky preconditioner
L = ichol(GtG,struct('michol','off'));
[x1,fl1,rr1,it1,rv1] = pcg(GtG,Gtd,1e-4,100,L,L');

% Below: comparing convergence of pcg() with and without IC
% preconditioner (along with comparable residuals from Part D's
% CGLS implementation

figure(13);
semilogy(0:it0,rv0/norm(Gtd),'.k');
hold on;
semilogy(0:it1,rv1/norm(Gtd),'-*k');
semilogy(0:100,resNormStackCGLSnormEq/norm(Gtd),'--k');
legend('No preconditioner (pcg())','IC (pcg())',...
    'CGLS method from Part D');
xlabel('iteration number');
ylabel('relative NORMAL residual ( ||GtGm - Gtd || / || Gtd || )');
hold off;

% image output from non-preconditioned pcg() solution

img = reshape(x0,100,100);
figure(14); imagesc(img); colormap('gray'); colorbar
xlabel('x location in image (out of 100)');
ylabel('y location in image (out of 100)');

% image output from preconditioned (IC) pcg() solution

img = reshape(x1,100,100);
figure(15); imagesc(img); colormap('gray'); colorbar
xlabel('x location in image (out of 100)');
ylabel('y location in image (out of 100)');

%%% END PROBLEM 1