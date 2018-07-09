function intStack = intStackSetup(intLocs,xgrid,ufxtK,bigSize,N)

    numInterfaces = size(intLocs,1);
    
    intStack = zeros(numInterfaces,5);
    
    for m = 2:N
        for intCounter = 1:numInterfaces
            if xgrid(m,1) > intLocs(intCounter,1) && ...
                    xgrid(m-1,1) < intLocs(intCounter,1)
                intStack(intCounter,1) = intLocs(intCounter,1);
                intStack(intCounter,2) = 1/(ufxtK(2*bigSize+m-1,1));
                intStack(intCounter,3) = 1/(ufxtK(2*bigSize+m,1));
                intStack(intCounter,4) = (ufxtK(2*bigSize+m-1,1)*...
                    ufxtK(2*bigSize+m-1+N,1))^0.5;
                intStack(intCounter,5) = (ufxtK(2*bigSize+m,1)*...
                    ufxtK(2*bigSize+m+N,1))^0.5;
            end
        end
    end

end

