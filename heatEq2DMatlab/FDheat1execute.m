N1D = 25;
HOflag = 1;

const = max(1*round(((N1D-1)/49)^2),1);

endtime = 0.11;               % end time of model
k = 0.0001*1/const;
storesteps = (1:100*const:1000*const+1)';

[datastore xnodes ynodes] = ...
    FDheat1(N1D,endtime,k,HOflag,storesteps);