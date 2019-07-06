function [ x, lambda ] = iterateSlp( lambda0, adj, errtol , maxIter )
%Iterates SLP starting from lambda0 and tries to find the closest nonlinear
%eigenvalue
err=1;
iter=0;
lambda = lambda0;
opts.issym=1;
 while err > errtol && iter < maxIter ,
        iter=iter+1;
        BH = buildBH( lambda, adj );
        BHprime = buildBHprime( lambda, adj );  
        [x, mu]=eigs(BH,BHprime,1,'sm',opts);
        err=abs(mu);
        lambda=lambda-mu;
 end




