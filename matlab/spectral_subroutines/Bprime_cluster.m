function [inferred_sigma,overlap]=Bprime_cluster(adj,sigma,q)
    % Does spectral clustering based on Bprime, and compute the overlap with the true community memberships.
    
    
    opts.maxit=50;
    opts.isreal=1;  
    opts.tol=1e-3;

    N=length(sigma);
    degrees=sum(adj);
    D=sparse(1:N,1:N,degrees,N,N,N);
    Bprime=[sparse(N,N),(D-speye(N,N)) ; -speye(N,N),adj];

    
    [vectorsBprime , eigenvalues]=eigs(Bprime,q,'lm',opts);


    eigenvalues=diag(eigenvalues);
    idReal=find(abs(imag(eigenvalues))<0.0001);
    eigenvalues=eigenvalues(idReal);
    vectorsBprime=vectorsBprime(:,idReal);
    [val,id]=sort(abs(eigenvalues),'descend');
    eigenvalues=eigenvalues(id);
    rhoB=eigenvalues(1);
    nCommBprime=sum((abs(eigenvalues)>sqrt(rhoB)));

    if(nCommBprime<q)
        vectorsBprime=vectorsBprime(:,id(1:nCommBprime));
    else
        vectorsBprime=vectorsBprime(:,id(1:q));
    end



    vectorsBprime=vectorsBprime(N+1:2*N,:); 
    %Column normalization to 1
    vectorsBprime= bsxfun(@times,1./sqrt(sum(vectorsBprime.^2)),vectorsBprime); 




    overlap=0;

    if(q==2)
        if(nCommBprime>1)
            inferred_sigma=(vectorsBprime(:,2)>0)+1;
            accuracy=sum((inferred_sigma==sigma))/N;
            accuracy=max(accuracy,1-accuracy);
            overlap=2*(accuracy-1/2);

            
        else
            inferred_sigma=ones(N,1);
            overlap=0;
        end
    end


    if(q>2)
        
        inferred_sigma=kmeans(vectorsBprime,q,'Replicates',1);

        
        costMat=zeros(q,q);
        for a=1:q
            for b=1:q
                id=find(inferred_sigma==a);
                costMat(a,b)=-(sum(sigma(id)==b))^2;
            end
        end
        [permutation,cost] = munkres(costMat);


        inferred_sigma=permutation(inferred_sigma(:))';

        accuracy=sum((inferred_sigma==sigma))/N;

        overlap=(accuracy-1/q)/(1-1/q);


    end
    
    % keyboard

end