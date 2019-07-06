function [inferred_sigma,overlap]=Lap_cluster(adj,sigma,q,deltac)
 opts.maxit=50;
 opts.isreal=1;
 opts.issym=1;  
 opts.tol=1e-3;

 N=length(sigma);

 Nold=N;

 [A,p] = largest_component(adj);    
 N=length(A(:,1));

 degrees=sum(A);
 D=sparse(1:N,1:N,degrees,N,N,N);


 Dminus_half=spfun(@(x)1./sqrt(x),D);
 Lap=Dminus_half*(D-A)*Dminus_half;


 if(deltac>0)
    [vectorsLap , eigenvalues]=eigs(Lap,q,'sa',opts);
    eigenvalues=diag(eigenvalues);
    [eigenvalues, id]=sort(eigenvalues,'ascend');

    vectorsLap=vectorsLap(:,id);
else
    [vectorsLap , eigenvalues]=eigs(Lap,q,'la',opts);
    eigenvalues=diag(eigenvalues);
    [eigenvalues, id]=sort(eigenvalues,'descend');

    vectorsLap=vectorsLap(:,id);
end


vectorsLap= bsxfun(@times,1./sqrt(sum(vectorsLap.^2)),vectorsLap);


overlap=0;

if(q==2)

    inferred_sigmaGC=(vectorsLap(:,2)>0)+1;

    inferred_sigma=zeros(Nold,1);
    inferred_sigma(p==1)=inferred_sigmaGC;
    %Assigning randomly the nodes that are not in the giant component
    inferred_sigma(p==0)=randi(2,Nold-N,1);

    accuracy=sum((inferred_sigma==sigma))/Nold;
    accuracy=max(accuracy,1-accuracy);
    overlap=2*(accuracy-1/2);

end



if(q>2)

    inferred_sigmaGC=kmeans(vectorsLap,q,'Replicates',1);

    inferred_sigma=zeros(Nold,1);
    inferred_sigma(p==1)=inferred_sigmaGC;
    %Assigning randomly the nodes that are not in the giant component
    inferred_sigma(p==0)=randi(q,Nold-N,1);


    costMat=zeros(q,q);
    for a=1:q
        for b=1:q
            id=find(inferred_sigma==a);
            costMat(a,b)=-(sum(sigma(id)==b))^2;
        end
    end
    [permutation,cost] = munkres(costMat);

    inferred_sigma=permutation(inferred_sigma(:))';



    accuracy=sum((inferred_sigma==sigma))/Nold;
    
    overlap=(accuracy-1/q)/(1-1/q);

end



end







