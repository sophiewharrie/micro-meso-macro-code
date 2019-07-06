function [inferred_sigma,overlap]=BH_cluster_real_world(adj,sigma,q)

    opts.maxit=500;
    opts.isreal=1;
    opts.issym=1;  
    opts.tol=1e-10;

    deg=sum(adj);
    N=length(deg);
    guessForFirstEigen=mean(deg.^2)/mean(deg)-1; 

    errtol=1e-10;
    maxIter=10;
    [ x, rhoB ] = iterateSlp( guessForFirstEigen, adj , errtol , maxIter );

    r=sqrt(rhoB);
    


    maxComm=floor(sqrt(N));



    BHplus=buildBH(r, adj );
    BHminus=buildBH(-r, adj );

    [vectorsBHp , eigenvaluesp]=eigs(BHplus,maxComm,'sa',opts);

    eigenvaluesp=diag(eigenvaluesp);
    idplus=find(eigenvaluesp<0);


    [vectorsBHm , eigenvaluesm]=eigs(BHminus,maxComm,'sa',opts);



    eigenvaluesm=diag(eigenvaluesm);
    idminus=find(eigenvaluesm<0);

    nCommBH=length(idplus)+length(idminus);

    eigenvalues=[eigenvaluesp;eigenvaluesm];
    vectorsBH=[vectorsBHp,vectorsBHm];
    [val,id]=sort(eigenvalues,'ascend');



    if(nCommBH<=q)
        vectorsBH=vectorsBH(:,id(1:nCommBH));
    else
        vectorsBH=vectorsBH(:,id(1:q));
    end



% Normalizing all eigenvectors to 1
vectorsBH= bsxfun(@times,1./sqrt(sum(vectorsBH.^2)),vectorsBH);

eigenvalues=eigenvalues(id);


overlapBH=0;

if(q==2)
    inferred_sigma=(vectorsBH(:,2)>0)+1;
    accuracy=sum((inferred_sigma==sigma))/N;
    accuracy=max(accuracy,1-accuracy);
    overlapBH=2*(accuracy-1/2);
end



if(q>2)
    inferred_sigma=kmeans(vectorsBH,q,'Replicates',1000);

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

    overlapBH=(accuracy-1/q)/(1-1/q);


end


%fprintf('%d communities were actually found using clustering with the Bethe Hessian\n',nCommBH);
%fprintf('The overlap computed with %d communities is %f\n',q,overlapBH);


overlap=overlapBH;


end


