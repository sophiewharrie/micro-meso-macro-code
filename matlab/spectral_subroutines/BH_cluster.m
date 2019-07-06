function [inferred_sigma,overlap]=BH_cluster(adj,sigma,q,deltac)



    opts.maxit=50;
    opts.isreal=1;
    opts.issym=1;  
    opts.tol=1e-3;

    N=length(sigma);
    degrees=sum(adj);
    D=sparse(1:N,1:N,degrees,N,N,N);
    %Regularizing parameter :
    r=sqrt(mean(degrees.^2)/mean(degrees)-1);

    BHplus=(r^2-1)*speye(N,N)-r*adj+D;
    BHminus=(r^2-1)*speye(N,N)+r*adj+D;
    % keyboard

    if(deltac>0)
        [vectorsBHp , eigenvaluesp]=eigs(BHplus,q,'sa',opts);
        vectorsBHm=[];
        eigenvaluesm=[];
    else
        [vectorsBHp , eigenvaluesp]=eigs(BHplus,1,'sa',opts);
        [vectorsBHm , eigenvaluesm]=eigs(BHminus,q-1,'sa',opts);
    end

    eigenvaluesp=diag(eigenvaluesp);
    idplus=find(eigenvaluesp<0);

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



    overlap=0;

    if(q==2)
        if(nCommBH>1)
            inferred_sigma=(vectorsBH(:,2)>0)+1;
            accuracy=sum((inferred_sigma==sigma))/N;
            accuracy=max(accuracy,1-accuracy);
            overlap=2*(accuracy-1/2);

        else
            inferred_sigma=ones(N,1);
            overlap=0
        end
    end



    if(q>2)

        inferred_sigma=kmeans(vectorsBH,q,'Replicates',1);

        % We maximize the overlap over all permutations of the inferred communities using a hungarian
        % algorithm. 

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
