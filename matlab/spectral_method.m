function inferred_sigma = spectral_method(filename)
    % Spectral community detection method adapted from code written
    % by Alaa Saade (available at http://mode_net.krzakala.org/).
 	% INPUTS :
 	% name is a string containing the name of the network to load, e.g. 'dolphins.gml'.
    % OUTPUTS :
    % inferred_sigma is the community assignement inferred by the spectral clustering method 
    warning('off','all')
    path(path,'spectral_subroutines');

    % read in network
    [E,sigma]=read_gml(filename);

    N=max(max(E));
    q=N; % edit
    nnz=2*length(E)+N;
    si=nnz-N;
    i0=E(:,1);
    j0=E(:,2);
    i=[i0;j0];
    j=[j0;i0];
    
    adj=sparse(i,j,ones(2*length(E),1),N,N);
  
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

    deltac = 1; % edit
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

    q=length(idplus)+length(idminus);
    
    % workaround for issue with code crashing when spectral method detects only 1 community
    if(q==1)
        inferred_sigma = ones(N,1);
        overlap = 0.0; % set it to anything because it's not used
        return
    end
    
    warning('off','MATLAB:eigs:NoEigsConverged');
    warning('off','MATLAB:eigs:NotAllEigsConverged');
    warning('off','stats:kmeans:FailedToConvergeRep');
    warning('off','stats:kmeans:EmptyClusterInBatchUpdate');
    
    [inferred_sigma,overlap]=BH_cluster_real_world(adj,sigma,q);
    
end
    