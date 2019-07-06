function [ Tprime ] = buildTprime( lambda, adj )

deg=sum(adj);
N=length(deg);
D=sparse(1:N,1:N,deg,N,N,N);

Tprime=2*lambda*speye(N,N)-adj;

end

