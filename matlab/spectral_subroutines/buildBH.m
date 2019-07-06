function [ T ] = buildT( lambda, adj )

deg=sum(adj);
N=length(deg);
D=sparse(1:N,1:N,deg,N,N,N);

T=(lambda^2-1)*speye(N,N)-lambda*adj+D;

end

