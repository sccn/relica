function X = relica_bootstrap(X,ltrial)

if nargin == 1
    N=size(X,2);
    index=round(rand(N,1)*N+.5);
    X=X(:,index);
else
    N=size(X,2)/ltrial;
    T = reshape(X,size(X,1),ltrial,N);
    index=round(rand(N,1)*N+.5);
    T=T(:,:,index);
    X = reshape(T,size(T,1),prod(size(T(1,:))));
end