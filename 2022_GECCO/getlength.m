function l = getlength(ad,X)
% Return length (cost) of a list of tsp tours given distance matrix between vertices
[m,n]=size(X);
l=zeros(m,1);
for i=1:m
    x=X(i,:);
    y=(x-1)*n;
    ind=[x(1:end-1)+y(2:end),x(end)+y(1)];
    l(i)=sum(ad(ind));
end
end