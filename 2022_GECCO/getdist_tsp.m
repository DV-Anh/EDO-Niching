function f=getdist_tsp(X)
% Return edge distance matrix from a list of tsp tours
[m,n]=size(X);
mask=false(m,n*(n-1)/2);
for i=1:m
    x=sort([X(i,:);X(i,2:end),X(i,1)]);mask(i,(x(2,:)-1).*(x(2,:)-2)/2+x(1,:))=true;
end
f=zeros(m,m);
for i=1:m-1
    for j=i+1:m
        f(i,j)=n-sum(mask(i,:)&mask(j,:));f(j,i)=f(i,j);
    end
end
end