function mask=GMM(dis,m)
% Greedy heuristic for Maximimzing Minimum distance
if m<2
    error('Must select more than 2 elements. Cannot select %d elements.',m)
end
n=size(dis,1);
if m>n
    error('Cannot select %d elements from %d elements.',m,n)
end
mask=false(n,1);
%Choose 2 most distant elements
[maxa,a]=max(dis);[~,b]=max(maxa);mask(b(1))=true;mask(a(1,b(1)))=true;
for i=3:m
    ind=find(~mask);
    maxa=sort(dis(mask,ind')); %smallest distance from members to non-members with tie-breakers
    [~,a]=sortrows(maxa','desc'); %greedy
    mask(ind(a(1)))=true;
end
end