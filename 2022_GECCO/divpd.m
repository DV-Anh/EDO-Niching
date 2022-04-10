function [P,F,e,o,obj,last] = divpd(adjacent,m,l,thres,init_pop)
% PD approach, constrained
% adjacent = adjacent matrix
% m = population size
% l = iteration limit
% opt = optimal tour
n=size(adjacent,1);

P=initgen();
L=zeros(m+1,1);flags=zeros(m+1,1);
for i=1:m
    L(i)=getlength(i);flags(i)=max(L(i)-thres,0);
end

F=[];ed=false;fitness=[];last=0;ind_o=[];
new_dist();

counter=m;
is_better=false;
while counter<l% && (sum(f)>goal || max(f)>min(f)+1)
    P(end,:)=mutation(parentselect());
    L(end)=getlength(m+1);flags(end)=max(L(end)-thres,0);
    update_dist();
    get_fit();
    maxind=getmaxFit_ind();
    if maxind<=m
        if is_better
            last=counter+1;
        end
        P(maxind,:)=P(end,:);
        ed(maxind,:)=false;ed(maxind,ind_o)=true;
        F([maxind,end],:)=F([end,maxind],:);
        F(:,[maxind,end])=F(:,[end,maxind]);
        L(maxind)=L(end);flags(maxind)=flags(end);
    else
    end
    counter=counter+1;
end
P=P(1:m,:);
% F=F(1:m,1:m);
% e=e(1:counter);
e=0;o=0;obj=0;

function c = parentselect()
c=randi(m);
end

function a = mutation(ind)
a=P(ind,:);
fr=randi(n);to=rem(fr+randi(n-3),n)+1;if to>fr to=to-1; else te=to;to=fr-1;fr=te; end;a(fr:to)=a(to:-1:fr); %2-OPT
% fr=randi(n);to=rem(fr+randi(n-4),n)+1;te=a(fr);if to>fr a(fr:to-1)=a(fr+1:to);else a(to+1:fr)=a(to:fr-1); end;a(to)=te; %3-OPT
%fr=randi(n);to=rem(fr+randi(n-5)+1,n)+1;a([to,fr])=a([fr,to]); %4-OPT
end

function get_fit()
fitness=[flags,sort(F,2,'descend')];
end

function d = getmaxFit_ind()
% [~,d]=sortrows(fitness,'desc');d=d(1);
d=m+1;is_better=false;
for I=1:m
    is_eq=true;
    for j=1:m+1
        if fitness(I,j)>fitness(d,j)
            d=I;is_better=true;is_eq=false;
            break;
        elseif fitness(I,j)<fitness(d,j)
            is_eq=false;break;
        end
    end
    if is_eq;d=I;end
end
end

function l = getlength(ind)
y=(P(ind,:)-1)*n;
ind=[P(ind,1:end-1)+y(2:end),P(ind,end)+y(1)];
l=sum(adjacent(ind));
end

function new_dist()
ind=1:m;
cy=(P(ind,:)-1).*(P(ind,:)-2)/2;
ind=max([P(ind,1:end-1)+cy(:,2:end),P(ind,end)+cy(:,1)],[cy(:,1:end-1)+P(ind,2:end),cy(:,end)+P(ind,1)]);
ed=false(m,n*(n-1)/2);
for I=1:m
    ed(I,ind(I,:))=true;
end
F=zeros(m+1);
for I=1:m-1
    for j=I+1:m
        s=sum(ed(I,:)&ed(j,:));F(I,j)=s;F(j,I)=s;
    end
end
end

function update_dist()
cy=(P(end,:)-1).*(P(end,:)-2)/2;
ind_o=max([P(end,1:end-1)+cy(2:end),P(end,end)+cy(1)],[cy(1:end-1)+P(end,2:end),cy(end)+P(end,1)]);
%ed(end,:)=false;ed(end,ind_o)=true;
for I=1:m
    s=sum(ed(I,ind_o));F(I,end)=s;F(end,I)=s;
end
end

function p = initgen()
if size(init_pop,1)>0
    p=[init_pop;zeros(1,n)];
else
    p=zeros(m+1,n);
    for I=1:m
        p(I,:)=randperm(n);
    end
end
end
end