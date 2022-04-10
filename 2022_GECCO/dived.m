function [P,F,e,o,obj,last] = dived(adjacent,m,l,thres,init_pop)
% ED approach, constrained
% adjacent = adjacent matrix
% m = population size
% l = iteration limit
% thres = alpha threshold
n=size(adjacent,1);
P=initgen();
L=zeros(m+1,1);flags=zeros(m+1,1);
for i=1:m
    L(i)=getlength(i);flags(i)=max(L(i)-thres,0);
end

last=0;
F=zeros(n,n);Ind_a=zeros(m+1,n);init_freq();
fitness=zeros(m+1,n+1);

counter=m;
is_better=false;
while counter<l% && max(f)>min(f)+1
    P(end,:)=mutation(parentselect());
    L(end)=getlength(m+1);flags(end)=max(L(end)-thres,0);
    add_freq();
    get_fit();
    maxind=getmaxFit_ind();
    sub_freq();
    if maxind<=m
        if is_better
            last=counter+1;
        end
        P(maxind,:)=P(end,:);Ind_a(maxind,:)=Ind_a(end,:);
        L(maxind)=L(end);flags(maxind)=flags(end);
    else
    end
    counter=counter+1;
end
P=P(1:m,:);
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

function init_freq()
ind=1:m;
cy=(P(ind,:)-1)*n;
Ind=[P(ind,1:end-1)+cy(:,2:end),P(ind,end)+cy(:,1),cy(:,1:end-1)+P(ind,2:end),cy(:,end)+P(ind,1)];
F=zeros(n);
for I=1:m
    F(Ind(I,:))=F(Ind(I,:))+1;
end
Ind_a(ind,:)=Ind(:,1:n);
end

function add_freq()
ind=m+1;
cy=(P(ind,:)-1)*n;
Ind=[P(ind,1:end-1)+cy(:,2:end),P(ind,end)+cy(:,1),cy(:,1:end-1)+P(ind,2:end),cy(:,end)+P(ind,1)];
F(Ind)=F(Ind)+1;
Ind_a(end,:)=Ind(1:n);
end

function sub_freq()
ind=maxind;
cy=(P(ind,:)-1)*n;
Ind=[P(ind,1:end-1)+cy(:,2:end),P(ind,end)+cy(:,1),cy(:,1:end-1)+P(ind,2:end),cy(:,end)+P(ind,1)];
F(Ind)=F(Ind)-1;
end

function get_fit()
fitness=[flags,sort(F(Ind_a),2,'desc')];
end

function d = getmaxFit_ind()
% [~,d]=sortrows(fitness,'desc');d=d(1);
d=m+1;is_better=false;
for I=1:m
    is_eq=true;
    for j=1:n+1
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