function [p1_pops,eval,gen,sat_count,divs]=div_tsp_p1(adj,thres,output_size,init_pop,budget)
n=size(adj,1);
edge_no=n*(n-1)/2;
rem=0;eval=0;
if size(init_pop,2)==n
    P=init_pop;
    m=size(P,1);
else
    m=max(init_pop,output_size);
    P=initialize();
end
m=size(P,1);
n_min=4;n_max=12;%min(12,floor(m/output_size*1.5));
if mod(n_max,2);n_max=n_max+1;end %min/max size of neighborhood
beta=ceil(m/n_min);
crossover_rate=.9;mutation_rate=.01;localsearch_bud=ceil(n*(n-3));
fit=getlength(adj,P);
eval=eval+numel(fit);
sg=zeros(m,1);not_local_opt=true(m,1);
CES=zeros(edge_no,1); %critical edge set
no_2opt=n*(n-3)/2;ind_2opt=2*ones(no_2opt,2);
ind_2opt(1:n-3,2)=3:(n-1);k=n-2;
for i=3:(n-1);ind_2opt(k:(k+n-i-1),2)=(i+1):n;ind_2opt(k:(k+n-i-1),1)=i;k=k+n-i;end
sat_count=zeros(1,budget);sat_count(1:eval)=sum(fit<=thres);gen_eval=eval;

% Store running diversity (optional)
% divs=zeros(1,budget);divs(1:eval)=diversity(P);
no_skip=true;gen=1;
while eval<budget
    Parent=P;Parent_fit=fit;P_not_local_opt=not_local_opt;
    new_off=false(m,1);
    mate_pool=neighborStrat();
    group_no=numel(mate_pool);
	
	% Diversity enhancement subroutine
	divEnhance();
	
	% Update Critical Edge Set, unused by default
    %updateCES();
	
	% Creating offspring
    generate_offspring();
    sat_count(gen_eval+1:eval-1)=sum(Parent_fit<=thres);
%     divs(gen_eval+1:eval-1)=diversity(Parent);
    if no_skip;survivalSelect();else;break;end
    P=Parent;fit=Parent_fit;not_local_opt=P_not_local_opt;
    gen=gen+1;
	
    sat_count(eval)=sum(fit<=thres);
%     divs(eval)=diversity(P);
    gen_eval=eval;
	
	% If enough satisfactory solutions are found, terminates
    if sat_count(eval)>=output_size;break;end
end
p1_pops=P;
%preserved=preserve();
function p=initialize()
    p=zeros(m,n);
	
	% Random initialization, unused by default
%     for I=1:m
%         p(I,:)=randperm(n);
%     end
	
	% Half-greedy initialization
    st=floor(n/2);add=0;
    for I=1:m
        p(I,1:st)=randperm(n,st);
        mask=true(1,n);mask(p(I,1:st))=false;
        for j=st+1:n-1
            Ind=find(mask);
            [~,ind]=min(adj(p(I,j-1),Ind));
            add=add+sum(mask);
            ind=Ind(ind(1));p(I,j)=ind;mask(ind)=false;
        end
        left=find(mask);p(I,end)=left(1);
    end
    rem=rem+mod(-add,n);
    eval=eval+ceil(add/n)-floor(rem/n);
    rem=mod(rem,n);
end

% Post-process function from original NMA, unused by default
function p=preserve()
    l=min(fit);p=false(m,1);
    for I=1:group_no
        ind=mate_pool{I};group_size=numel(ind);
        [L,s]=sort(fit(ind));p(s(1))=true;
        for J=1:group_size
            if L(J)<=l;p(ind(s(J)))=true;else
                if J==1&&fit(ind(s(J)))<=1.1*l
                    pre=find(p);p_no=numel(pre);
                    if p_no<=0;p(ind(s(J)))=true;continue;end
                    D=zeros(p_no,1);
                    for K=1:p_no
                        D(K)=distance(pre(K),ind(s(J)));
                    end
                    MIN=min(D);if MIN/n>.2;p(ind(s(J)))=true;end
                end
            end
        end
    end
end

% Group forming function
function pool=neighborStrat()
    [fit_sorted,fit_rank]=sort(fit(1:m));fit_rank_inv(fit_rank)=1:m;
    length_min=fit_sorted(1);
    length_beta=fit_sorted(beta);
    tmpP=true(m,1);pool=cell(ceil(m/n_min),1);group_count=1;tmpP_size=m;
    while tmpP_size>0
        tmpP_ranks=fit_rank(tmpP);leader=tmpP_ranks(1); %by index on fit_rank
        if tmpP_size>=n_max
            gamma=min((fit(leader)-length_min)/(length_beta-length_min),1);
            M=floor(gamma*(n_max-n_min)+n_min);
            if mod(M,2);M=M+1;end
        else
            M=tmpP_size;
        end
        %Calculate sharing dist
        share_dist=Inf(tmpP_size-1,1);
        for I=1:tmpP_size-1
%             if fit(tmpP_ranks(I+1))>thres %Avoid putting satisfactory sols into same group
            share_dist(I)=distance(tmpP_ranks(I+1),leader);
%             end
        end
        [~,share_dist_sortInd]=sort(share_dist);
        share_dist_sortInd=share_dist_sortInd(1:M-1);share_dist_sortInd=share_dist_sortInd+1;
        neighbor_group=[leader;tmpP_ranks(share_dist_sortInd(randperm(M-1)))]; %shuffle & keep leader at top
        pool{group_count}=neighbor_group;group_count=group_count+1;
        tmpP(fit_rank_inv(neighbor_group))=false;tmpP_size=tmpP_size-M;
    end
    pool=pool(1:group_count-1);
end

% Diversity enhancement subroutine
function divEnhance()
    for I=1:group_no
        group=mate_pool{I};group_size=numel(group);
		% Considering non-leaders only
        for j=2:group_size
			% If distance is 0 (duplicate), perform migration
            if distance(group(1),group(j))<=0
				% Sample group to migrate, must have more than 1 solutions
                ing_to_migrate=mod(I+randi(group_no-1)-1,group_no)+1;
                group_to_migrate=mate_pool{ing_to_migrate};
                while numel(group_to_migrate)<=1
                    ing_to_migrate=mod(I+randi(group_no-1)-1,group_no)+1;
                    group_to_migrate=mate_pool{ing_to_migrate};
                end
				
				% Sample non-leader position (leader is stored at 1)
                ind_to_migrate=randi(numel(group_to_migrate)-1)+1;
				
				% Swap groups
                tmp=group_to_migrate(ind_to_migrate);
                group_to_migrate(ind_to_migrate)=group(j);
                mate_pool{ing_to_migrate}=group_to_migrate;
                group(j)=tmp;
            end
        end
    end
end

function updateCES()
    CES(CES>n)=0;
    CES_pos=CES>0;CES(CES_pos)=CES(CES_pos)+1;
    common_leader_edges=true(edge_no,1);count=0;
    for I=1:group_no
        l=mate_pool{I}(1);
        if sg(l)<n;continue;end
        x_b=false(edge_no,1);
        x=[P(l,:);P(l,2:end),P(l,1)];
        x=sort(x);x_b((x(2,:)-1).*(x(2,:)-2)/2+x(1,:))=true;
        common_leader_edges=common_leader_edges&x_b;count=count+1;
    end
    if sum(common_leader_edges)<=0||count<=0;return;end
    CES(common_leader_edges)=1;CES_pos=CES_pos|common_leader_edges;
    for I=1:n %check conflict: more than 2 edges adjacent to a vertex
        edges=[1:(I-1),I*ones(1,n-I)];edges=edges.*(edges-1)/2+[I*ones(1,I-1),I+1:n]-1;
        edges_flag=false(edge_no);edges_flag(edges)=true;
        edges_flag=find(edges_flag&CES_pos);s=numel(edges_flag);
        if s>=2
            [~,added_order_sorted]=sort(CES(edges_flag),'descend'); %oldest edges removed
            remove_ind=edges_flag(added_order_sorted(1:s-2));
            CES(remove_ind)=0;CES_pos(remove_ind)=false;
        end
    end
end

function generate_offspring()
    for I=1:group_no
        group=mate_pool{I};
        group_size=numel(group);
        do_crossover=find(rand(floor(group_size/2),1)<=crossover_rate);
        do_crossover_ind=do_crossover*2;
        new_off(group(do_crossover_ind))=true;
        new_off(group(do_crossover_ind-1))=true;
        %Crossover
        for jj=1:numel(do_crossover) %CEA-PMX
            j=do_crossover(jj)*2-1;
            off=crossover(group(j),group(j+1));
            P([group(j),group(j+1)],:)=off;
        end
    end
    %Mutation
    do_mutate=find(rand(m,1)<=mutation_rate);
    P(do_mutate,:)=mutation(do_mutate);
    new_off(do_mutate)=true;
    %Local search
    new_no=sum(new_off);if eval+new_no>budget;no_skip=false;return;end
    fit(new_off)=getlength(adj,P(new_off,:));eval=eval+new_no;
    not_local_opt(new_off)=true;
    %assert(max(abs(fit-getlength(adj,P))./fit)<0.00001,'4')
    for I=1:group_no %do local search on unsatisfactory solutions only
        ind=mate_pool{I};group_size=numel(ind);
        unsat=ind(fit(ind)>thres&not_local_opt(ind));unsat_size=numel(unsat);
        if unsat_size<=0;continue;end
        [~,do_localsearch]=sort(fit(unsat));
		do_localsearch=ind(do_localsearch(1:min(ceil(group_size/2),unsat_size)));
        localSearch(do_localsearch);
    end
end

function mutants=mutation(ind)
    mutants=zeros(numel(ind),n);
    for mi=1:numel(ind)
        x=P(ind(mi),:);tx=[x;x(2:end),x(1)];tx=sort(tx);tx=(tx(2,:)-1).*(tx(2,:)-2)/2+tx(1,:);
        allowed_select=true(n,1);ti=1;
        while ti<=n
            if CES(tx(ti)) %preserve critical edges
                allowed_select(ti)=false;
                allowed_select(mod(ti,n)+1)=false;
                ti=ti+2;
            else
                ti=ti+1;
            end
        end
        allowed_ind=find(allowed_select);
        if numel(allowed_ind)<=0;continue;end
        s1=allowed_ind(randi(numel(allowed_ind)));
        allowed_ind=allowed_ind(mod(allowed_ind-s1,n)>2&mod(s1-allowed_ind,n)>2); %non-adjacent picks (4-opt)
        if numel(allowed_ind)<=0;continue;end
        s2=allowed_ind(randi(numel(allowed_ind)));
        x([s1,s2])=x([s2,s1]);mutants(mi,:)=x;
    end
end

function off=crossover(a,b)
    off=zeros(2,n);off_b=false(2,n);rand_pos=sort(randperm(n,2));rand_r=rand_pos(1):rand_pos(2);
    P(a,:)=circshift(P(a,:),mod(randi(n),n));
    off(:,rand_r)=P([a,b],rand_r); %random segment
    off_b(1,P(a,rand_r))=true;off_b(2,P(b,rand_r))=true;
    inv=zeros(2,n);inv(1,P(a,:))=1:n;inv(2,P(b,:))=1:n;
    for K=1:n
        if off(1,K)<=0
            insert=P(b,K);
            while(off_b(1,insert))
                insert=P(b,inv(1,insert));
            end
            off(1,K)=insert;off_b(1,insert)=true;
        end
        if off(2,K)<=0
            insert=P(a,K);
            while(off_b(2,insert))
                insert=P(a,inv(2,insert));
            end
            off(2,K)=insert;off_b(2,insert)=true;
        end
    end
end

function localSearch(ind)
    for il=1:numel(ind)
        I=ind(il);
        x=[P(I,:);P(I,2:end),P(I,1)];x=sort(x);x=(x(2,:)-1).*(x(2,:)-2)/2+x(1,:);
        crits=CES(x);allowed_ind=crits<=0;rem_bud=min(localsearch_bud,floor(n*(budget-eval)/4+rem/4));count=0;
        while count<rem_bud %first improvement local search
            is_stop=false;is_stuck=true;
            ind_2opt=ind_2opt(randperm(no_2opt),:); %scan in random order
            for j=1:no_2opt
                if count>=rem_bud||fit(I)<=thres
                    is_stop=true;is_stuck=false;break;
                end
                jl=ind_2opt(j,1)-1;jj=ind_2opt(j,2);
                if ~(allowed_ind(jl)&&allowed_ind(jj));continue;end
                next=mod(jj,n)+1;
                e=[P(I,jl),P(I,jl+1),P(I,jl),P(I,jj);P(I,jj),P(I,next),P(I,jl+1),P(I,next)];
                e=sort(e);ed=(e(2,:)-1).*(e(2,:)-2)/2+e(1,:);
                plus=adj(e(1,1),e(2,1))+adj(e(1,2),e(2,2))-adj(e(1,3),e(2,3))-adj(e(1,4),e(2,4));
                count=count+1;
                if plus<0
                    is_stuck=false;
                    fw=jl+1:jj;P(I,fw)=P(I,flip(fw));fit(I)=fit(I)+plus;
                    allowed_ind(fw(1:end-1))=allowed_ind(fw(end-1:-1:1));
                    allowed_ind(jj)=CES(ed(2))<=0;allowed_ind(jl)=CES(ed(1))<=0;
                end
            end
            if is_stuck;not_local_opt(I)=false;break;end
            if is_stop;break;end
        end
        add=count*4;rem=rem+mod(-add,n);
        eval=eval+ceil(add/n)-floor(rem/n);
        rem=mod(rem,n);
    end
end

function survivalSelect()
    %prepare edge mask
    mask_o=false(m,edge_no);mask_p=false(m,edge_no);
    for I=1:m
        ti=[P(I,:);P(I,2:end),P(I,1)];ti=sort(ti);mask_o(I,(ti(2,:)-1).*(ti(2,:)-2)/2+ti(1,:))=true;
        ti=[Parent(I,:);Parent(I,2:end),Parent(I,1)];ti=sort(ti);mask_p(I,(ti(2,:)-1).*(ti(2,:)-2)/2+ti(1,:))=true;
    end
    %compete with closest parent in the same group
    replaced_off_ind=zeros(m,1);
    for I=1:group_no
        ind=mate_pool{I};group_size=numel(ind);[~,s_ind]=sort(fit(ind));
        %Only replace unsatisfactory parents, unused by default
%         bad=find(Parent_fit(ind)>thres);no_bad=numel(bad);if no_bad<=0;continue;end
        for j=1:group_size
            off=ind(s_ind(j));
            sim=zeros(group_size,2);
            for h=1:group_size
                sim(h,:)=Parent_fit(ind(h));
            end
            [~,par]=sortrows(sim,'descend');
            par=ind(par(1));
            if Parent_fit(par)>fit(off)
                Parent_fit(par)=fit(off);
                replaced_off_ind(par)=off;
            end
        end
    end
    replaced=find(replaced_off_ind>=1);sg=sg+1;
    for I=1:numel(replaced)
        r=replaced(I);sg(r)=0; %reset stagnation
        Parent(r,:)=P(replaced_off_ind(r),:);
        P_not_local_opt(r)=not_local_opt(replaced_off_ind(r));
    end
    %assert(max(abs(Parent_fit-getlength(adj,Parent))./Parent_fit)<0.00001,'5')
end

function d=distance(i,j)
    ti_b=false(edge_no,1);
    ti=[P(i,:);P(i,2:end),P(i,1)];ti=sort(ti);ti_b((ti(2,:)-1).*(ti(2,:)-2)/2+ti(1,:))=true;
    tj_b=false(edge_no,1);
    tj=[P(j,:);P(j,2:end),P(j,1)];tj=sort(tj);tj_b((tj(2,:)-1).*(tj(2,:)-2)/2+tj(1,:))=true;
    d=n-sum(ti_b&tj_b);
end
end