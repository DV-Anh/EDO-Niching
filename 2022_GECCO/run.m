t=[.05,.1,.2]; % Thresholds
L=30; % Repetitions

% Load instances
load('tsp_instances.mat');

% Variable holding run results
results=cell(numel(tsp),numel(t),L);

% Go through instances
for i=1:numel(tsp)
    name=tsp{i}.name;
    g=tsp{i}.graph;OPT=tsp{i}.opt;
    n=size(g,1);
    u=floor(n/4);
    opt=getlength(g,OPT);
    init_pop=repmat(OPT,u,1);
    iter=floor(4*(u)*(n^2));
    for h=1:numel(t)
        th=t(h);
        thres=opt*(1+th);
        for k=1:L
			% Setup
            results{i,h,k}.name=name;
            results{i,h,k}.u=u;results{i,h,k}.iter=iter;
            results{i,h,k}.thres=thres;results{i,h,k}.thres_ratio=th;
			
			% Run baselines
			[P,~,~,~,~,e]=dived(g,results{i,h,k}.u,results{i,h,k}.iter,results{i,h,k}.thres,init_pop);
			results{i,h,k}.ED.pop=P;results{i,h,k}.ED.iter=e;
			[P,~,~,~,~,e]=divpd(g,results{i,h,k}.u,results{i,h,k}.iter,results{i,h,k}.thres,init_pop);
			results{i,h,k}.PD.pop=P;results{i,h,k}.PD.iter=e;
			
			% Run modified NMA
			[P,e]=div_tsp_p1(g,results{i,h,k}.thres,results{i,h,k}.u,results{i,h,k}.u*3,results{i,h,k}.iter);
			results{i,h,k}.niche.pop=P;results{i,h,k}.niche.iter=e;
			l=getlength(g,P);P=P(l<=results{i,h,k}.thres,:);
			M=getdist_tsp(P);M=GMM(M,results{i,h,k}.u);X=P(M,:);
			
			% Run phase 2
			[P,~,~,~,~,e]=dived(g,results{i,h,k}.u,results{i,h,k}.iter-results{i,h,k}.niche.iter,results{i,h,k}.thres,X);
			results{i,h,k}.nicheED.pop=P;results{i,h,k}.nicheED.iter=e+results{i,h,k}.niche.iter;
			[P,~,~,~,~,e]=divpd(g,results{i,h,k}.u,results{i,h,k}.iter-results{i,h,k}.niche.iter,results{i,h,k}.thres,X);
			results{i,h,k}.nichePD.pop=P;results{i,h,k}.nichePD.iter=e+results{i,h,k}.niche.iter;
        end
    end
end
save(['results_tsplib_niche','.mat'], 'results', '-v7.3','-nocompression');