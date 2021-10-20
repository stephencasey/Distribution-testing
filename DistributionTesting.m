clearvars -except mapsize mapn xmin
clear

% site values: N=1 C1=2 C2=3 C3=4 C4=5 C5=6 C6=7 C7=8 C8=8 C9=9 C10=10
% Transiton=11 Impounded=12 PSU4=13 10x10=14

excelfilename = 'C:\Users\thorn\OneDrive\Desktop\RSL_Github\1995\DistributionFitting.xlsx';

delete(excelfilename)
mapn=33;            % Number of maps
ntrials=1000;       % Number of random sets for P-score
mapsize=6000;       % Map size
xmin=100;           % Hard limit for xmin 

tstart2=tic;

for site=1:mapn
    
    filename= ['C:\Users\thorn\OneDrive\Desktop\RSL_Github\1995\OutputData\dist_array' int2str(site) '.mat'];
    
    load(filename)
    tstart=tic;
    mapelements=mapsize^2;

    % data_n values:  areas=1 ridge_width=2 ridge_length =3 slough_width=4
    % slough_length=5 TI_width=6 TI_length=7 no_TI_width=8 no_TI_length=9
    
    % dist_n values:  pareto=1 GP=2 exp(cont.)=3 logn=4 invgauss=5 gamma=6
    % powexp=7 weibull=8 log-wbl=9 powhyp=10 exp(discrete)=11
    
    for data_n=[1]
        disp(site);
        disp(data_n);
        dist_n=0;
        
        % Set x=data and truncate below xmin
        x=dist_array(:,data_n);
        if data_n==2
            x=sloughareas;
        end
        
        if max(x) == 0
            disp('No data')
            break
        end
        
        xbound=xmin-1;
        x=x(x>=xmin);
        nx=length(x);
        unique_x=unique(x);
        
        %% Check for pareto lower bound
        
        if data_n==1
            x1=unique_x;
            x2=sort(x);
            nx=length(x1);
            exp_power1=zeros(nx-10,1);
            ks_power1=zeros(nx-10,1);
            
            for n=1:nx-10
                xm=x1(n);
                x2(x2 <= xm) = [];
                exp_power1(n)=1+(mean(log(x2./xm)))^-1;
                cdf_model_power=unique(1-(x2./xm).^(-exp_power1(n)+1));
                cdf_data=ecdfcalc(x2);
                ks_power1(n)=max(abs(cdf_data-cdf_model_power));
            end
            
            [ks_power nmin]=min(ks_power1);
            xmin_pow=x1(nmin);
            xmin=xmin^2;
            xbound=xmin-1;
            x=x(x>=xmin);
            nx=length(x);
            unique_x=unique(x);
            
        end
        
        
        %% Calculate Empirical CDF
        emp_x=ecdfcalc(x);
        
        aic=zeros(1,9);
        pscore1=zeros(1,9);
        pscore2=zeros(1,9);
        like=zeros(nx,9);
        rsquared=zeros(1,10);
        x_diff1=zeros(length(unique_x),10);
        x_diff2=zeros(length(unique_x),10);
        
        emp_x_all=zeros(length(x),1);
        for n=1:length(x)
            emp_x_all(n)=emp_x(x(n)==unique_x);
        end;
        
        tstart=tic;
        %% 1.pareto
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=1+nx*(1/(sum(log(x./xbound))));
        pareto_pdf=@(x,alpha) ((alpha-1)*xbound^(alpha-1))*x.^-alpha;
        flike=@(fit) -sum(log(pareto_pdf(x,fit)));
        pareto_fit=fminsearch(flike, fit1)
        pareto_cdf=@(x,pareto_fit) 1-(x./xbound).^(-pareto_fit+1);
        
        pareto_ccdf= @(x,pareto_fit) (1-pareto_cdf(x,pareto_fit))./(1-pareto_cdf(xmin,pareto_fit));
        
        model_x=pareto_cdf(unique_x,pareto_fit);
        ksstat1=max(abs(emp_x-model_x));
        ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
        
        x_diff1(:,dist_n)=emp_x-model_x;
        x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
        
        r=xbound*(1-rand(nx,ntrials)).^(-1/(pareto_fit-1));
        emp_r=ecdfcalc(r);
        
        model_r=zeros(nx,ntrials);
        for j=1:ntrials
            current_r=r(:,j);
            flike_r=@(fit) -sum(log(pareto_pdf(current_r,fit(1))));
            r_fit=fminsearch(flike_r,pareto_fit);
            model_r(1:length(unique(current_r)),j)=pareto_cdf(unique(current_r),r_fit);
        end
        
        ksstat_rand1=max(abs(emp_r-model_r));
        ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
        pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
        pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
        
        like(:,dist_n)=(log(pareto_pdf(x,pareto_fit)));
        aic(dist_n)=2*flike(pareto_fit)+2*length(pareto_fit);
        rsquared(dist_n)=1-sum((emp_x_all-pareto_cdf(x,pareto_fit)).^2)./sum((pareto_cdf(x,pareto_fit)-mean(pareto_cdf(x,pareto_fit))).^2);
        
        disp('Pareto')
        toc(tstart)
        tstart=tic;
        %% 2.generalized pareto
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=gpfit(x-xbound);
        
        flike=@(fit) -sum(log(gppdf(x,fit(1),fit(2),xbound)));
        gp_fit=fminsearch(flike, fit1);
        gp_cdf=@(x,gp_fit) gpcdf(x,gp_fit(1),gp_fit(2),xbound);
        gp_ccdf= @(x,gp_fit) (1-gp_cdf(x,gp_fit))./(1-gp_cdf(xmin,gp_fit));
        model_x=gp_cdf(unique_x,gp_fit);
        ksstat1=max(abs(emp_x-model_x));
        ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
        
        x_diff1(:,dist_n)=emp_x-model_x;
        x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
        
        r=gprnd(gp_fit(1),gp_fit(2),xbound,nx,ntrials);
        emp_r=ecdfcalc(r);
        oldopts=optimset('fminsearch');
        opts=optimset(oldopts,'MaxFunEvals',15000,'MaxIter',13000);
        
        model_r=zeros(nx,ntrials);
        for j=1:ntrials
            current_r=r(:,j);
            flike_r=@(fit) -sum(log(gppdf(current_r,fit(1),fit(2),xbound)));
            r_fit=fminsearch(flike_r,gpfit(current_r-xbound),opts);
            model_r(1:length(unique(current_r)),j)=gp_cdf(unique(current_r),r_fit);
        end
        
        ksstat_rand1=max(abs(emp_r-model_r));
        ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
        pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
        pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
        
        like(:,dist_n)=(log(gppdf(x,gp_fit(1),gp_fit(2),xbound)));
        aic(dist_n)=2*flike(gp_fit)+2*length(gp_fit);
        rsquared(dist_n)=1-sum((emp_x_all-gp_cdf(x,gp_fit)).^2)./sum((gp_cdf(x,gp_fit)-mean(gp_cdf(x,gp_fit))).^2);
        
        disp('GP')
        toc(tstart)
        tstart=tic;
        %% 3.exponential (continuous)
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        exp_trunc = @(x,sigma) exppdf(x,sigma)./(1-expcdf(xbound,sigma));
        fit1=expfit(x-xbound);
        flike=@(fit) -sum(log(exp_trunc(x,fit)));
        exp_fit=fminsearch(flike, fit1);
        exp_cdf = @(x,exp_fit) (expcdf(x,exp_fit)-expcdf(xbound,exp_fit))./(1-expcdf(xbound,exp_fit));
        exp_ccdf= @(x,exp_fit) (1-exp_cdf(x,exp_fit))./(1-exp_cdf(xmin,exp_fit));
        
        model_x=exp_cdf(unique_x,exp_fit);
        ksstat1=max(abs(emp_x-model_x));
        ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
        
        x_diff1(:,dist_n)=emp_x-model_x;
        x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
        
        %Random number generation
        %generates numbers greater than xbound
        r=xbound-(exp_fit)*log(rand(nx,ntrials));
        emp_r=ecdfcalc(r);
        
        model_r=zeros(1,ntrials);
        for j=1:ntrials
            current_r=r(:,j);
            flike_r=@(fit) -sum(log(exp_trunc(current_r,fit)));
            oldopts=optimset('fminsearch');
            opts=optimset(oldopts,'MaxFunEvals',15000,'MaxIter',13000);
            r_fit=fminsearch(flike_r,exp_fit,opts);
            model_r(1:length(unique(current_r)),j)=exp_cdf(unique(current_r),r_fit);
        end
        
        %KS and Modified KS (Clauset et al.)
        ksstat_rand1=max(abs(emp_r-model_r));
        ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
        pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
        pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
        
        like(:,dist_n)=(log(exp_trunc(x,exp_fit)));
        aic(dist_n)=2*flike(exp_fit)+2*length(exp_fit);
        rsquared(dist_n)=1-sum((emp_x_all-exp_cdf(x,exp_fit)).^2)./sum((exp_cdf(x,exp_fit)-mean(exp_cdf(x,exp_fit))).^2);
        
        disp('Exponential')
        toc(tstart)
        tstart=tic;
        %% 4.lognormal
        dist_n=dist_n+1;
        
        
        if ismember(site,[])
            fit1=mle(x-xbound,'distribution','lognormal');
        else
            fit1=mle(x,'distribution','lognormal');
        end
        logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        
        if (data_n==1 && site~=26) || data_n==2
%             if site==1
                logn_fit=mle(x, 'pdf',logn_trunc, 'start',[-13 6], 'lowerbound',[-inf 0],'optimfun','fmincon');
%             else
%                 logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0],'optimfun','fmincon');
%             end
        else
            logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        end
        
        %
        % if ismember(site,[5 6 8 18])
        % if xmin>1000
        % fit1=mle(x-xmin,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % oldopts=statset('mlecustom');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0],'options',opts);
        % else
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        % end
        %
        % elseif ismember(site,[3 10 16 20 21])
        % if xmin>1000
        % fit2=[0 2];
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % oldopts=statset('mlecustom');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit2, 'lowerbound',[-inf 0],'options',opts);
        % else
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        % end
        %
        % elseif ismember(site,[4])
        % if xmin>1000
        % fit2=[-5 3];
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % oldopts=statset('mlecustom');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit2, 'lowerbound',[-inf 0],'options',opts);
        % else
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        % end
        %
        % elseif ismember(site,[7])
        % if xmin>1000
        % fit2=[10 1];
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % oldopts=statset('mlecustom');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit2, 'lowerbound',[-inf 0],'options',opts);
        % else
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        % end
        %
        % else
        % if xmin>1000
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % oldopts=statset('mlecustom');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0], 'options',opts);
        % else
        % fit1=mle(x,'distribution','lognormal');
        % logn_trunc = @(x,muu,sigma) lognpdf(x,muu,sigma)./(1-logncdf(xbound,muu,sigma));
        % logn_fit=mle(x, 'pdf',logn_trunc, 'start',fit1, 'lowerbound',[-inf 0]);
        % end
        % end
        
        
        flike=@(fit) -sum(log(logn_trunc(x,fit(1),fit(2))));
        logn_cdf = @(x,logn_fit) (logncdf(x,logn_fit(1),logn_fit(2))-logncdf(xbound,logn_fit(1),logn_fit(2)))./(1-logncdf(xbound,logn_fit(1),logn_fit(2)));
        logn_ccdf= @(x,logn_fit) (1-logn_cdf(x,logn_fit))./(1-logn_cdf(xmin,logn_fit));
        
        if logncdf(xbound,logn_fit(1),logn_fit(2)) < .9999
            model_x=logn_cdf(unique_x,logn_fit);
            ksstat1=max(abs(emp_x-model_x));
            ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
            
            x_diff1(:,dist_n)=emp_x-model_x;
            x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
            
            
            r1=logncdf(xbound,logn_fit(1),logn_fit(2))+(logncdf((density(2)+density(3)+density(4))*mapelements,logn_fit(1),logn_fit(2))-logncdf(xbound,logn_fit(1),logn_fit(2))).*rand(nx,ntrials);
            r=logninv(r1,logn_fit(1),logn_fit(2));
            
            % r=0;
            % j=1;
            % while sum(r(:)>0)<nx*ntrials
            % r1=lognrnd(logn_fit(1),logn_fit(2),2*nx*ntrials,1);
            % r1=r1(r1>xbound);
            % r(find(r==0,1):find(r==0,1)+length(r1)-1,1)=r1;
            % r=padarray(r,[1 0],'post');
            % end;
            % r=r(1:nx*ntrials);
            % r=reshape(r,nx,ntrials);
            
            emp_r=ecdfcalc(r);
            
            model_r=zeros(1,ntrials);
            for j=1:ntrials
                current_r=r(:,j);
                flike_r=@(fit) -sum(log(logn_trunc(current_r,fit(1),fit(2))));
                if data_n==1
                    if ismember(site,[])
                        oldopts=optimset('fminsearch');
                        opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',10000);
                        r_fit=fminsearch(flike_r,mle(current_r,'distribution','lognormal'),opts);
                    else
                        oldopts=optimset('fminsearch');
                        opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',10000);
                        r_fit=mle(current_r, 'pdf',logn_trunc, 'start',mle(current_r,'distribution','lognormal'), 'lowerbound',[-inf 0], 'optimfun','fmincon','options',opts);
                    end
                    
                else
                    oldopts=optimset('fminsearch');
                    opts=optimset(oldopts,'MaxFunEvals',25000,'MaxIter',25000);
                    r_fit=fminsearch(flike_r,logn_fit,opts);
                end
                
                model_r(1:length(unique(current_r)),j)=logn_cdf(unique(current_r),r_fit);
            end
            
            %KS and Modified KS (Clauset et al.)
            ksstat_rand1=max(abs(emp_r-model_r));
            ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
            pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
            pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
            
        else
            pscore1(dist_n) = -1;
            pscore2(dist_n) = -1;
        end
        
        like(:,dist_n)=(log(logn_trunc(x,logn_fit(1),logn_fit(2))));
        aic(dist_n)=2*flike(logn_fit)+2*length(logn_fit);
        rsquared(dist_n)=1-sum((emp_x_all-logn_cdf(x,logn_fit)).^2)./sum((logn_cdf(x,logn_fit)-mean(logn_cdf(x,logn_fit))).^2);
        
        disp('Log-normal')
        toc(tstart)
        tstart=tic;
        %% 5.inverse guassian
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=mle(x-xbound,'distribution','inversegaussian');
        invgaus_trunc = @(x,muu,sigma) pdf('inversegaussian',x,muu,sigma)./(1-cdf('inversegaussian',xbound,muu,sigma));
        flike=@(fit) -sum(log(invgaus_trunc(x,fit(1),fit(2))));
        oldopts=optimset('fminsearch');
        opts=optimset(oldopts,'MaxFunEvals',15000,'MaxIter',13000);
        invgaus_fit=fminsearch(flike,fit1,opts);
        invgaus_cdf = @(x,invgaus_fit) (cdf('inversegaussian',x,invgaus_fit(1),invgaus_fit(2))-cdf('inversegaussian',xbound,invgaus_fit(1),invgaus_fit(2)))./(1-cdf('inversegaussian',xbound,invgaus_fit(1),invgaus_fit(2)));
        invgaus_ccdf= @(x,invgaus_fit) (1-invgaus_cdf(x,invgaus_fit))./(1-invgaus_cdf(xmin,invgaus_fit));
        
        if cdf('inversegaussian',xbound,invgaus_fit(1),invgaus_fit(2)) < .99
            model_x=invgaus_cdf(unique_x,invgaus_fit);
            ksstat1=max(abs(emp_x-model_x));
            ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
            
            x_diff1(:,dist_n)=emp_x-model_x;
            x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
            
            r=0;
            j=1;
            while sum(r(:)>0)<nx*ntrials
                r1=random(ProbDistUnivParam('inversegaussian',[invgaus_fit(1) invgaus_fit(2)]),[2*nx*ntrials,1]);
                r1=r1(r1>xbound);
                r(find(r==0,1):find(r==0,1)+length(r1)-1,1)=r1;
                r=padarray(r,[1 0],'post');
            end;
            r=r(1:nx*ntrials);
            r=reshape(r,nx,ntrials);
            
            emp_r=ecdfcalc(r);
            
            model_r=zeros(1,ntrials);
            for j=1:ntrials
                current_r=r(:,j);
                flike_r=@(fit) -sum(log(invgaus_trunc(current_r,fit(1),fit(2))));
                r_fit=fminsearch(flike_r,invgaus_fit,opts);
                model_r(1:length(unique(current_r)),j)=invgaus_cdf(unique(current_r),r_fit);
            end
            
            ksstat_rand1=max(abs(emp_r-model_r));
            ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
            pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
            pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
            
        else
            pscore1(dist_n) = -1;
            pscore2(dist_n) = -1;
        end
        
        like(:,dist_n)=(log(invgaus_trunc(x,invgaus_fit(1),invgaus_fit(2))));
        aic(dist_n)=2*flike(invgaus_fit)+2*length(invgaus_fit);
        rsquared(dist_n)=1-sum((emp_x_all-invgaus_cdf(x,invgaus_fit)).^2)./sum((invgaus_cdf(x,invgaus_fit)-mean(invgaus_cdf(x,invgaus_fit))).^2);
        
        disp('Inverse Gaussian')
        toc(tstart)
        tstart=tic;
        %% 6.gamma
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=gamfit(x-xbound);
        gam_trunc = @(x,a,b) pdf('gamma',x,a,b)./(1-cdf('gamma',xbound,a,b));
        flike=@(fit) -sum(log(gam_trunc(x,fit(1),fit(2))));
        oldopts=optimset('fminsearch');
        opts=optimset(oldopts,'MaxFunEvals',15000,'MaxIter',10000);
        gam_fit=fminsearch(flike, fit1, opts);
        gam_cdf = @(x,gam_fit) (cdf('gamma',x,gam_fit(1),gam_fit(2))-cdf('gamma',xbound,gam_fit(1),gam_fit(2)))./(1-cdf('gamma',xbound,gam_fit(1),gam_fit(2)));
        gam_ccdf= @(x,gam_fit) (1-gam_cdf(x,gam_fit))./(1-gam_cdf(xmin,gam_fit));
        
        
        if cdf('gamma',xbound,gam_fit(1),gam_fit(2)) < .99
            
            model_x=gam_cdf(unique_x,gam_fit);
            ksstat1=max(abs(emp_x-model_x));
            ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
            
            x_diff1(:,dist_n)=emp_x-model_x;
            x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
            
            r=0;
            j=1;
            while sum(r(:)>0)<nx*ntrials
                r1=gamrnd(gam_fit(1),gam_fit(2),2*nx*ntrials,1);
                r1=r1(r1>xbound);
                r(find(r==0,1):find(r==0,1)+length(r1)-1,1)=r1;
                r=padarray(r,[1 0],'post');
            end;
            r=r(1:nx*ntrials);
            r=reshape(r,nx,ntrials);
            
            emp_r=ecdfcalc(r);
            
            model_r=zeros(1,ntrials);
            for j=1:ntrials
                current_r=r(:,j);
                flike_r=@(fit) -sum(log(gam_trunc(current_r,fit(1),fit(2))));
                r_fit=fminsearch(flike_r,gam_fit);
                model_r(1:length(unique(current_r)),j)=gam_cdf(unique(current_r),r_fit);
            end
            
            ksstat_rand1=max(abs(emp_r-model_r));
            ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
            pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
            pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
            
        else
            pscore1(dist_n) = -1;
            pscore2(dist_n) = -1;
        end
        
        like(:,dist_n)=(log(gam_trunc(x,gam_fit(1),gam_fit(2))));
        aic(dist_n)=2*flike(gam_fit)+2*length(gam_fit);
        rsquared(dist_n)=1-sum((emp_x_all-gam_cdf(x,gam_fit)).^2)./sum((gam_cdf(x,gam_fit)-mean(gam_cdf(x,gam_fit))).^2);
        
        disp('Gamma')
        toc(tstart)
        tstart=tic;


        %% 7. Power Law with exponential cutoff
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=[pareto_fit-.05 1/(max(x)/2)];
        if fit1(1)<1
            error('Raise the initial fit for powexp and logwbl')
        end
        powexp_trunc = @(x,alpha,lambda) ((lambda^(1-alpha))/gamma_incomplete(xbound*lambda,1 - alpha))*exp(-x.*lambda).*x.^(-alpha);
        oldopts=statset('mlecustom');
        opts=optimset(oldopts,'MaxFunEvals',15000,'MaxIter',15000);
        if data_n==1 || data_n==2
            powexp_fit=mle(x, 'pdf',powexp_trunc, 'start',fit1, 'lowerbound',[1 0],'options',opts);
        else
            powexp_fit=mle(x, 'pdf',powexp_trunc, 'start',fit1, 'lowerbound',[1 0],'options',opts, 'optimfun','fmincon');
        end
        flike=@(fit) -sum(log(powexp_trunc(x,fit(1),fit(2))));
        powexp_cdf =@(x,powexp_fit) (1-(((powexp_fit(2)^(1-powexp_fit(1)))/gamma_incomplete(xbound*powexp_fit(2),1 - powexp_fit(1)))).*(x').^(1 - powexp_fit(1)).*(powexp_fit(2)*(x')).^(powexp_fit(1) - 1).*gamma_incomplete(powexp_fit(2).*x,1 - powexp_fit(1)))';
        powexp_ccdf= @(x,powexp_fit) (1-powexp_cdf(x,powexp_fit))./(1-powexp_cdf(xmin,powexp_fit));
        
        if 1/powexp_fit(2) > xbound && 1/powexp_fit(2) <= max(x)
            model_x=powexp_cdf(unique_x,powexp_fit);
            ksstat1=max(abs(emp_x-model_x));
            ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
            
            x_diff1(:,dist_n)=emp_x-model_x;
            x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
            
            r=0;
            j=1;
            while sum(r(:)>0)<nx*ntrials
                r1=xbound-(1/powexp_fit(2)).*log(rand(2*nx*ntrials,1));
                r1=r1(rand(2*nx*ntrials,1) < (r1./xbound).^-powexp_fit(1));
                r(find(r==0,1):find(r==0,1)+length(r1)-1,1)=r1;
                r=padarray(r,[1 0],'post');
            end;
            r=r(1:nx*ntrials);
            r=reshape(r,nx,ntrials);
            
            emp_r=ecdfcalc(r);
            
            model_r=zeros(1,ntrials);
            for j=1:ntrials
                current_r=r(:,j);
                flike_r=@(fit) -sum(log(powexp_trunc(current_r,fit(1),fit(2))));
                r_fit=fminsearch(flike_r,powexp_fit);
                model_r(1:length(unique(current_r)),j)=powexp_cdf(unique(current_r),r_fit);
            end
            
            ksstat_rand1=max(abs(emp_r-model_r));
            ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
            pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
            pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
            
        else
            pscore1(dist_n) = -1;
            pscore2(dist_n) = -1;
        end
        
        like(:,dist_n)=(log(powexp_trunc(x,powexp_fit(1),powexp_fit(2))));
        aic(dist_n)=2*flike(powexp_fit)+2*length(powexp_fit);
        rsquared(dist_n)=1-sum((emp_x_all-powexp_cdf(x,powexp_fit)).^2)./sum((powexp_cdf(x,powexp_fit)-mean(powexp_cdf(x,powexp_fit))).^2);
        
        
        disp('Pareto with exponential cutoff')
        toc(tstart)
        tstart=tic;
        %% 8. Weibull (Stretched Exponential)
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=wblfit(x);
        wbl_trunc = @(x,a,b) wblpdf(x,a,b)./(1-wblcdf(xbound,a,b));
        flike=@(fit) -sum(log(wbl_trunc(x,fit(1),fit(2))));
        % oldopts=optimset('fminsearch');
        % opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',5000);
        oldopts=statset('mlecustom');
        opts=optimset(oldopts,'MaxFunEvals',10000,'MaxIter',8000);
        wbl_fit=mle(x,'pdf',wbl_trunc,'start',fit1,'lowerbound',[0 0], 'options',opts);
        
        % wbl_fit=fminsearch(flike,fit1,opts);
        wbl_cdf = @(x,wbl_fit) (wblcdf(x,wbl_fit(1),wbl_fit(2))-wblcdf(xbound,wbl_fit(1),wbl_fit(2)))/(1-wblcdf(xbound,wbl_fit(1),wbl_fit(2)));
        wbl_ccdf= @(x,wbl_fit) (1-wbl_cdf(x,wbl_fit))./(1-wbl_cdf(xmin,wbl_fit));
%         
%         if wblcdf(xbound,wbl_fit(1),wbl_fit(2)) < .99
%             model_x=wbl_cdf(unique_x,wbl_fit);
%             ksstat1=max(abs(emp_x-model_x));
%             ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
%             
%             x_diff1(:,dist_n)=emp_x-model_x;
%             x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
%             
%             r=(xbound^wbl_fit(2)-(1/(wbl_fit(1)^-wbl_fit(2)))*log(1-rand(nx,ntrials))).^(1/wbl_fit(2));
%             emp_r=ecdfcalc(r);
%             
%             model_r=zeros(1,ntrials);
%             for j=1:ntrials
%                 current_r=r(:,j);
%                 flike_r=@(fit) -sum(log(wbl_trunc(current_r,fit(1),fit(2))));
%                 if data_n==1 || data_n==2
%                     r_fit=mle(current_r,'pdf',wbl_trunc,'start',wbl_fit,'lowerbound',[0 0],'options',opts,'optimfun','fmincon');
%                 else
%                     r_fit=fminsearch(flike_r,wbl_fit,opts);
%                 end
%                 model_r(1:length(unique(current_r)),j)=wbl_cdf(unique(current_r),r_fit);
%             end
%             
%             ksstat_rand1=max(abs(emp_r-model_r));
%             ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
%             pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
%             pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
%             
%         else
            pscore1(dist_n) = -1;
            pscore2(dist_n) = -1;
%         end
%         
        like(:,dist_n)=(log(wbl_trunc(x,wbl_fit(1),wbl_fit(2))));
        aic(dist_n)=2*flike(wbl_fit)+2*length(wbl_fit);
        rsquared(dist_n)=1-sum((emp_x_all-wbl_cdf(x,wbl_fit)).^2)./sum((wbl_cdf(x,wbl_fit)-mean(wbl_cdf(x,wbl_fit))).^2);
        
        disp('Weibull')
        toc(tstart)
        tstart=tic;
        %% 9. Log-Wiebull (gumbel)
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=[pareto_fit-.5 1];
        logw_trunc = @(x,b,c) (b*c./x).*((log(x./xbound)).^(c-1)).*exp(-b.*(log(x./xbound).^c));
        flike=@(fit) -sum(log(logw_trunc(x,fit(1),fit(2))));
        logw_fit=mle(x, 'pdf',logw_trunc, 'start',fit1, 'lowerbound',[0 0]);
        logw_cdf = @(x,logw_fit) 1-exp(-logw_fit(1)*(log(x./xbound)).^logw_fit(2));
        logw_ccdf= @(x,logw_fit) (1-logw_cdf(x,logw_fit))./(1-logw_cdf(xmin,logw_fit));
        
        like(:,dist_n)=(log(logw_trunc(x,logw_fit(1),logw_fit(2))));
        aic(dist_n)=2*flike(logw_fit)+2*length(logw_fit);
        rsquared(dist_n)=1-sum((emp_x_all-logw_cdf(x,logw_fit)).^2)./sum((logw_cdf(x,logw_fit)-mean(logw_cdf(x,logw_fit))).^2);
        
        disp('Log-weibull')
        toc(tstart)
        %% 10 Power Law with Hyperexponential Cutoff
        dist_n=dist_n+1;
        
 
        
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        fit1=[powexp_fit(1) powexp_fit(2)];
        hypexp_trunc = @(x,alpha,lambda) 2*lambda^(1-alpha)./(x.^alpha.*exp(lambda^2.*x.^2).*gamma_incomplete(lambda^2*xbound^2,1/2-alpha/2));
        oldopts=statset('mlecustom');
        opts=optimset(oldopts,'MaxFunEvals',25000,'MaxIter',25000);
%         if ismember(data_n,[1 3]) || (data_n==8) % && ismember(site,[2 3]))
            hypexp_fit=mle(x, 'pdf',hypexp_trunc, 'start',fit1, 'lowerbound',[1 0],'options',opts);
%         else
%             hypexp_fit=mle(x, 'pdf',hypexp_trunc, 'start',fit1, 'lowerbound',[1 0],'options',opts, 'optimfun','fmincon');
%         end
        flike=@(fit) -sum(log(hypexp_trunc(x,fit(1),fit(2))));
        hypexp_cdf = @(x,hypexp_fit) 1-gamma_incomplete(x.^2*hypexp_fit(2)^2,1/2-hypexp_fit(1)/2)'./gamma_incomplete(xbound^2.*hypexp_fit(2)^2,1/2-hypexp_fit(1)/2);
        hypexp_ccdf= @(x,hypexp_fit) (1-hypexp_cdf(x,hypexp_fit))./(1-hypexp_cdf(xmin,hypexp_fit));
        
        pscore1(dist_n) = -1;
        pscore2(dist_n) = -1;
        
        like(:,dist_n)=(log(hypexp_trunc(x,hypexp_fit(1),hypexp_fit(2))));
        aic(dist_n)=2*flike(hypexp_fit)+2*length(hypexp_fit);
        rsquared(dist_n)=1-sum((emp_x_all-hypexp_cdf(x,hypexp_fit)).^2)./sum((hypexp_cdf(x,hypexp_fit)-mean(hypexp_cdf(x,hypexp_fit))).^2);
        
        % if ismember(data_n,3)
        % figure(site)
        % loglog(unique(x),ccdfcalc(x),'.k');
        % hold on
        % fplot(@(x) powexp_ccdf(x,powexp_fit),[min(x) max(x)],'g');
        % fplot(@(x) hypexp_ccdf(x,hypexp_fit),[min(x) max(x)],'r');
        % hold off
        % drawnow
        % end
        
        disp('Pareto with hyper-exponential cutoff')
        toc(tstart)
        tstart=tic;
        
        %% 11.exponential (discrete)
        
        dist_n=dist_n+1;
        clear fit1 flike model_x ksstat1 ksstat2 r emp_r model_r current_r flike_r r_fit model_r ksstat_rand1 ksstat_rand2 r1
        
        discexp=@(x,fit) (1-exp(-fit))*exp(fit*xmin)*exp(-fit*x);
        discfit1=log(1+nx/sum(x-xmin));
        flike=@(fit) -sum(log(discexp(x,fit)));
        discexp_fit=fminsearch(flike, discfit1);
        
        if data_n==4 || data_n==5
            
            discexpcdf1=0;
            discexp_cdf_all=zeros(max(x)-min(x)+1,1);
            for i=xmin:max(x)*10
                discexp_cdf_all(i)=discexp(i,discexp_fit)+discexpcdf1;
                discexpcdf1=discexp_cdf_all(i);
            end;
            discexp_all=discexp_cdf_all;
            model_x=discexp_cdf_all(unique_x);
            
            discexp_ccdf= (1-model_x./(1-model_x(2)));
            
            ksstat1=max(abs(emp_x-model_x));
            ksstat2=max(abs(emp_x-model_x)./sqrt(model_x.*(1-model_x)));
            
            x_diff1(:,dist_n)=emp_x-model_x;
            x_diff2(:,dist_n)=(emp_x-model_x)./sqrt(model_x.*(1-model_x));
            
            r=zeros(nx,ntrials);
            for j=1:ntrials
                for n=1:nx
                    r(n,j)=find(discexp_all>rand,1);
                end
            end
            
            if max(r)>=max(x)*9
                error('increase max for discexp_all');
            end
            
            emp_r=ecdfcalc(r);
            model_r=zeros(1,ntrials);
            for j=1:ntrials
                current_r=r(:,j);
                flike_r=@(fit) -sum(log(discexp(current_r,fit)));
                r_fit=fminsearch(flike_r,discexp_fit);
                model_r(1:length(unique(current_r)),j)=discexp_cdf_all(unique(current_r));
            end
            
            ksstat_rand1=max(abs(emp_r-model_r));
            ksstat_rand2=max(abs(emp_r-model_r)./sqrt(model_r.*(1-model_r)));
            pscore1(dist_n)=sum(ksstat_rand1>ksstat1)./ntrials;
            pscore2(dist_n)=sum(ksstat_rand2>ksstat2)./ntrials;
            
            rsquared(dist_n)=1-sum((emp_x_all-discexp_cdf_all(x)).^2)./sum((discexp_cdf_all(x)-mean(discexp_cdf_all(x))).^2);
        end
        
        disp('Discrete exponential')
        toc(tstart)
        tstart=tic;
        
        %% Calculate log likelihoods
        
        like_ratio=zeros(10,10);
        like_variance=zeros(10,10);
        like_p=zeros(10,10);
        delta=zeros(10,10);
        weight=zeros(10,10);
        
        for i=1:10
            for j=1:10
                like_ratio(i,j)=sum(like(:,i)-like(:,j));
                if like_ratio(i,j)<=0
                    like_ratio(i,j)=NaN;
                end
                like_variance(i,j)=(1/nx)*sum(((like(:,i)-like(:,j))-(mean(like(:,i))-mean(like(:,j)))).^2);
                like_p(i,j)=1-abs(erf(like_ratio(i,j)/(sqrt(like_variance(i,j))*sqrt(2*nx))));
                delta(i,j)=aic(i)-min([aic(i) aic(j)]);
            end
        end
        for i=1:10
            for j=1:10
                weight(i,j)=exp(-delta(i,j)/2)/(exp(-delta(i,j)/2)+exp(-delta(j,i)/2));
            end
        end
        
        %% Write to excel
        
        loglikeresults=-sum(like);
        
        if sum(aic==-inf) > 0
            disp('fitting didnt work for distribution:');
            baddist=find(aic==-inf)
            loglikeresults(baddist)=inf;
            aic(baddist)=inf;
        end
        
        [m best_like_index] = min(loglikeresults);
        [m best_aic_index] = min(aic);
        best_p = like_p(best_like_index,:);
        best_p(best_like_index)=1;
        best_p(best_p<.0001)=0;
        best_aic = weight(:,best_aic_index)';
        best_aic(best_aic_index)=1;
        best_aic(best_aic<.0001)=0;
        
        header={'Areas' 'Pareto' 'GP' 'Exp' 'Logn' 'Inv Gauss' 'Gamma' 'Pow-Exp' 'Weibull' 'Log-Wbl' 'Pow-Hyp' ' ' 'Min x' 'Max x' 'Mean x' 'Med x' 'hypexp Cutoff' 'max>&min<' 'GP K' 'Wbl B' 'Heavy Tail' 'Logwbl C' 'Gamma A' ' ' ' ' 'GP(k)' 'GP(a)' 'Pareto' 'Exp' 'Logn(mu)' 'Logn(s)' 'igaus(mu)' 'igaus(s)' 'gam(a)' 'gam(b)' 'pexp(a)' 'pexp(l)' 'wbl(a)' 'wbl(b)' 'lwbl(b)' 'lwbl(c)' 'dexp(l)' 'hypexp(a)' 'hypexp(l)'};
        xlswrite(excelfilename,header,'likelihood','A1');
        excel_index = ['B' int2str((2+site)+(data_n-1)*25)];
        xlswrite(excelfilename,best_p,'likelihood',excel_index);
        header={'Areas' 'Pareto' 'GP' 'Exp' 'Logn' 'Inv Gauss' 'Gamma' 'Pow-Exp' 'Weibull' 'Log-Wbl'};
        xlswrite(excelfilename,header,'aic','A1');
        xlswrite(excelfilename,best_aic,'aic',excel_index);
        xlswrite(excelfilename,header,'p_score1','A1');
        xlswrite(excelfilename,pscore1,'p_score1',excel_index);
        xlswrite(excelfilename,header,'p_score2','A1');
        xlswrite(excelfilename,pscore2,'p_score2',excel_index);
        header={'Areas' 'Pareto' 'GP' 'Exp' 'Logn' 'Inv Gauss' 'Gamma' 'Pow-Exp' 'Weibull' 'Log-Wbl' 'DiscExp'};
        xlswrite(excelfilename,header,'R-squared','A1');
        xlswrite(excelfilename,rsquared,'R-squared',excel_index);
        
        distparams=[min(x) max(x) mean(x) median(x) 1/hypexp_fit(2) min(x)<1/hypexp_fit(2)&max(x)>1/hypexp_fit(2) gp_fit(1) wbl_fit(2) gp_fit(1)>0&wbl_fit(2)<1 logw_fit(2) gam_fit(1)];
        excel_index = ['M' int2str((2+site)+(data_n-1)*25)];
        xlswrite(excelfilename,distparams,'likelihood',excel_index);
        
        distparams=[gp_fit pareto_fit exp_fit logn_fit invgaus_fit gam_fit powexp_fit wbl_fit logw_fit 1/discexp_fit hypexp_fit];
        excel_index = ['Z' int2str((2+site)+(data_n-1)*25)];
        xlswrite(excelfilename,distparams,'likelihood',excel_index);
        
        if data_n==1
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'Areas_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'Areas_CDF',excel_index);
            excel_index = [idx2A1(2*site+1) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'Areas_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'Areas_CDF',excel_index);
            excel_index = [idx2A1(2*site+1+52) '2'];
            xlswrite(excelfilename,pareto_ccdf(unique_x,pareto_fit),'Areas_CDF',excel_index);
        
            
        end
        
        if data_n==4
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WSloughs_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'WSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+1) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'WSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+51) '2'];
            xlswrite(excelfilename,unique_x,'WSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+1+51) '2'];
            xlswrite(excelfilename,exp_ccdf(unique_x,exp_fit),'WSloughs_CDF',excel_index);
            
            header={'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Discexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Logn:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Gamma:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Weibull:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WSloughs_1','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,3),'WSloughs_1',excel_index);
            excel_index = ['A' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,10),'WSloughs_1',excel_index);
            excel_index = ['B' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,4),'WSloughs_1',excel_index);
            excel_index = ['C' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,6),'WSloughs_1',excel_index);
            excel_index = ['D' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,8),'WSloughs_1',excel_index);
            
            header={'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Discexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Logn:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Gamma:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Weibull:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WSloughs_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,3),'WSloughs_2',excel_index);
            excel_index = ['A' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,10),'WSloughs_2',excel_index);
            excel_index = ['B' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,4),'WSloughs_2',excel_index);
            excel_index = ['C' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,6),'WSloughs_2',excel_index);
            excel_index = ['D' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,8),'WSloughs_2',excel_index);
        end
        
        if data_n==5
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LSloughs_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'LSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+1) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'LSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'LSloughs_CDF',excel_index);
            excel_index = [idx2A1(2*site+52+1) '2'];
            xlswrite(excelfilename,exp_ccdf(unique_x,exp_fit),'LSloughs_CDF',excel_index);
            
            header={'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Discexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Logn:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Gamma:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Weibull:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LSloughs_1','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,3),'LSloughs_1',excel_index);
            excel_index = ['A' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,10),'LSloughs_1',excel_index);
            excel_index = ['B' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,4),'LSloughs_1',excel_index);
            excel_index = ['C' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,6),'LSloughs_1',excel_index);
            excel_index = ['D' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,8),'LSloughs_1',excel_index);
            
            header={'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Discexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Logn:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Gamma:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Weibull:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LSloughs_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,3),'LSloughs_2',excel_index);
            excel_index = ['A' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,10),'LSloughs_2',excel_index);
            excel_index = ['B' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,4),'LSloughs_2',excel_index);
            excel_index = ['C' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,6),'LSloughs_2',excel_index);
            excel_index = ['D' idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff2(:,8),'LSloughs_2',excel_index);
        end
        
        if data_n==2
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WRidge_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'WRidge_CDF',excel_index);
            excel_index = [idx2A1(2*site+1) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'WRidge_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'WRidge_CDF',excel_index);
            excel_index = [idx2A1(1+2*site+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'WRidge_CDF',excel_index);
            
            header={'PowExp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WRidge_1','A1');
            xlswrite(excelfilename,header,'WRidge_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'WRidge_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'Wridge_2',excel_index);
        end
        
        if data_n==3
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Powexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Exp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LRidge_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'LRidge_CDF',excel_index);
            excel_index = [idx2A1(1+site*2) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'LRidge_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'LRidge_CDF',excel_index);
            excel_index = [idx2A1(1+site*2+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'LRidge_CDF',excel_index);
            excel_index = [idx2A1(1+site*2+102) '2'];
            xlswrite(excelfilename,exp_ccdf(unique_x,exp_fit),'LRidge_CDF',excel_index);
            
            header={'Pexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LRidge_1','A1');
            xlswrite(excelfilename,header,'LRidge_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'LRidge_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'Lridge_2',excel_index);
        end
        
        if data_n==8
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WRidgeOK_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'WRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(1+site*2) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'WRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'WRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(1+site*2+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'WRidgeOK_CDF',excel_index);
            
            header={'Pexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WRidgeOK_1','A1');
            xlswrite(excelfilename,header,'WRidgeOK_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'WRidgeOK_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'WridgeOK_2',excel_index);
        end
        
        if data_n==9
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LRidgeOK_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'LRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(1+2*site) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'LRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'LRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(1+2*site+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'LRidgeOK_CDF',excel_index);
            excel_index = [idx2A1(2*site+1+102) '2'];
            xlswrite(excelfilename,exp_ccdf(unique_x,exp_fit),'LRidgeOK_CDF',excel_index);
            
            header={'Pexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LRidgeOK_1','A1');
            xlswrite(excelfilename,header,'LRidgeOK_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'LRidgeOK_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'LridgeOK_2',excel_index);
        end
        
        if data_n==10
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WNN_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'WNN_CDF',excel_index);
            excel_index = [idx2A1(1+site*2) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'WNN_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'WNN_CDF',excel_index);
            excel_index = [idx2A1(1+site*2+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'WNN_CDF',excel_index);
            
            header={'Pexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'WNN_1','A1');
            xlswrite(excelfilename,header,'WNN_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'WNN_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'WNN_2',excel_index);
        end
        
        if data_n==11
            header={'X:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'ECDF:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 'Pareto:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LNN_CDF','A1');
            excel_index = [idx2A1(2*site) '2'];
            xlswrite(excelfilename,unique_x,'LNN_CDF',excel_index);
            excel_index = [idx2A1(1+site*2) '2'];
            xlswrite(excelfilename,ccdfcalc(x),'LNN_CDF',excel_index);
            excel_index = [idx2A1(2*site+52) '2'];
            xlswrite(excelfilename,unique_x,'LNN_CDF',excel_index);
            excel_index = [idx2A1(1+site*2+52) '2'];
            xlswrite(excelfilename,powexp_ccdf(unique_x,powexp_fit),'LNN_CDF',excel_index);
            
            header={'Pexp:' 'C1' 'C1' 'C2' 'C2' 'C3' 'C3' 'C4' 'C4' 'C5' 'C5' 'C6' 'C6' 'C7' 'C7' 'C8' 'C8' 'C9' 'C9' 'C10' 'C10' 'C11' 'C11' 'C12' 'C12' 'C13' 'C13' 'C14' 'C14' 'C15' 'C15' 'South1' 'South1' 'South2' 'South2' 'South3' 'South3' 'East' 'East' 'North' 'North' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
            xlswrite(excelfilename,header,'LNN_1','A1');
            xlswrite(excelfilename,header,'LNN_2','A1');
            excel_index = [idx2A1(1+site) '2'];
            xlswrite(excelfilename,x_diff1(:,7),'LNN_1',excel_index);
            xlswrite(excelfilename,x_diff2(:,7),'LNN_2',excel_index);
        end
        
        
        header={' ' 'Slough' 'Ridge' 'TI Aff' 'TI' 'Other'};
        xlswrite(excelfilename,header,'Density','A1')
        header={' ' 'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C7' 'C8' 'C9' 'C10' 'C11' 'C12' 'C13' 'C14' 'C15' 'South1' 'South2' 'South3' 'East' 'North'}';
        xlswrite(excelfilename,header,'Density','A1')
        excel_index = ['B' int2str(site+1)];
        xlswrite(excelfilename,density,'Density',excel_index);
        toc(tstart2);
    end;
    
    %% Plot distributions
    
    % Fat Tailed
    
    % figure(1)
    % loglog(unique_x,1-emp_x,'.k');
    % hold on
    % fplot(@(x) 1-powexp_cdf(x,powexp_fit),[min(x) max(x)],'g');
    % hold off
    
    % figure(2)
    % semilogx(unique(r),ecdfcalc(r),'.k');
    % hold on
    % fplot(@(r) exp_cdf(r,exp_fit,xmin),[min(r) max(r)],'g');
    % hold off
    %
    % %Empirical Data KS plots
    % figure(3)
    % plot((emp_r-model_r)./sqrt(model_r.*(1-model_r)))
    % figure(4)
    % plot((emp_r-model_r))
    %
    % %Random Model KS plots
    % figure(6)
    % plot((emp_rand-model_rand)./sqrt(model_rand.*(1-model_rand)));
    % figure(7)
    % plot((emp_rand-model_rand));
    
end
