% Copyright
%%%%%%%%%%%
% This main code was written by O.Damette and S.Goutte
% Please cite this paper if you use this code:

% Olivier Damette & Stephane Goutte & Qing Pei, 2020. 
% "Climate and nomadic migration in a nonlinear world: 
% evidence of the historical China," Climatic Change, Springer, vol.163(4),pages 2055-207.
%
%
% Needs the Patton Copula Toolbox which can be downloaded
% from http://public.econ.duke.edu/~ap172/Patton_copula_toolbox.zip
% This toolbox is a collection of Matlab functions on copulas for financial
% time series. The main papers from that research are listed below.
% http://fmg.lse.ac.uk/~patton
% References
% - Granger, C.W.J, T. Ter‰svirta, and A.J. Patton, 2006 Common Factors in Conditional
% Distributions for Bivariate Time Series, Journal of Econometrics, 132(1), 43-57.   
% - Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence
% for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168.
% - Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic
% Review, 47(2), 527-556.
% - Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different
% Lengths, Journal of Applied Econometrics, 21(2), 147-173.  


%Cleaning of data.
clc;
clear all;
close all

%Starting Time clock
tic

%Reading of an example data
Data_Example;
Label={'Year','RSL I7G','RSL SD I7G','SST','SD SST','Surface T°C Average annual T°C anamoly','Surface T°C SD','NAO annual value','Sunspots 10 year mean','NAO Index (Hurrell Station Based Annual)'};


U=empiricalCDF(Data); % Empirical transformation

%Choix des variables à croiser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=2:3
    for kkk=4:5
        Var1=kk;%Loop to test all bivariate copulas
        Var2=kkk;% 4 à 10

        %Extraction of variables
        u=U(:,Var1);
        v=U(:,Var2);

        % Extraction of labels
        seriesname={Label{Var1} Label{Var2}};

        % Parameters of the algorithm (optimization)
        transmarg=[u,v];
        options = optimset('Algorithm','interior-point','Display','iter','Hessian','bfgs','MaxFunEvals',1000);
        options = optimset(options,'FinDiffType','central','MaxIter',1500,'TolCon',10^-12,'TolFun',10^-5,'TolX',10^-5);
        dados=transmarg;

        %Panel of copulas
        %%%%%%%%%%%%%%%%%

        % 1. Normal Copula
        kappabarN = corrcoef12(norminv(u),norminv(v));
        LL1 = NormalCopula_CL(kappabarN,[u,v]);

        % 3. Clayton copula
        lower = 0.0001;
        theta0 = 1; %Initial value
        [kappabarCL LL3] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,dados);

        % 4. Gumbel copula
        lower = 1.1;
        theta0 = 1.11;
        [kappabarGL LL4] = fmincon('gumbelCL',theta0,[],[],[],[],lower,[],[],options,dados);

        % 5. Symmetrised Joe-Clayton LLcopula
        lower = [0 , 0 ];
        upper = [ 1 , 1];
        theta0 = [0.25;0.25];
        [kappabarSJC LL5] = fmincon('sym_jc_CL',theta0,[],[],[],[],lower,upper,[],options,dados);

        % 8. Student's t copula
        lower = [-0.9 , 2.1 ];
        upper = [ 0.9 , 100 ];
        theta0 = [kappabarN;10];
        [kappabarT LL8] = fmincon('tcopulaCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);

        kappa0=[1 -1 0]; %initial values for time-varying copulas

        lower = -25*ones(3,1);
        upper =  25*ones(3,1);

        % 11. Time-varying normal Copula
        %kappa0 = [log((1+kappabarN)/(1-kappabarN));0;0];
        [kappa11 LL11] = fmincon('bivnorm_tvp1_CL',kappa0,[],[],[],[],lower,upper,[],options,[u,v],kappabarN);
        [LL10, kappaN] = bivnorm_tvp1_CL(kappa11,[u,v],kappabarN);

        % 13. Time-varying Clayton copula
        [kappa13 LL13] = fmincon('Clayton_tvp1_CL',kappa0,[],[],[],[],lower,upper,[],options,[u,v],kappabarCL);
        [LL13 kappaCL] = Clayton_tvp1_CL(kappa13,[u,v],kappabarCL);

        % % 14. Time-varying Gumbel copula
        [kappa14 LL14] = fmincon('Gumbel_tvp1_CL',kappa0,[],[],[],[],lower,upper,[],options,[u,v],kappabarGL);
        [LL14 kappaGL] = Gumbel_tvp1_CL(kappa14,[u,v],kappabarGL);

        kappa0=[1 -1 0 1 0 -1]; %initial values for time-varying copulas with 2 parameters

        lower = -40*ones(6,1);
        upper =  40*ones(6,1);

        % 15. Time-varying Symetrized Joe-Clayton copula
        [kappa15 LL15] = fmincon('sym_jc_tvp_CL',kappa0,[],[],[],[],lower,upper,[],options,[u,v],kappabarSJC);
        [LL15 tauU tauL] = sym_jc_tvp_CL(kappa15,[u,v],kappabarSJC);

        % Robust standard errors calculation
        numgrad=MyFuncScores('sym_jc_tvp_CL',kappa15',U,kappabarSJC);
        H = hessian_2sided('sym_jc_tvp_CL',kappa15',U,kappabarSJC);
        T = length(U);
        VCV = (inv(H/T)*cov(numgrad)*inv(H/T))/T; % covariance matrix
        robstd1 = sqrt(diag(VCV)); % robust standard errors
        %
        % 18. Time-varying T-copula
        [kappa18 LL18] = fmincon('bivt_tvp1_CL',kappa0,[],[],[],[],lower,upper,[],options,[u,v],kappabarT);
        [LL18, kappaT, nuT] = bivt_tvp1_CL(kappa18,[u,v],kappabarT);

        toc/60/60;
        formatSpec = 'Elapsed time hours %4.2f hours.\n';
        fprintf(formatSpec,ans)

        %Backup of the results 
        save(File(kk));

        %Graphical Output
        %%%%%%%%%%%%%%%%%
        n=size(kappaT,1);
        figure()
        hold on

        subplot(4,1,1)
        hold on
        grid
        plot(kappaN)
        legend('Time-varying Normal Copula')
        title(['Copula: ' seriesname{1} '-' seriesname{2} ] )
        axis([0,n,min(kappaN)-0.1,1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})
        subplot(4,1,2)
        hold on
        grid
        plot(kappaCL)
        legend('Time-varying Clayton copula')
        axis([0,n,min(kappaCL)-0.1,max(kappaCL)+0.1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})
        subplot(4,1,3)
        hold on
        grid
        plot(kappaT)
        legend('Time-varying T-copula')
        axis([0,n,min(kappaT)-0.1,1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})
        subplot(4,1,4)
        hold on
        grid
        plot(kappaGL)
        legend('Time-varying Gumbel copula')
        axis([0,n,min(kappaGL)-0.1,max(kappaGL)+0.1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})


        figure()
        hold on
        grid
        plot(kappaGL)
        %plot(1:length(u),kappabarGL*ones(1,length(u)),'r')
        legend('Time-varying Extreme Positive Dependence')
        %legend('Time-varying Extreme Positive Dependence','Static Extreme Positive Dependence')
        axis([0,n,min(kappaGL)-0.1,max(kappaGL)+0.1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})



        figure()
        hold on
        grid
        plot(kappaCL)
        %plot(1:length(u),kappabarGL*ones(1,length(u)),'r')
        legend('Time-varying Extreme Negative Dependence')
        %legend('Time-varying Extreme Negative Dependence','Static Extreme Positive Dependence')
        axis([0,n,min(kappaCL)-0.1,max(kappaCL)+0.1]);
        xticks([1 5 15 25 35 45 55 64])
        xticklabels({'1956','1960','1970','1980','1990','2000','2010','2019'})



    end
end

