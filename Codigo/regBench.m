% Regressão: Benchmarks

%% Dataset

load dataset;

% data='concrete';
% dataset=R_concrete;

% data='housing';
% dataset=R_housing;

% data='nasa';
% dataset=R_nasa;

% data='rmsd';
% dataset=R_rmsd;

data='yacht';
dataset=R_yacht;

%% Dataset: organização

porcaoTR=0.7000;
[L,C]=size(dataset);
linhasTR=round(porcaoTR*L); % "linhasTR=L", se é CV sem teste
GP.xtrain=dataset(1:linhasTR,1:(C-1));
GP.ytrain=dataset(1:linhasTR,C);
if porcaoTR~=1 % Há teste!
    GP.xtest=dataset((linhasTR+1):L,1:(C-1));
    GP.ytest=dataset((linhasTR+1):L,C);
end

clearvars -except data porcaoTR nCV L C linhasTR GP

%% K-fold CV

m=1; % Não mexer
nCV=1;
k=10;

for j=1:nCV % #experimentos
    for i=1:k % k-fold CV
        if i==1

            % Definição do teste (validação)
            inicio=1;
            fim=round((i/k)*linhasTR);

            gp.userdata.xtest=GP.xtrain(inicio:fim,:);
            gp.userdata.xtrain=GP.xtrain((fim+1):linhasTR,:);

            gp.userdata.ytest=GP.ytrain(inicio:fim,:);
            gp.userdata.ytrain=GP.ytrain((fim+1):linhasTR,:);

        else
            inicio=fim+1;
            fim=round((i/k)*linhasTR);

            gp.userdata.xtest=GP.xtrain(inicio:fim,:);
            gp.userdata.xtrain=[GP.xtrain(1:(inicio-1),:) ; GP.xtrain((fim+1):linhasTR,:)];

            gp.userdata.ytest=GP.ytrain(inicio:fim,:);
            gp.userdata.ytrain=[GP.ytrain(1:(inicio-1),:) ; GP.ytrain((fim+1):linhasTR,:)];

        end
        
        %% Treino
        alfa=0.05;
        xtreino=gp.userdata.xtrain;
        ytreino=gp.userdata.ytrain;

        % gene_outputs
        Xtreino=x2fx(xtreino,'quadratic');
%         Xtreino=x2fx(Xtreino,'quadratic');
        % Estimação 1
        mdl=fitlm(Xtreino,ytreino);

        % Regressores Significantes
        pvalue=mdl.Coefficients.pValue;
        pvalue=pvalue(2:end,:);
        reg_signif=pvalue<alfa;

        % gene_outputs signif e Estimação 2
        Xtreino=Xtreino(:,reg_signif);
        mdl=fitlm(Xtreino,ytreino);
        theta=mdl.Coefficients.Estimate;
        
        GP.results(i,1)=mdl.Rsquared.Adjusted;
        
        %% Validação
        xtest=gp.userdata.xtest;
        ytest=gp.userdata.ytest;
        
        % gene_outputs
        Xtest=x2fx(xtest,'quadratic');
%         Xtest=x2fx(Xtest,'quadratic');
        Xtest=Xtest(:,reg_signif);
        Xtest=[ones(length(ytest),1) Xtest];

        ypredteste=Xtest*theta;

        y=ytest;
        [num_data_points,num_reg]=size(Xtest);
        num_reg=num_reg-1;
        GP.results(i,4)=num_reg;

        STE=sum((ypredteste-mean(y)).^2);
        STQ=sum((y-mean(y)).^2);
        r2=STE/STQ;
        r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));
        
        GP.results(i,2)=r2ajustado;

        %% Teste
        xtest=GP.xtest;
        ytest=GP.ytest;
        
        % gene_outputs
        Xtest=x2fx(xtest,'quadratic');
%         Xtest=x2fx(Xtest,'quadratic');
        Xtest=Xtest(:,reg_signif);
        Xtest=[ones(length(ytest),1) Xtest];

        ypredteste=Xtest*theta;

        y=ytest;
        [num_data_points,num_reg]=size(Xtest);
        num_reg=num_reg-1;

        STE=sum((ypredteste-mean(y)).^2);
        STQ=sum((y-mean(y)).^2);
        r2=STE/STQ;
        r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));
        
        GP.results(i,3)=r2ajustado;
    
    end %  nCV #experimentos
end

clearvars -except GP data
filename=[data 'BENCH' '.mat' ];
save(filename);




%% RL: Treino

% alfa=0.05;
% 
% % gene_outputs
% Xtreino=x2fx(xtreino,'quadratic');
% % Estimação 1
% mdl=fitlm(Xtreino,ytreino);
% 
% % Regressores Significantes
% pvalue=mdl.Coefficients.pValue;
% pvalue=pvalue(2:end,:);
% reg_signif=pvalue<alfa;
% 
% % gene_outputs signif e Estimação 2
% Xtreino=Xtreino(:,reg_signif);
% mdl=fitlm(Xtreino,ytreino);
% theta=mdl.Coefficients.Estimate;

%% RL: Teste

% % gene_outputs
% Xtest=x2fx(xtest,'quadratic');
% Xtest=Xtest(:,reg_signif);
% Xtest=[ones(length(ytest),1) Xtest];
% 
% ypredteste=Xtest*theta;
% 
% y=ytest;
% [num_data_points,num_reg]=size(Xtest);
% num_reg=num_reg-1;
% 
% STE=sum((ypredteste-mean(y)).^2);
% STQ=sum((y-mean(y)).^2);
% r2=STE/STQ;
% r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

