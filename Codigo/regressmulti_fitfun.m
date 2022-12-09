function [fitness,gp,ypredtrain]=regressmulti_fitfun(evalstr,gp)
%% Treino ou Teste

if strcmp(gp.runcontrol.dataset,'treino')==1
    y=gp.userdata.ytrain;
elseif strcmp(gp.runcontrol.dataset,'teste')==1
    y=gp.userdata.ytest;
end

%% Preparação

alfa=0.05; % Teste bicaudal: usar alfa/2 = 0.025
num_data_points=size(y,1); % recebe #linhas de y
num_genes=length(evalstr);

%% Design Matrix

gene_outputs=ones(num_data_points,num_genes+1);
for i=1:num_genes
    ind=i+1;
    eval(['gene_outputs(:,ind)=' evalstr{i} ';'])
    if  any(isnan(gene_outputs(:,ind))) || any(isinf(gene_outputs(:,ind)))
        fitness=Inf;
        return
    end
end
gene_outputs=gene_outputs(:,2:end); % SEM coluna de 1's

%% Treino
if strcmp(gp.runcontrol.dataset,'treino')==1
 
    % Estimação 1
    if gp.runcontrol.tipo==0 % REGRESSÃO
        mdl=fitlm(gene_outputs,y);
        pvalue=mdl.Coefficients.pValue;
        
        
        
        
        
%         % Preparação
%         gene_outputs=[ones(length(y),1) gene_outputs]; % +[1]
%         prj=gene_outputs'*gene_outputs; % X'X
%         produto_invertido=pinv(prj); % Inversa(X'X)
%         beta=produto_invertido*gene_outputs'*y; % Beta
%         ypredtrain=gene_outputs*beta; % Predição
%         
%         % Variância de White
%         erros=y-ypredtrain;
%         erros_ao2=(erros).^2;
%         matriz_sigma=diag(erros_ao2);
%         matriz_variancia_beta=produto_invertido*(gene_outputs')*matriz_sigma*gene_outputs*produto_invertido;
%         variancia_beta=diag(matriz_variancia_beta);
%         
%         % Estatísticas de Teste t-Student e t(df,alfa/2)
%         % beta_temp=gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes};
%         beta_temp=mdl.Coefficients.Estimate;
%         estat_t=beta_temp./sqrt(variancia_beta);
%         estat_ttemp=abs(estat_t);
%         num_reg=gp.geracoes.num_reg(gp.geracoes.pop_size,gp.geracoes.num_geracoes);
%         limite_t_pos=tinv(1-(alfa/2),num_data_points-num_reg-1);
%         
%         th_HOMO=(pvalue<=0.05); % 0 = aceita H0
%         th_HETERO=(estat_ttemp>=limite_t_pos);
%         res=sum(th_HOMO==th_HETERO)./length(th_HOMO)
%         gp.TH.tabela(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=sum(th_HOMO==th_HETERO)./length(th_HOMO);
%                 
%         gene_outputs=gene_outputs(:,2:end); % SEM coluna de 1's
        
        
        
        
        
    elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
        [theta,~,stats]=mnrfit(gene_outputs,y);
        pvalue=stats.p;
    end
    
    if (gp.runcontrol.tipo==0)||(gp.runcontrol.num_classes==2)
        % Regressores Significantes
        pvalue=pvalue(2:end,:);
        reg_signif=pvalue<alfa;
        gp.geracoes.geneout_teste{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=reg_signif;
        gene_outputs=gene_outputs(:,reg_signif);

        % Regressores Significantes: Quem são
        individuos=gp.geracoes.individuos{gp.geracoes.pop_size,gp.geracoes.num_geracoes}; % R: todos
        gp.geracoes.regressores{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=individuos(reg_signif); % R: signif
        a1=strjoin(individuos(reg_signif));
        a2={'y','c',a1};
        gp.geracoes.eviews{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=strjoin(a2);

        % Regressores Significantes: NÚMERO
        numreg_signif=sum(reg_signif);
        gp.geracoes.numreg_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=numreg_signif;

        % Regressores Significantes: PROPORÇÃO
        prop_numreg_signif=numreg_signif/length(pvalue);
        gp.geracoes.prop_numreg_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=prop_numreg_signif;

        % Regressores Significantes: MÉDIA
        gp.geracoes.media_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mean(pvalue(reg_signif));
    end

    % Estimação 2
    if gp.runcontrol.tipo==0 % REGRESSÃO
        mdl=fitlm(gene_outputs,y);
        theta=mdl.Coefficients.Estimate;
        pvalue=mdl.Coefficients.pValue;
    elseif (gp.runcontrol.tipo==1)&&(gp.runcontrol.num_classes==2) % CLASSIFICAÇÃO
        [theta,~,stats]=mnrfit(gene_outputs,y);
        pvalue=stats.p;
    end

        % Alocação dos "#classes-1" vetores de betas
        gp.fitness.returnvalues{gp.state.current_individual}=theta;
        gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=theta;

    % Predição
    if gp.runcontrol.tipo==0 % REGRESSÃO
        ypredtrain=mdl.Fitted;
    elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
        ypred_prob=mnrval(theta,gene_outputs);
        [~,ypredtrain]=max(ypred_prob,[],2);
    end

    % TH e parâmetros de significância: Verificação
    if (gp.runcontrol.tipo==0)||(gp.runcontrol.num_classes==2)
        reg_signif=pvalue(2:end,:)<alfa;
        numreg_signif=sum(reg_signif);
        prop_numreg_signif=numreg_signif/(length(pvalue)-1);
        if prop_numreg_signif==1
            gp.geracoes.verif=gp.geracoes.verif+1;
            gp.geracoes.verif_2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=1;
        end
    end

    % Fitness
    if gp.runcontrol.tipo==0 % REGRESSÃO
        gp.geracoes.rmse(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mdl.RMSE;
        rmse=mdl.RMSE;
        fitness=mdl.RMSE;
        gp.geracoes.r2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mdl.Rsquared.Ordinary;
        gp.geracoes.r2ajustado(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mdl.Rsquared.Adjusted;
        disp('-------------------------------------------------------------------------------------');
    elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
        acertos=sum(y==ypredtrain);
        p_acertos=acertos/length(ypredtrain);
        % c=confusionmat(y,ypredtrain);
        % ber=0.5*((c(1,2)/(c(1,1)+c(1,2)))+(c(2,1)/(c(2,1)+c(2,2)))) % balanced error rate
        
        fitness=1-p_acertos;
        gp.geracoes.rmse(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=fitness;
        disp('-------------------------------------------------------------------------------------');
    end
  
elseif strcmp(gp.runcontrol.dataset,'teste')==1 % Teste
%% Teste

    % Estimação
    % theta=gp.fitness.returnvalues{gp.state.current_individual};
    theta=gp.geracoes.betas{gp.geracoes.ind_menoresRMSE(end),end};

    % Predição
    aux=gp.geracoes.geneout_teste(gp.geracoes.ind_menoresRMSE(1,end),end);
    aux=cell2mat(aux);
    gene_outputs=gene_outputs(:,aux);
    gene_outputs=[ones(length(y),1) gene_outputs];
    
    if gp.runcontrol.tipo==0 % REGRESSÃO
        ypredteste=gene_outputs*theta;
    elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
        ypred_prob=mnrval(theta,gene_outputs(:,2:end));
        [~,ypredteste]=max(ypred_prob,[],2);
    end

    % Fitness
    if gp.runcontrol.tipo==0 % REGRESSÃO
        rmse=sqrt(mean((y-ypredteste).^2));
        fitness=rmse;

            gp.geracoes.treinoEteste(2,1)=gp.geracoes.rmse(gp.geracoes.ind_menoresRMSE(1,end),end); % RMSE Treino
            gp.geracoes.treinoEteste(2,2)=fitness; % RMSE Teste

            num_reg=gp.geracoes.num_reg(gp.geracoes.ind_menoresRMSE(end),end);
            STE=sum((ypredteste-mean(y)).^2);
            STQ=sum((y-mean(y)).^2);
            r2=STE/STQ;
            r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

            gp.geracoes.treinoEteste(1,1)=gp.geracoes.r2ajustado(gp.geracoes.ind_menoresRMSE(1,end),end); % R2ajustado Treino
            gp.geracoes.treinoEteste(1,2)=r2ajustado; % R2ajustado Teste
            
    elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
        acertos=sum(y==ypredteste);
        p_acertos=acertos/length(ypredteste);
        % c=confusionmat(y,ypredteste);
        % ber=0.5*((c(1,2)/(c(1,1)+c(1,2)))+(c(2,1)/(c(2,1)+c(2,2))))
        fitness=p_acertos;

            gp.geracoes.treinoEteste(1,1)=1-gp.geracoes.rmse(gp.geracoes.ind_menoresRMSE(1,end),end); % Treino
            gp.geracoes.treinoEteste(1,2)=fitness; % Teste
            
            rmse=sqrt(mean((y-ypredteste).^2));
            gp.geracoes.treinoEteste(2,2)=rmse; % RMSE Teste
            
    end      
    
end

