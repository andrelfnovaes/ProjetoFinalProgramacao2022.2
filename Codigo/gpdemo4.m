function []=gpdemo4(data,dataset,porcaoTR,nCV,limiar)

% GPDEMO4 GPTIPS demo of multiple gene regression on a concrete compressive strength data set.

disp('Programação Genética Econométrica à Problemas de Regressão e Classificação');
% pause;

k=10; % k-fold CV

if strcmp(data,'sonar')==1
    k=13;
end

%% Dataset: organização

[L,C]=size(dataset);
linhasTR=round(porcaoTR*L); % "linhasTR=L", se é CV sem teste
GP.xtrain=dataset(1:linhasTR,1:(C-1));
GP.ytrain=dataset(1:linhasTR,C);
if porcaoTR~=1 % Há teste!
    GP.xtest=dataset((linhasTR+1):L,1:(C-1));
    GP.ytest=dataset((linhasTR+1):L,C);
end

clearvars -except data porcaoTR nCV limiar k L C linhasTR GP

%% K-fold CV

m=1; % Não mexer

for j=1:nCV % #experimentos
    for i=1:k % k-fold CV
        tic;
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

        % Call GP program using the configuration in gpdemo4_config.m
        gp=rungp('gpdemo4_config',gp);
        t=toc;
        gp.geracoes.tempo=t;

        if porcaoTR~=1  % Há teste!

            %%%%%%%%%%% TESTE %%%%%%%%%%%

            % Evalstr: somente regressores ES
            evalstr=gp.geracoes.regressores{gp.geracoes.ind_menoresRMSE(end),end};

            % Modificações necessárias
            pat='x(\d+)';
            evalstr=regexprep(evalstr,pat,'GP.xtest(:,$1)');
            pat='^';
            evalstr=strrep(evalstr,pat,'.^');
            pat='*';
            evalstr=strrep(evalstr,pat,'.*');

            y=GP.ytest;

            num_data_points=size(y,1); % #linhas de y
            num_genes=length(evalstr); % #regressores
            gene_outputs=ones(num_data_points,num_genes+1);
            for I=1:num_genes
                ind=I+1;
                eval(['gene_outputs(:,ind)=' evalstr{I} ';'])
                if  any(isnan(gene_outputs(:,ind))) || any(isinf(gene_outputs(:,ind)))
                    fitness=Inf;
                    return
                end
            end
            gene_outputs=gene_outputs(:,2:end); % SEM coluna de 1's

            % Estimação
            theta=gp.geracoes.betas{gp.geracoes.ind_menoresRMSE(end),end};

            % Predição
            gene_outputs=[ones(length(y),1) gene_outputs]; % coluna de 1's + colunas de ES

            if gp.runcontrol.tipo==0 % REGRESSÃO
                ypredteste=gene_outputs*theta;
            elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
                ypred_prob=mnrval(theta,gene_outputs(:,2:end));
                [~,ypredteste]=max(ypred_prob,[],2);
            end

            % Fitness
            if gp.runcontrol.tipo==0 % REGRESSÃO

                num_reg=gp.geracoes.num_reg(gp.geracoes.ind_menoresRMSE(end),end);
                STE=sum((ypredteste-mean(y)).^2);
                STQ=sum((y-mean(y)).^2);
                r2=STE/STQ;
                r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));
                fitness=r2ajustado;

            elseif gp.runcontrol.tipo==1 % CLASSIFICAÇÃO
                acertos=sum(y==ypredteste);
                p_acertos=acertos/length(ypredteste);
                % c=confusionmat(y,ypredteste);
                % ber=0.5*((c(1,2)/(c(1,1)+c(1,2)))+(c(2,1)/(c(2,1)+c(2,2))))
                fitness=p_acertos;

            end

            %%%%%%%%%%% %%%%% %%%%%%%%%%%

        end

        % Registro
        GP.treinoEvalEteste(m,[1 2])=gp.geracoes.treinoEteste(1,[1 2]);
        if porcaoTR~=1  % Há teste!
            GP.treinoEvalEteste(m,3)=fitness;
        end
        GP.treinoEvalEteste(m,4)=gp.geracoes.reg_signif_menoresRMSE(1,length(gp.geracoes.reg_signif_menoresRMSE));
        filename=[data num2str(m) '.mat' ];
        save(filename);

        % Verificação: achou o que quer?
        if porcaoTR~=1  % Há teste!
            GP.teste=GP.treinoEvalEteste(m,3);
            if GP.treinoEvalEteste(m,3)>=limiar
                disp('ACHOU! ACHOU! ACHOU!');
                filename=[data 'GP' '.mat' ];
                save(filename);
                return;
            end
        else            % NÃO Há teste!
            if rem(m,k)==0
                media=mean(GP.treinoEvalEteste((m-(k-1)):m,2));
                GP.dp=std(GP.treinoEvalEteste((m-(k-1)):m,2));
                GP.media=media;
                if media>=(limiar)
                    disp('ACHOU! ACHOU! ACHOU!');
                    filename=[data 'GP' '.mat' ];
                    save(filename);
                    return;
                end
            end 
        end
        m=m+1;
    end % k-fold CV
    
    filename=[data 'NAO' num2str(j) '.mat' ];
    save(filename,'GP');
    
%     % Deletar
%     for o=(m-k):(m-1)
%         filename=[data num2str(o) '.mat'];
%         delete(filename);
%     end
    
end %  nCV #experimentos

end % function
