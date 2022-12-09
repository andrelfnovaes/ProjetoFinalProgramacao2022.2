function [gp,ypredtrain,fitness_test,ypredtest,pvals]=regressmulti_fitfun_teste(evalstr,gp)
%REGRESSMULTI_FITFUN GPTIPS fitness function to perform multigene

y=gp.userdata.ytest;
num_data_points=size(y,1);
num_genes=length(evalstr);
gene_outputs=ones(num_data_points,num_genes+1);

for i=1:num_genes
    ind=i+1;
    eval(['gene_outputs(:,ind)=' evalstr{i} ';'])
    if  any(isnan(gene_outputs(:,ind))) || any(isinf(gene_outputs(:,ind)))
        fitness=Inf;
        return
    end
end

%% Estimação
theta=gp.fitness.returnvalues{gp.state.current_individual};

%% Predição
ypredtrain=gene_outputs*theta;

%% Fitness

%calculate RMSE (Root Mean Squared prediction Error)
rmse=sqrt(mean((y-ypredtrain).^2));
% gp.geracoes.rmse(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=rmse;
rmse_teste=rmse;
gp.geracoes.treinoEteste(3,2)=rmse_teste; % RMSE Teste
gp.geracoes.treinoEteste(3,1)=gp.geracoes.rmse(gp.geracoes.ind_menoresRMSE(1,end),end); % RMSE Treino

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;

gp.geracoes.treinoEteste(2,2)=r2; % R2 Teste
gp.geracoes.treinoEteste(2,1)=gp.geracoes.r2(gp.geracoes.ind_menoresRMSE(1,end),end); % R2 Treino

num_reg=gp.geracoes.num_reg(gp.geracoes.ind_menoresRMSE(end),end);

r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));
r2ajustado_teste=r2ajustado;
gp.geracoes.treinoEteste(1,2)=r2ajustado_teste; % R2ajustado Teste
gp.geracoes.treinoEteste(1,1)=gp.geracoes.r2ajustado(gp.geracoes.ind_menoresRMSE(1,end),end); % R2ajustado Treino

mad=mean(abs(y-ypredtrain));
gp.geracoes.treinoEteste(4,2)=mad; % MAD Teste
gp.geracoes.treinoEteste(4,1)=gp.geracoes.mad(gp.geracoes.ind_menoresRMSE(1,end),end); % MAD Treino

mape=mean(abs((y-ypredtrain)./y))*100;
gp.geracoes.treinoEteste(5,2)=mape; % MAPE Teste
gp.geracoes.treinoEteste(5,1)=gp.geracoes.mape(gp.geracoes.ind_menoresRMSE(1,end),end); % MAPE Treino

%% Comentei os TESTES 1 e 2

% %% TESTE 1: Jarque-Bera (Normalidade)
% alfa=0.01;
% % Se jbtest(erro)=0, ACEITA Normalidade(H0)
% erro=gp.userdata.ytrain-ypredtrain;
% % hist(erro)
% [gp.geracoes.teste1(gp.geracoes.pop_size,gp.geracoes.num_geracoes),pvalor_T1]=jbtest(erro,alfa);
% gp.geracoes.pvalor_T1(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=pvalor_T1;
% 
% % TROCA: Se ACEITA Normalidade(H0), recebe '1'.
% if gp.geracoes.pvalor_T1(gp.geracoes.pop_size,gp.geracoes.num_geracoes)>=alfa
%     gp.geracoes.teste1(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=1;
% else
%     gp.geracoes.teste1(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=0;
% end
% 
% %% TESTE 2: Breuch-Pagan (Heteroce)
% % Se pvalor_T2 > alfa, ACEITA HOMOcedasticidade
% % Se ACEITA HOMOcedasticidade, recebe '1'.
% pvalor_T2=bpagan(gp.userdata.ytrain,gene_outputs);
% gp.geracoes.pvalor_T2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=pvalor_T2;
% 
% if gp.geracoes.pvalor_T2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)>=(alfa/2)
%     gp.geracoes.teste2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=1;   
% else
%     gp.geracoes.teste2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=0;
% end

%%

%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%% FIM FIM FIM FIM %%%%%
%%

%---- below is code for post run evaluation of individuals, it is not executed during a GPTIPS run----------

if gp.state.run_completed
    gp.userdata.showgraphs=true;
elseif gp.state.count==1
    gp.userdata.showgraphs=false;
end
% Entra aqui em runtree
if gp.userdata.showgraphs %%% PULA! Só entra se é o fim! %%%

    %compute variation explained (R2_train) for training data
    varexp_train=100*(1- sum( (gp.userdata.ytrain-ypredtrain).^2 )/sum( (gp.userdata.ytrain-mean(gp.userdata.ytrain)).^2 ) );

    plot_validation=0;


    %first,check if validation data is present. If so, need to plot that too
    if (isfield(gp.userdata,'xval')) && (isfield(gp.userdata,'yval')) && ...
            ~isempty(gp.userdata.xval) && ~isempty(gp.userdata.yval)


        plot_validation=1;

        evalstr=strrep(evalstr,'.xtrain','.xval');

        if gp.userdata.scale
            yval=gp.userdata.yvalS;
        else
            yval=gp.userdata.yval;
        end


        num_data_points=length(yval);

        %set up a matrix to store the tree outputs plus a bias column of ones
        gene_outputs_val=zeros(num_data_points,num_genes+1);
        gene_outputs_val(:,1)=ones;


        %eval each tree
        for i=1:num_genes
            ind=i+1;
            eval(['gene_outputs_val(:,ind)=' evalstr{i} ';'])
        end % atualiza GENE_OUTPUTS com os cálculos de EVALSTR

        ypredval=gene_outputs_val*theta; %create the prediction on the validation data

        %unscale for fitness and stats reporting
        if gp.userdata.scale
            ypredval= (ypredval.*gp.userdata.sigmay)+gp.userdata.muy;
        end

        fitness_val=sqrt(mean((gp.userdata.yval-ypredval).^2)); % RMSE de val 

        %compute variation explained for validation data
        varexp_val=100*(1- sum( (gp.userdata.yval-ypredval).^2 )/sum( (gp.userdata.yval-mean(gp.userdata.yval)).^2 ) );

        evalstr=strrep(evalstr,'.xval','.xtrain');
    end






    %generate prediction on test data (if present)

    if (isfield(gp.userdata,'xtest')) && (isfield(gp.userdata,'ytest')) && ...
            ~isempty(gp.userdata.xtest) && ~isempty(gp.userdata.ytest)


        evalstr=strrep(evalstr,'.xtrain','.xtest');


        if gp.userdata.scale
            ytest=gp.userdata.ytestS;
        else
            ytest=gp.userdata.ytest;
        end

        num_data_points=length(ytest);

        %set up a matrix to store the tree outputs plus a bias column of ones
        gene_outputs_test=zeros(num_data_points,num_genes+1);
        gene_outputs_test(:,1)=ones;


        %eval each tree
        for i=1:num_genes
            ind=i+1;
            eval(['gene_outputs_test(:,ind)=' evalstr{i} ';'])
        end

        ypredtest=gene_outputs_test*theta; %create the prediction  on the testing data

        %now unscale  for plotting and stats reporting
        if gp.userdata.scale
            ypredtest= (ypredtest.*gp.userdata.sigmay)+gp.userdata.muy;
        end

        fitness_test=sqrt(mean((gp.userdata.ytest-ypredtest).^2));

        %compute variation explained for test data
        varexp_test=100*(1- sum( (gp.userdata.ytest-ypredtest).^2 )/sum( (gp.userdata.ytest-mean(gp.userdata.ytest)).^2 ) );

        %%% Gráficos %%%
        % model prediction
        f=figure('name','GPTIPS Multigene regression. Model prediction of individual.','numbertitle','off');
        subplot(2+plot_validation,1,1)
        plot(ypredtrain,'r')
        hold on
        plot(gp.userdata.ytrain)
        ylabel('y')
        xlabel('Data point')
        legend('Predicted y (training values)','Actual y (training values)')
        title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])
        hold off


        subplot(2+plot_validation,1,2)
        plot(ypredtest,'r')
        hold on
        plot(gp.userdata.ytest)
        ylabel('y')
        xlabel('Data point')
        legend('Predicted y (test values)','Actual y (test values)')
        title(['RMS test set error: ' num2str(fitness_test) ' Variation explained: ' num2str(varexp_test) ' %'])
        hold off


        %scatterplot
        s=figure('name','GPTIPS Multigene regression. Prediction scatterplot of individual.','numbertitle','off');
        subplot(2+plot_validation,1,1)
        minval=min([gp.userdata.ytrain;ypredtrain]);
        maxval=max([gp.userdata.ytrain;ypredtrain]);
        scatter(gp.userdata.ytrain,ypredtrain)
        axis ([minval maxval minval maxval])
        l1=line([minval maxval], [minval maxval]);
        set(l1,'color','black')
        box on;grid on
        ylabel('Predicted')
        xlabel('Actual')
        title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])



        subplot(2+plot_validation,1,2)
        minval=min([gp.userdata.ytest;ypredtest]);
        maxval=max([gp.userdata.ytest;ypredtest]);
        scatter(gp.userdata.ytest,ypredtest)
        axis ([minval maxval minval maxval])
        l2=line([minval maxval], [minval maxval]);
        set(l2,'color','black')
        box on;grid on
        ylabel('Predicted')
        xlabel('Actual')
        title(['RMS test set error: ' num2str(fitness_test) ' Variation explained: ' num2str(varexp_test) ' %'])

        if plot_validation

            figure(f)
            subplot(3,1,3)
            plot(ypredval,'r')
            hold on
            plot(gp.userdata.yval)
            ylabel('y')
            xlabel('Data point')
            legend('Predicted y (validation values)','Actual y (validation values)')
            title(['RMS validation set error: ' num2str(fitness_val) ' Variation explained: ' num2str(varexp_val) ' %'])
            hold off

            figure(s)
            subplot(3,1,3)
            minval=min([gp.userdata.yval;ypredval]);
            maxval=max([gp.userdata.yval;ypredval]);
            scatter(gp.userdata.yval,ypredval)
            axis ([minval maxval minval maxval])
            l3=line([minval maxval], [minval maxval]);
            set(l3,'color','black')
            box on;grid on
            title(['RMS validation set error: ' num2str(fitness_val) ' Variation explained: ' num2str(varexp_val) ' %'])
            ylabel('Predicted')
            xlabel('Actual')



        end


    else %if no test data just show training data


        if plot_validation


            f=figure('name','GPTIPS Multigene regression. Model prediction of individual.','numbertitle','off');
            train_ax=subplot(2,1,1);
            plot(train_ax,ypredtrain,'r')
            hold on
            plot(train_ax,gp.userdata.ytrain)
            ylabel('y')
            xlabel('Data point')
            legend('Predicted y (training values)','Actual y (training values)')
            title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])
            hold off

            val_ax=subplot(2,1,2);
            plot(val_ax,ypredval,'r')
            hold on
            plot(val_ax,gp.userdata.yval)
            ylabel('y')
            xlabel('Data point')
            legend('Predicted y (validation values)','Actual y (validation values)')
            title(['RMS validation set error: ' num2str(fitness_val) ' Variation explained: ' num2str(varexp_val) ' %'])
            hold off

            %scatterplot
            figure('name','GPTIPS Multigene regression. Prediction scatterplot of individual.','numbertitle','off')
            subplot(2,1,1)
            minval=min([gp.userdata.ytrain;ypredtrain]);
            maxval=max([gp.userdata.ytrain;ypredtrain]);
            scatter(gp.userdata.ytrain,ypredtrain)
            axis ([minval maxval minval maxval])
            l1=line([minval maxval], [minval maxval]);
            set(l1,'color','black')
            box on;grid on
            ylabel('Predicted')
            xlabel('Actual')
            title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])



            subplot(2,1,2)
            minval=min([gp.userdata.yval;ypredval]);
            maxval=max([gp.userdata.yval;ypredval]);
            scatter(gp.userdata.yval,ypredval)
            axis ([minval maxval minval maxval])
            l3=line([minval maxval], [minval maxval]);
            set(l3,'color','black')
            box on;grid on
            title(['RMS validation set error: ' num2str(fitness_val) ' Variation explained: ' num2str(varexp_val) ' %'])
            ylabel('Predicted')
            xlabel('Actual')

            disp('No test set data found.')



        else



            f=figure('name','GPTIPS Multigene regression. Model prediction of individual.','numbertitle','off');
            plot(ypredtrain,'r')
            hold on
            plot(gp.userdata.ytrain)
            ylabel('y')
            xlabel('Data point')
            legend('Predicted y (training values)','Actual y (training values)')
            title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])
            hold off


            %scatterplot
            figure('name','GPTIPS Multigene regression. Prediction scatterplot of individual.','numbertitle','off')
            minval=min([gp.userdata.ytrain;ypredtrain]);
            maxval=max([gp.userdata.ytrain;ypredtrain]);
            scatter(gp.userdata.ytrain,ypredtrain)
            axis ([minval maxval minval maxval])
            l1=line([minval maxval], [minval maxval]);
            set(l1,'color','black')
            box on;grid on
            ylabel('Predicted')
            xlabel('Actual')
            title(['RMS training set error: ' num2str(fitness) ' Variation explained: ' num2str(varexp_train) ' %'])


            disp('No test set data found.')


        end

        fitness_test=[];
        ypredtest=[];

    end




    %Display statistical analysis of model term significance (if stats
    %toolbox is present and graphs are enabled)

    if license('test','statistics_toolbox') % TENHO O PACOTE!
        % Regress tree outputs (and bias) against y train data and get stats
        stats=regstats(y,gene_outputs(:,2:end));
        pvals=stats.tstat.pval;
    else
        pvals=[];
    end

    if  license('test','statistics_toolbox')



        %generate x labels for bar graphs
        gene_labels={'Bias'};
        for i=1:num_genes
            gene_labels{i+1}=['Gene ' int2str(i)];
        end

        %plot gene weights and offset
        statfig=figure;
        coeffs_ax=subplot(2,1,1);
        set(statfig,'name','Statistical properties of multigene model (on training data)','numbertitle','off')
        bar(coeffs_ax,stats.beta); shading faceted
        set(coeffs_ax,'xtick',1:(num_genes+1))
        set(coeffs_ax,'xticklabel',gene_labels)
        title(coeffs_ax,'Gene weights')



        %plot p-vals
        pvals_ax=subplot(2,1,2);
        bar(pvals_ax,stats.tstat.pval); shading faceted
        set(pvals_ax,'xtick',1:(num_genes+1))
        set(pvals_ax,'xticklabel',gene_labels)
        title(pvals_ax,'P value (low = significant)')
        xlabel(['R squared = ' num2str(stats.rsquare) ' Adj. R squared = ' num2str(stats.adjrsquare)])

    end
    figure(f)
else
    ypredtest=[];
    pvals=[];
    fitness_test=[];
end
