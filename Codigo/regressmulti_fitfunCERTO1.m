function [fitness,gp,ypredtrain,fitness_test,ypredtest,pvals]=regressmulti_fitfun(evalstr,gp)
%REGRESSMULTI_FITFUN GPTIPS fitness function to perform multigene
%non-linear symbolic regression on data comprising one output y and
%multiple inputs x1, ..xn.
%
%   Fitness function for multigene symbolic regression.
%
%   [FITNESS,GP]=REGRESSMULTI_FITFUN(EVALSTR,GP) returns the FITNESS of
%   the symbolic expression(s) in the cell array EVALSTR using information
%   contained in the GP data structure. In this case FITNESS is the root
%   mean squared prediction error on the training data set. (FITNESS = raiz(RMSE))
%
%   [FITNESS,GP,YPREDTRAIN,FITNESS_TEST,YPREDTEST,PVALS]=REGRESSMULTI_FITFUN(EVALSTR,GP)
%   may be used post-run to compute the fitness value FITNESS_TEST on the test data set
%   as well as the prediction of the model on the training data YPREDTRAIN and the
%   testing data YPREDTEST. The statistical p-values are returned as PVALS
%   (PVALS only computed if the Statistics Toolbox is present, otherwise an
%   empty variable is returned).
%
%   Remarks:
%   Each observation of the response variable y is assumed to be a
%   non-linear function of the corresponding observations of the predictor
%   variables x1,..xn.
%
%   TRAINING data:
%   The user's GPTIPS configuration file should populate the following
%   required fields for the training data assuming 'Ntrain' observations on
%   the input and output data:
%   GP.USERDATA.XTRAIN should be a (Ntrain X n) matrix where the ith column
%   contains the Ntrain observations of the ith input variable xi.
%   GP.USERDATA.YTRAIN should be a (Ntrain x 1) vector containing the
%   corresponding observations of the response variable y.
%
%   TESTING data:
%   The following fields are optional and may be used, post-run, to see how
%   well evolved models generalise to an unseen test data set with Ntest
%   observations. They do not affect the model building process.
%   GP.USERDATA.XTEST should be a (Ntest X n) matrix where the ith column
%   contains the Ntest observations of the ith input variable xi.
%   GP.USERDATA.YTEST should be a (Ntest x 1) vector containing the
%   corresponding observations of the response variable y.
%
%
%   How multigene symbolic regression works:
%   In multigene symbolic regression, each prediction of y is formed by the
%   weighted output of each of the trees/genes in the multigene individual
%   plus a bias term. The number (M) and structure of the trees is evolved
%   automatically during a GPTIPS run (subject to user defined constraints).
%
%   i.e. ypredtrain = c0 + c1*tree1 + ... + cM*treeM
%
%   where c0 = bias term
%         c1,..,cM are the weights
%         M is the number of genes/trees comprising the current individual
%
%   The weights (i.e. regression coefficients) are automatically determined
%   by a least squares procedure for each multigene individual and are
%   stored in GP.FITNESS.RETURNVALUES for future use.
%
%
%   Note:
%   Because the GP structure is modified within this
%   function (i.e. the field GP.FITNESS.RETURNVALUES is used to store
%   the computed weighting coefficients for each gene) the GP structure
%   must be returned as an output argument.
%
%   This fitness function is used for multigene symbolic regression for
%   GPDEMO2, GPDEMO3 and GPDEMO4 (the configuration files for these are
%   GPDEMO2_CONFIG.M, GPDEMO3_CONFIG.M and GPDEMO4_CONFIG.M respectively)
%   but it can and should be used for the user's own non-linear regression
%   problems.
%
%
%   (c) Dominic Searson 2009
%
%   v1.0
%
%   See also REGRESSMULTI_FITFUN_VALIDATE, GPDEMO2_CONFIG, GPDEMO3_CONFIG,
%   GPDEMO4_CONFIG, GPDEMO2, GPDEM03, GPDEMO4




% process evalstr with regex to allow direct access to data matrices
pat='x(\d+)';

% if gp.userdata.scale % Se foi escalonado(normalizado)...
%     evalstr=regexprep(evalstr,pat,'gp.userdata.xtrainS(:,$1)') % regexprep = strrep (há 'pat' em 'evalstr'?) / FAZ UMA BAITA SUBSTITUIÇÃO AQUI! Entendi.
%     y=gp.userdata.ytrainS % ytrainS = ytrain SCALED!
% else
% %     evalstr=regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)')
% 
% %     temp_evalstr=evalstr % no fim, evalstr=temp_evalstr
% %     evalstr=gppretty(gp,evalstr)
% %     evalstr
% 
%     y=gp.userdata.ytrain
% end

y=gp.userdata.ytrain;
num_data_points=size(y,1); % recebe #linhas de y
num_genes=length(evalstr);


%set up a matrix to store the tree outputs plus *A BIAS COLUMN OF ONES*
gene_outputs=ones(num_data_points,num_genes+1); % Tudo bem acima: é notação matricial


%eval each gene in the current individual
for i=1:num_genes
    ind=i+1;
    eval(['gene_outputs(:,ind)=' evalstr{i} ';']) % DESCOBRI!
    % LINHA ACIMA: atualiza GENE_OUTPUTS com os cálculos de EVALSTR 
    %check for nonsensical answers and break out early with an 'inf' if so
    if  any(isnan(gene_outputs(:,ind))) || any(isinf(gene_outputs(:,ind)))
        fitness=Inf;
        return
    end
%     save gp
end


if ~gp.state.run_completed %i.e. only calc. weighting coeffs during an actual run


    % if data sampling is enabled, only fit regression coeffs on random (per
    % generation) sub-sample of training data: ?????
    if gp.userdata.datasampling

        %if NEW GENERATION, then specify a new random subset of training data
        if gp.state.current_individual==1
            rand_vec = rand(num_data_points,1);
            gp.userdata.y_select = (rand_vec<=0.75); % Por que 0,75?
        end

        %get subset of training data: gene_outputs REDUZIDO!
        gene_outputs_sampled=gene_outputs(gp.userdata.y_select,:);

        %prepare LS matrix (Least Squared)
        prj=gene_outputs_sampled'*gene_outputs_sampled;

        %compute gene weights on subset only
        try % pinv = PSEUDO INVERSA (?)
            theta=pinv(prj)*gene_outputs_sampled'*y(gp.userdata.y_select);
        catch
            fitness=Inf;
            return
        end

    else % AQUI! AQUI! AQUI! AQUI! AQUI! AQUI! AQUI! AQUI!

        %prepare LS matrix
        prj=gene_outputs'*gene_outputs;
        produto_invertido=pinv(prj);

        %calculate coeffs using SVD least squares on full training data set
        try
            mqo=produto_invertido*gene_outputs'*y;
            % mv=mvregress(gene_outputs,y,'algorithm','mvn');
            theta=mqo; % theta=mqo ; theta=mean(mv,mqo)
            % razao=round(mqo./mv)
                        
        catch
            fitness=Inf;
            return
        end

    end

    %assign poor fitness if any NaN or Inf
    if any(isinf(theta)) || any(isnan(theta))
        fitness=Inf;
        return
    end

    %write coeffs to returnvalues field for storage
    gp.fitness.returnvalues{gp.state.current_individual}=theta;
    gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=theta;

else % if post-run, get stored coeffs from return value field

    theta=gp.fitness.returnvalues{gp.state.current_individual};
end

%calc. PREDICTION of full training data set using the estimated weights
ypredtrain=gene_outputs*theta;
num_reg=gp.geracoes.num_reg(gp.geracoes.pop_size,gp.geracoes.num_geracoes);





% whichstats = {'yhat','tstat'}
% stats = regstats(y,gene_outputs(:,2:end),'linear',whichstats)
% yhat = stats.yhat
% tstat = stats.tstat
% pvals=stats.tstat.pval
% mean(abs((ypredtrain-yhat)./yhat))
% SQR=sum(erros_ao2)
% sigma2=SQR/(num_data_points-num_reg-1)

% Matriz de Variâncias-Covariâncias de Beta
erros=y-ypredtrain;
erros_ao2=(erros).^2;
matriz_sigma=diag(erros_ao2);
matriz_variancia_beta=produto_invertido*(gene_outputs')*matriz_sigma*gene_outputs*produto_invertido;
variancia_beta=diag(matriz_variancia_beta);

% Estatísticas de Teste t-Student e t(df,alfa/2)
alfa=0.05; % Teste bicaudal: usar alfa/2 = 0.025
beta_temp=gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes};
estat_t=beta_temp./sqrt(variancia_beta);
estat_ttemp=abs(estat_t);
limite_t_pos=tinv(1-(alfa/2),num_data_points-num_reg-1);

% NÚMERO de regressores significantes (quantidade)
reg_signif=(estat_ttemp>limite_t_pos);
numreg_signif=size(find(reg_signif(2:end,:)==1))
gp.geracoes.numreg_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=numreg_signif(1,1);

% PROPORÇÃO de regressores significantes
prop_numreg_signif=numreg_signif(1,1)/num_reg
gp.geracoes.prop_numreg_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=prop_numreg_signif;

% MÉDIA de Significância
index=find(reg_signif(2:end,:)==1);
index=index+1;
media_signif=(sum(estat_ttemp(index)-limite_t_pos))/numreg_signif(1,1);
gp.geracoes.media_signif(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=media_signif;

% Somente Regressores Significativos
reg_signif2=reg_signif';
reg_signif2=reg_signif2(2:end);
gp.geracoes.index_reg_signif{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=reg_signif2;
evalstr_signif=evalstr(reg_signif2); % Acho que NÃO vou precisar.
gp.geracoes.evalstr_signif{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=evalstr_signif; % Acho que NÃO vou precisar.

% Regressores INsignificantes
reg_INsignif=(estat_ttemp<limite_t_pos);
reg_INsignif2=reg_INsignif';
reg_INsignif2=reg_INsignif2(2:end);
gp.geracoes.index_reg_INsignif{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=reg_INsignif2;

% UNSCALE before reporting fitness, if required
if gp.userdata.scale % Rever: sigmay e muy
    ypredtrain= (ypredtrain.*gp.userdata.sigmay)+gp.userdata.muy;
end

%% Fitness

%calculate RMSE (Root Mean Squared prediction Error)
rmse=sqrt(mean((gp.userdata.ytrain-ypredtrain).^2))
gp.geracoes.rmse(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=rmse;

% RMSE: Root Mean Square Error
% MSE:       Mean Square Error
% MAD: Mean Absolute Deviation
% MAPE: Mean Absolute Percent Error (é o único adimensional)

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ
disp('-------------------------------------------------------------------------------------');

mad=mean(abs(gp.userdata.ytrain-ypredtrain));
gp.geracoes.mad(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mad;
mape=mean(abs((gp.userdata.ytrain-ypredtrain)./gp.userdata.ytrain))*100;
gp.geracoes.mape(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=mape;

% if r2<=0.90 % r2<=0.90 rmse<=10
    gp.geracoes.r2(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=r2;
    % mdl = fitlm(gene_outputs(:,2:end),y);
    % mdl.Rsquared.Ordinary;

%     n_reg=size(evalstr);
%     num_reg=n_reg(1,2);
%     gp.geracoes.num_reg(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=num_reg;
    r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));
    gp.geracoes.r2ajustado(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=r2ajustado;
% else
%     r2ajustado=0;
%     gp.geracoes.r2ajustado(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=0;
% end

fitness=rmse;

%% COMENTEI OS TESTES

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
