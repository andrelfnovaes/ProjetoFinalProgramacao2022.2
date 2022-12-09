%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Subset Selection %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark1: Stepwise Regression
%%

load dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% DATASET REAL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data='nasa'; % 'concrete', 'housing', 'nasa'

if strcmp(data,'concrete')==1 % data=='concrete' 
    gp.userdata.xtrain=Concrete_Data(rand_TR_conc,1:8); % tr_ind = 0 ou 1 (rand)
    gp.userdata.ytrain=Concrete_Data(rand_TR_conc,9);
    gp.userdata.xtest=Concrete_Data(rand_TE_conc,1:8); % te_ind + tr_ind = vetor de 1's
    gp.userdata.ytest=Concrete_Data(rand_TE_conc,9); % te_ind = 0 ou 1 (rand)
elseif strcmp(data,'housing')==1 % data=='housing'
    gp.userdata.xtrain=Housing_Data(rand_TR_house,1:13);
    gp.userdata.ytrain=Housing_Data(rand_TR_house,14);
    gp.userdata.xtest=Housing_Data(rand_TE_house,1:13);
    gp.userdata.ytest=Housing_Data(rand_TE_house,14);
elseif strcmp(data,'nasa')==1 % data=='nasa'
    gp.userdata.xtrain=Nasa_Data(rand_TR_nasa,1:5);
    gp.userdata.ytrain=Nasa_Data(rand_TR_nasa,6);
    gp.userdata.xtest=Nasa_Data(rand_TE_nasa,1:5);
    gp.userdata.ytest=Nasa_Data(rand_TE_nasa,6);       
end

%% Treino

X=gp.userdata.xtrain;
y=gp.userdata.ytrain;

[betahat,se,pval,inmodel,stats,history]=stepwisefit(X,y,'penter',.05,'premove',0.10,'display','on');

% RMSE e R2
b1step.treino.rmse=stats.rmse;
b1step.treino.r2=1-(stats.SSresid/stats.SStotal);

% R2ajustado
r2=b1step.treino.r2;
num_data_points=size(y,1);
num_reg=sum(inmodel);

b1step.treino.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
betahat_modelo=betahat(inmodel);
X_modelo=X(:,inmodel);
theta=vertcat(stats.intercept,betahat_modelo);
gene_outputs=[ones(length(y),1) X_modelo];
ypredtrain=gene_outputs*theta;

b1step.treino.mad=mean(abs(y-ypredtrain));
b1step.treino.mape=mean(abs((y-ypredtrain)./y))*100;

%% Teste

X=gp.userdata.xtest;
y=gp.userdata.ytest;

X_modelo=X(:,inmodel);
gene_outputs=[ones(length(y),1) X_modelo];
ypredteste=gene_outputs*theta;

% RMSE e R2
rmse=sqrt(mean((y-ypredteste).^2));
b1step.teste.rmse=rmse;

STE=sum((ypredteste-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b1step.teste.r2=r2;

% R2ajustado
num_data_points=size(y,1);
b1step.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b1step.teste.mad=mean(abs(y-ypredteste));
b1step.teste.mape=mean(abs((y-ypredteste)./y))*100;

%% Treino e Teste

b1step.treinoEteste(1,1)=b1step.treino.r2ajustado;
b1step.treinoEteste(1,2)=b1step.teste.r2ajustado;

b1step.treinoEteste(2,1)=b1step.treino.r2;
b1step.treinoEteste(2,2)=b1step.teste.r2;

b1step.treinoEteste(3,1)=b1step.treino.rmse;
b1step.treinoEteste(3,2)=b1step.teste.rmse;

b1step.treinoEteste(4,1)=b1step.treino.mad;
b1step.treinoEteste(4,2)=b1step.teste.mad;

b1step.treinoEteste(5,1)=b1step.treino.mape;
b1step.treinoEteste(5,2)=b1step.teste.mape;

save b1step;
