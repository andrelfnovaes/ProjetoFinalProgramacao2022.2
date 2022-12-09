%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Shrinkage Methods %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark2: Ridge Regression
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
parameter=1; % Quanto maior, maior o efeito Shrinkage
results=ridge(y,X,parameter);

% RMSE
beta_0=mean(y)-(results.beta'*mean(X)');
beta=vertcat(beta_0, results.beta); % beta = antigo theta

gene_outputs=[ones(length(y),1) X];
ypredtrain=gene_outputs*beta;

rmse=sqrt(mean((y-ypredtrain).^2));
b2ridge.treino.rmse=rmse;

% R2 e R2ajustado
b2ridge.treino.r2=results.rsqr;
b2ridge.treino.r2ajustado=results.rbar;

% MAD e MAPE
b2ridge.treino.mad=mean(abs(y-ypredtrain));
b2ridge.treino.mape=mean(abs((y-ypredtrain)./y))*100;

%% Teste

X=gp.userdata.xtest;
y=gp.userdata.ytest;

% RMSE
gene_outputs=[ones(length(y),1) X];
ypredteste=gene_outputs*beta;

rmse=sqrt(mean((y-ypredteste).^2));
b2ridge.teste.rmse=rmse;

% R2 e R2ajustado
STE=sum((ypredteste-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b2ridge.teste.r2=r2;

num_data_points=size(y,1);
num_reg=results.nvar;
b2ridge.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b2ridge.teste.mad=mean(abs(y-ypredteste));
b2ridge.teste.mape=mean(abs((y-ypredteste)./y))*100;

%% Treino e Teste

b2ridge.treinoEteste(1,1)=b2ridge.treino.r2ajustado;
b2ridge.treinoEteste(1,2)=b2ridge.teste.r2ajustado;

b2ridge.treinoEteste(2,1)=b2ridge.treino.r2;
b2ridge.treinoEteste(2,2)=b2ridge.teste.r2;

b2ridge.treinoEteste(3,1)=b2ridge.treino.rmse;
b2ridge.treinoEteste(3,2)=b2ridge.teste.rmse;

b2ridge.treinoEteste(4,1)=b2ridge.treino.mad;
b2ridge.treinoEteste(4,2)=b2ridge.teste.mad;

b2ridge.treinoEteste(5,1)=b2ridge.treino.mape;
b2ridge.treinoEteste(5,2)=b2ridge.teste.mape;

save b2ridge;
