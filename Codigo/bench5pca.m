%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Derived Input Directions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark5: Principal Component Analysis/Regression (PCA/PCR)
%%

load dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% DATASET REAL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data='concrete'; % 'concrete', 'housing', 'nasa'

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
num_pcomp=4;
[n,p]=size(X);

[PCALoadings,PCAScores,PCAVar]=pca(X,'Economy',false);
betaPCR=regress(y-mean(y), PCAScores(:,1:num_pcomp));

betaPCR=PCALoadings(:,1:num_pcomp)*betaPCR;
betaPCR=[mean(y) - mean(X)*betaPCR; betaPCR];
yfitPCR=[ones(n,1) X]*betaPCR;

% RMSE
ypredtrain=yfitPCR;
rmse=sqrt(mean((y-ypredtrain).^2));
b5pca.treino.rmse=rmse;

% R2
STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b5pca.treino.r2=r2;

% R2ajustado
num_data_points=size(y,1);
num_reg=length(betaPCR)-1;
b5pca.treino.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b5pca.treino.mad=mean(abs(y-ypredtrain));
b5pca.treino.mape=mean(abs((y-ypredtrain)./y))*100;

%% Teste

X=gp.userdata.xtest;
y=gp.userdata.ytest;

% RMSE
ypredteste=[ones(length(y),1) X]*betaPCR;
rmse=sqrt(mean((y-ypredteste).^2));
b5pca.teste.rmse=rmse;

% R2
STE=sum((ypredteste-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b5pca.teste.r2=r2;

% R2ajustado
num_data_points=size(y,1);
b5pca.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b5pca.teste.mad=mean(abs(y-ypredteste));
b5pca.teste.mape=mean(abs((y-ypredteste)./y))*100;

%% Treino e Teste

b5pca.treinoEteste(1,1)=b5pca.treino.r2ajustado;
b5pca.treinoEteste(1,2)=b5pca.teste.r2ajustado;

b5pca.treinoEteste(2,1)=b5pca.treino.r2;
b5pca.treinoEteste(2,2)=b5pca.teste.r2;

b5pca.treinoEteste(3,1)=b5pca.treino.rmse;
b5pca.treinoEteste(3,2)=b5pca.teste.rmse;

b5pca.treinoEteste(4,1)=b5pca.treino.mad;
b5pca.treinoEteste(4,2)=b5pca.teste.mad;

b5pca.treinoEteste(5,1)=b5pca.treino.mape;
b5pca.treinoEteste(5,2)=b5pca.teste.mape;

save b5pca;
