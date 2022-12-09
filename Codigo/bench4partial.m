%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Derived Input Directions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark4: Partial Least Squares Regression (PLSR)
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
    num_comp=8;
    num_pcs=4;
elseif strcmp(data,'housing')==1 % data=='housing'
    gp.userdata.xtrain=Housing_Data(rand_TR_house,1:13);
    gp.userdata.ytrain=Housing_Data(rand_TR_house,14);
    gp.userdata.xtest=Housing_Data(rand_TE_house,1:13);
    gp.userdata.ytest=Housing_Data(rand_TE_house,14);
    num_comp=13;
    num_pcs=6;
elseif strcmp(data,'nasa')==1 % data=='nasa'
    gp.userdata.xtrain=Nasa_Data(rand_TR_nasa,1:5);
    gp.userdata.ytrain=Nasa_Data(rand_TR_nasa,6);
    gp.userdata.xtest=Nasa_Data(rand_TE_nasa,1:5);
    gp.userdata.ytest=Nasa_Data(rand_TE_nasa,6);
    num_comp=5;
    num_pcs=4;
end

%% Treino

X=gp.userdata.xtrain;
y=gp.userdata.ytrain;

[n,p] = size(X);
% [Xloadings,Yloadings,Xscores,Yscores,betaPLS,PLSPctVar]=plsregress(X,y,num_reg);
[XL,YL,XS,YS,betaPLS,PLSPctVar]=plsregress(X,y,num_comp);

% Para ver quantos PCs devo usar
plot(1:num_comp,cumsum(100*PLSPctVar(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');

% Resposta baseada no #PCs escolhido anteriormente
[Xloadings,Yloadings,Xscores,Yscores,betaPLS]=plsregress(X,y,num_pcs);
yfitPLS=[ones(n,1) X]*betaPLS;

% RMSE
ypredtrain=yfitPLS;
rmse=sqrt(mean((y-ypredtrain).^2));
b4partial.treino.rmse=rmse;

% R2
STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b4partial.treino.r2=r2;

% R2ajustado
num_data_points=size(y,1);
num_reg=length(betaPLS)-1;
b4partial.treino.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b4partial.treino.mad=mean(abs(y-ypredtrain));
b4partial.treino.mape=mean(abs((y-ypredtrain)./y))*100;

%% Teste

X=gp.userdata.xtest;
y=gp.userdata.ytest;

% RMSE
ypredteste=[ones(length(y),1) X]*betaPLS;
rmse=sqrt(mean((y-ypredteste).^2));
b4partial.teste.rmse=rmse;

% R2
STE=sum((ypredteste-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
r2=STE/STQ;
b4partial.teste.r2=r2;

% R2ajustado
num_data_points=size(y,1);
b4partial.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

% MAD e MAPE
b4partial.teste.mad=mean(abs(y-ypredteste));
b4partial.teste.mape=mean(abs((y-ypredteste)./y))*100;

%% Treino e Teste

b4partial.treinoEteste(1,1)=b4partial.treino.r2ajustado;
b4partial.treinoEteste(1,2)=b4partial.teste.r2ajustado;

b4partial.treinoEteste(2,1)=b4partial.treino.r2;
b4partial.treinoEteste(2,2)=b4partial.teste.r2;

b4partial.treinoEteste(3,1)=b4partial.treino.rmse;
b4partial.treinoEteste(3,2)=b4partial.teste.rmse;

b4partial.treinoEteste(4,1)=b4partial.treino.mad;
b4partial.treinoEteste(4,2)=b4partial.teste.mad;

b4partial.treinoEteste(5,1)=b4partial.treino.mape;
b4partial.treinoEteste(5,2)=b4partial.teste.mape;

save b4partial;
