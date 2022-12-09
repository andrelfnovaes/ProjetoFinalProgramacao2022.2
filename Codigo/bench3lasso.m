%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Shrinkage Methods %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark3: Lasso
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
[B,FitInfo]=lasso(X,y,'NumLambda',20);








