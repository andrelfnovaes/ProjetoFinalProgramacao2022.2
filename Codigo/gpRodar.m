% gprodar

load dataset;
nCV=1; % #experimentos (1 exp = 1 k-fold)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% 
%%%%%%%%%%%% Classifica��o %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% appendicitis

%% wisconsin

% % Tentativa 1
% data='wisconsin';
% dataset=C_wisconsin;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.9750;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.9670;
% target 3 = 0.9650;

%% ljubljana ---------- RUIM ----------

% Tentativa 1
data='ljubljana';
dataset=C_ljubljana;
porcaoTR=1.0;           % 1 = N�O h� teste
limiar=0.7620;
gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.7300;
% target 3 = 0.7100;

%% hepatitis ---------- RUIM ----------

% % Tentativa 1
% data='hepatitis';
% dataset=C_hepatitis;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.9200;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.8800;
% target 3 = 0.8600;

%% cleveland ---------- RUIM ----------

% % Tentativa 1
% data='cleveland';
% dataset=C_cleveland;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.8500;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.8400;
% target 3 = 0.8300;

%% diabetes

% % Tentativa 1
% data='diabetes';
% dataset=C_diabetes;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.7770;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.7600;
% target 3 = 0.7550;

%% hypothyroid ---------- RUIM ----------

% % Tentativa 1
% data='hypothyroid';
% dataset=C_hypothyroid;
% porcaoTR=0.52388888;    % 1 = N�O h� teste
% limiar=0.9920;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.9900;
% target 3 = 0.9850;

%% ionosphere

% % Tentativa 1
% data='ionosphere';
% dataset=C_ionosphere;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.9800;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.9400;
% target 3 = 0.9300;

%% sonar (k=13)

% nCV=50;
% 
% % Tentativa 1
% data='sonar';
% dataset=C_sonar;
% porcaoTR=1.0;           % 1 = N�O h� teste
% limiar=0.8800;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

% target 2 = 0.8600;
% target 3 = 0.8400;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% 
%%%%%%%%%%%%%% Regress�o %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% concrete

% data='concrete';
% dataset=R_concrete;
% porcaoTR=0.7000;           % 1 = N�O h� teste
% limiar=10000;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

%% housing

% % Tentativa 1
% data='housing';
% dataset=R_housing;
% porcaoTR=0.7000;           % 1 = N�O h� teste
% limiar=10000;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

%% nasa

% % Tentativa 1
% data='nasa';
% dataset=R_nasa;
% porcaoTR=0.7000;           % 1 = N�O h� teste
% limiar=10000;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

%% rmsd

% % Tentativa 1
% data='rmsd';
% dataset=R_rmsd;
% porcaoTR=0.7000;           % 1 = N�O h� teste
% limiar=10000;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);

%% yacht

% % Tentativa 1
% data='yacht';
% dataset=R_yacht;
% porcaoTR=0.7000;           % 1 = N�O h� teste
% limiar=10000;
% gpdemo4(data,dataset,porcaoTR,nCV,limiar);
