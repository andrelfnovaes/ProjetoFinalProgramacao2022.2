function [gp]=fitness_correto(gp)
%% Normalização do RMSE

% (gp.geracoes.pop_size,gp.geracoes.num_geracoes)
gp.geracoes.rmse(:,gp.geracoes.num_geracoes)=gp.geracoes.fitness(:,gp.geracoes.num_geracoes);

rmse_max=max(gp.geracoes.rmse);
rmse_min=min(gp.geracoes.rmse);
range=rmse_max(:,gp.geracoes.num_geracoes)-rmse_min(:,gp.geracoes.num_geracoes);

fitness_correto=zeros(gp.runcontrol.pop_size,1);

j=gp.geracoes.num_geracoes;
for i=1:gp.runcontrol.pop_size    
    fitness_correto(i,1)=(gp.geracoes.rmse(i,j)-rmse_min(1,j))/range;
end

gp.geracoes.rmse_norm(:,gp.geracoes.num_geracoes)=fitness_correto

%% Agregação de Objetivos

% igual=1/3

peso_T1 = 0.01;
peso_T2 = peso_T1;
peso_rmse = 1-(peso_T1+peso_T2);

% gp.geracoes.fitness(:,gp.geracoes.num_geracoes)=(peso_T1*gp.geracoes.pvalor_T1(:,gp.geracoes.num_geracoes))+(peso_T2*gp.geracoes.pvalor_T2(:,gp.geracoes.num_geracoes))+(peso_rmse*gp.geracoes.rmse_norm(:,gp.geracoes.num_geracoes))

gp.geracoes.fitness(:,gp.geracoes.num_geracoes)=(peso_rmse*gp.geracoes.rmse_norm(:,gp.geracoes.num_geracoes))+((peso_T1+peso_T2)*(1-((gp.geracoes.pvalor_T1(:,gp.geracoes.num_geracoes)+gp.geracoes.pvalor_T2(:,gp.geracoes.num_geracoes))/2)));

%% Indivíduos de Gauss-Markov (identificação)

% testes=gp.geracoes.teste1(:,gp.geracoes.num_geracoes)+gp.geracoes.teste2(:,gp.geracoes.num_geracoes)
% gp.geracoes.IGM(:,gp.geracoes.num_geracoes)=find(testes==2);

end