function []=eval_final()

% Calcula Grandezas (R2ajustado, R2, RMSE, MAD, MAPE) para um
% dataset referenciado (Treino ou Teste), com todos os regressores
% ou somente com os Estatisticamente Significativos.

load gp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Quem é o indivíduo? %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
individuo='2.0*x5 - 1.0*x3 - 1.0*x4 - 1.0*x2 - 3.0*x1*x3 - 1.0*x2*x3 - 1.0*x2*x4 + 2.0*x2*x5 - 1.0*x3*x4 + 2.0*x4*x5 + x1*x3^2 - 1.0*x3*x5^2 - 2.0*x4^2*x5 - 1.0*x2^2 - 1.0*x1*x3*x5^2 + 2.0*x3*x4*x5^2 - 1.0*x1*x3*x5 + x2*x3*x5 - 2.0*x2*x4*x5 - 1.0*x1*x2*x3*x5 + x1*x3*x4*x5 + x2*x3*x4*x5 + 1.0';

[evalstr,num_reg]=regressores(individuo);
temp_evalstr=evalstr;

%% COMPLETO: Treino

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Treino %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=gp.userdata.ytrain;

pat='x(\d+)';
evalstr=regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)');
pat='^';
evalstr=strrep(evalstr,pat,'.^');
pat='*';
evalstr=strrep(evalstr,pat,'.*');
evalstr_treino=evalstr;

num_data_points=size(y,1); % recebe #linhas de y
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

prj=gene_outputs'*gene_outputs;
produto_invertido=pinv(prj);
mqo=produto_invertido*gene_outputs'*y;
theta=mqo;
ypredtrain=gene_outputs*theta;

% Matriz de Variâncias-Covariâncias de Beta
erros=y-ypredtrain;
erros_ao2=(erros).^2;
matriz_sigma=diag(erros_ao2);
matriz_variancia_beta=produto_invertido*(gene_outputs')*matriz_sigma*gene_outputs*produto_invertido;
variancia_beta=diag(matriz_variancia_beta);

% Estatísticas de Teste t-Student e t(df,alfa/2)
alfa=0.05; % Teste bicaudal: usar alfa/2 = 0.025
beta_temp=theta;
estat_t=beta_temp./sqrt(variancia_beta);
estat_ttemp=abs(estat_t);
limite_t_pos=tinv(1-(alfa/2),num_data_points-num_reg-1);

% ÍNDICE dos significantes (NÚMERO de regressores significantes)
reg_signif=(estat_ttemp>limite_t_pos);
reg_signif2=reg_signif';
reg_signif2=reg_signif2(2:end);

% R2ajustado, R2, RMSE, MAD, MAPE
gp.final.completo.treino.rmse=sqrt(mean((gp.userdata.ytrain-ypredtrain).^2))

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
gp.final.completo.treino.r2=STE/STQ

gp.final.completo.treino.mad=mean(abs(gp.userdata.ytrain-ypredtrain));
gp.final.completo.treino.mape=mean(abs((gp.userdata.ytrain-ypredtrain)./gp.userdata.ytrain))*100;

r2=gp.final.completo.treino.r2;
gp.final.completo.treino.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

%% COMPLETO: Teste

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Teste %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=gp.userdata.ytest;

[evalstr,num_reg]=regressores(individuo);
pat='x(\d+)';
evalstr=regexprep(evalstr,pat,'gp.userdata.xtest(:,$1)');
pat='^';
evalstr=strrep(evalstr,pat,'.^');
pat='*';
evalstr=strrep(evalstr,pat,'.*');

num_data_points=size(y,1); % recebe #linhas de y
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

prj=gene_outputs'*gene_outputs;
produto_invertido=pinv(prj);
mqo=produto_invertido*gene_outputs'*y;
theta=mqo;
ypredtrain=gene_outputs*theta;

% R2ajustado, R2, RMSE, MAD, MAPE
gp.final.completo.teste.rmse=sqrt(mean((gp.userdata.ytest-ypredtrain).^2))

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
gp.final.completo.teste.r2=STE/STQ

gp.final.completo.teste.mad=mean(abs(gp.userdata.ytest-ypredtrain));
gp.final.completo.teste.mape=mean(abs((gp.userdata.ytest-ypredtrain)./gp.userdata.ytest))*100;

r2=gp.final.completo.teste.r2;
gp.final.completo.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

%% INCOMPLETO: Treino

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Treino %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=gp.userdata.ytrain;

gp.final.ES.individuo_ES=temp_evalstr(reg_signif2);
evalstr=evalstr_treino(reg_signif2);

num_data_points=size(y,1); % recebe #linhas de y
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

prj=gene_outputs'*gene_outputs;
produto_invertido=pinv(prj);
mqo=produto_invertido*gene_outputs'*y;
theta=mqo;
ypredtrain=gene_outputs*theta;

% R2ajustado, R2, RMSE, MAD, MAPE
gp.final.ES.treino.rmse=sqrt(mean((gp.userdata.ytrain-ypredtrain).^2))

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
gp.final.ES.treino.r2=STE/STQ

gp.final.ES.treino.mad=mean(abs(gp.userdata.ytrain-ypredtrain));
gp.final.ES.treino.mape=mean(abs((gp.userdata.ytrain-ypredtrain)./gp.userdata.ytrain))*100;

r2=gp.final.ES.treino.r2;
gp.final.ES.treino.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));


%% INCOMPLETO: Teste

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Teste %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=gp.userdata.ytest;

[evalstr2,num_reg]=regressores(individuo);
evalstr=evalstr2(reg_signif2);

pat='x(\d+)';
evalstr=regexprep(evalstr,pat,'gp.userdata.xtest(:,$1)');
pat='^';
evalstr=strrep(evalstr,pat,'.^');
pat='*';
evalstr=strrep(evalstr,pat,'.*');

num_data_points=size(y,1); % recebe #linhas de y
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

prj=gene_outputs'*gene_outputs;
produto_invertido=pinv(prj);
mqo=produto_invertido*gene_outputs'*y;
theta=mqo;
ypredtrain=gene_outputs*theta;

% R2ajustado, R2, RMSE, MAD, MAPE
gp.final.ES.teste.rmse=sqrt(mean((gp.userdata.ytest-ypredtrain).^2))

STE=sum((ypredtrain-mean(y)).^2);
STQ=sum((y-mean(y)).^2);
gp.final.ES.teste.r2=STE/STQ

gp.final.ES.teste.mad=mean(abs(gp.userdata.ytest-ypredtrain));
gp.final.ES.teste.mape=mean(abs((gp.userdata.ytest-ypredtrain)./gp.userdata.ytest))*100;

r2=gp.final.ES.teste.r2;
gp.final.ES.teste.r2ajustado=r2-((1-r2)*(num_reg/(num_data_points-num_reg-1)));

save gp;
end