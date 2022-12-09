function [gp]=gpdemo4_config(gp)
%% Parâmetros

% REGRESSÃO (0) OU CLASSIFICAÇÃO (1) ?
gp.runcontrol.tipo=1;

% pop_size e num_gen
gp.runcontrol.pop_size=100;   % 100
gp.runcontrol.num_gen=100;    % 50
gp.treedef.max_depth=3;       % 4 ou 5; maximum depth of trees
gp.genes.max_genes=3;         % 3; absolute maximum #genes per individual

gp.treedef.max_mutate_depth=gp.treedef.max_depth;

%% ESCOLHAS 1

% MULTIGENE ?
gp.runcontrol.num_classes=2; % 2 (não mexer)
gp.genes.multigene=true;   % true = indivíduo MULTI , false = indivíduo SINGLE

% PRESSSÃO LEXICOGRÁFICA
gp.selection.tournament.lex_pressure=true; % True to use Luke & Panait's plain lexicographic tournament selection
                                           % popbuild -> selection                                        
% RODA Benchmark? (Não mexer)
gp.benchmarks.opcao=0; % 0 = Não
gp.benchmarks.estat_signif=1; % 1 = procurará só estat_signif
gp.benchmarks.achou=0; % Não mexer

%% Outros

% Display on Screen
gp.runcontrol.verbose=5;      % Set to n to display run information to screen every n generations  

% Selection method options
gp.selection.tournament.size=7; % 7 = Típico

% Fitness function specification
gp.fitness.fitfun=@regressmulti_fitfun;  % Function handle to name of the user's fitness function (filename with no .m extension).

gp.fitness.minimisation=true; % Set to TRUE if you want to minimise the fitness function (if false, it is maximised).

gp.fitness.terminate=false;  %terminate run if fitness below acheived
gp.fitness.terminate_value=0.05;

%%
gp.userdata.datasampling = false;
gp.userdata.user_fcn=@regressmulti_fitfun_validate; %enables hold out validation set

% INPUT CONFIGURATION
% This sets the number of inputs (i.e. the size of the terminal set NOT including constants)
% Pega #colunas de xtrain; size = [#linhas #colunas]
gp.nodes.inputs.num_inp=size(gp.userdata.xtrain,2);

% DEFINE FUNCTIONS
%   (Below are some definitions of functions that have been used for symbolic regression problems)
%
%         Function name                                     Number of arguments
%   (must be an mfile on the path)
          
gp.nodes.functions.name{1}='times';                 % recebe 5 de exist (função do Matlab)
gp.nodes.functions.name{2}='minus';                 % receberia 2 de exist de fosse criada    
gp.nodes.functions.name{3}='plus';                    
gp.nodes.functions.name{4}='rdivide';               % unprotected divide (may cause NaNs)
gp.nodes.functions.name{5}='psqroot';               % protected sqrt
gp.nodes.functions.name{6}='plog';                  % protected natural log
gp.nodes.functions.name{7}='square';                % .^2 square
gp.nodes.functions.name{8}='tanh';                  % tanh function
gp.nodes.functions.name{9}='pdivide';               % protected divide function
gp.nodes.functions.name{10}='iflte';                % IF-THEN-ELSE function
gp.nodes.functions.name{11}='sin';                   
gp.nodes.functions.name{12}='cos';                  
gp.nodes.functions.name{13}='exp';                     

% ACTIVE FUNCTIONS
% Manually setting a function node to inactive allows you to exclude a function node in a particular run.

gp.nodes.functions.active(1)=1; % times
gp.nodes.functions.active(2)=1; % minus
gp.nodes.functions.active(3)=1; % plus
gp.nodes.functions.active(4)=0; % UNprotected DIVIDE
gp.nodes.functions.active(5)=0;
gp.nodes.functions.active(6)=0;
gp.nodes.functions.active(7)=0;
gp.nodes.functions.active(8)=0;
gp.nodes.functions.active(9)=0; % protected DIVIDE
gp.nodes.functions.active(10)=0;
gp.nodes.functions.active(11)=0;
gp.nodes.functions.active(12)=0;
gp.nodes.functions.active(13)=0;
