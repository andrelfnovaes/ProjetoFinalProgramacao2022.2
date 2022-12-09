function gp=rungp(config_file,gp)

% if no config supplied then report error
if nargin<1
    disp('Cannot execute. To run GPTIPS a configuration m file must be ')
    disp('specified, e.g. gp=rungp(''quartic_poly'') ')
    return
end

% generate prepared data structure with some default parameter values
gp=gpdefaults(gp);

% run user configuration file
% FEVAL: feval(F,x1,...,xn) evaluates the function specified by a function
% handle or function name, F, at the given arguments, x1,...,xn.
gp=feval(config_file,gp);

% perform error checks
gp=gpcheck(gp);

% perform initialisation
gp=gpinit(gp);

% Set and store the the PRNG seed
gp.info.PRNGseed=sum(100*clock);
rand('twister',gp.info.PRNGseed)

%% Treino
gp.runcontrol.dataset='treino';
gp.geracoes.verif=0;

%%
% main generation loop
for count=1:gp.runcontrol.num_gen % for 1 : #gerações

    gp.geracoes.num_geracoes=count;

    if count==1;

        gp=initbuild(gp); % generate initial population
        gp=evalfitness(gp); % calculate fitnesses of population members

        if gp.benchmarks.achou==1
            return
        end

    else % Crossover, Mutação e GAP: Variantes
         % Pressão Lexográfica: Variantes

        num_inicial=round(gp.runcontrol.pop_size/1); %4
        
        % 1 UP
        if count>=1 && count<=round(0.30*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.05;

            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial;

        % 2 DOWN    
        elseif count>round(0.30*gp.runcontrol.num_gen) && count<=round(0.40*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.30;

            gp.geracoes.pressao.booleano=false;
            gp.geracoes.pressao.num_classes=num_inicial;


        % 3 UP
        elseif count>round(0.40*gp.runcontrol.num_gen) && count<=round(0.70*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.10;

            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial; % +2

        % 4 DOWN
        elseif count>round(0.70*gp.runcontrol.num_gen) && count<=round(0.80*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.35;

            gp.geracoes.pressao.booleano=false;
            gp.geracoes.pressao.num_classes=num_inicial;


        % 5 UP
        elseif count>round(0.80*gp.runcontrol.num_gen) && count<=round(0.90*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.15;

            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial;  % +4

        % 6 DOWN
        elseif count>round(0.90*gp.runcontrol.num_gen) && count<=round(1.00*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
%             gp.selection.elite_fraction=0.40;

            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial;  % +4

        end

        gp=popbuild(gp); % use crossover, mutation etc. to generate a new population
        gp=evalfitness(gp); % calculate fitnesses of population members

        if gp.benchmarks.achou==1
            return
        end

    end

end  %end generation loop

%% "Gráficos"

% Índices dos Melhores indivíduos (em todas as gerações)
menores=min(gp.geracoes.rmse);

if gp.runcontrol.tipo==0
    r2ajustado_menoresRMSE=[]; end
reg_menoresRMSE=[];
ind_menoresRMSE=[];
reg_signif_menoresRMSE=[];
propreg_signif_menoresRMSE=[];

for i=1:length(menores)
    ind_menor=find(gp.geracoes.rmse(:,i)==menores(i));
    if length(ind_menor)==1
        if gp.runcontrol.tipo==0
            r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)]; end
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
        reg_signif_menoresRMSE=[reg_signif_menoresRMSE gp.geracoes.numreg_signif(ind_menor,i)];
        propreg_signif_menoresRMSE=[propreg_signif_menoresRMSE gp.geracoes.prop_numreg_signif(ind_menor,i)];
    else
        ind_menor=ind_menor(1);
        if gp.runcontrol.tipo==0
            r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)]; end
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
        reg_signif_menoresRMSE=[reg_signif_menoresRMSE gp.geracoes.numreg_signif(ind_menor,i)];
        propreg_signif_menoresRMSE=[propreg_signif_menoresRMSE gp.geracoes.prop_numreg_signif(ind_menor,i)];
    end
end

if gp.runcontrol.tipo==0
    gp.geracoes.r2ajustado_menoresRMSE=r2ajustado_menoresRMSE; end
gp.geracoes.reg_menoresRMSE=reg_menoresRMSE;
gp.geracoes.ind_menoresRMSE=ind_menoresRMSE;
gp.geracoes.reg_signif_menoresRMSE=reg_signif_menoresRMSE;
gp.geracoes.propreg_signif_menoresRMSE=propreg_signif_menoresRMSE;

if gp.runcontrol.tipo==0
    gp.Bests.R2ajustado=gp.geracoes.r2ajustado_menoresRMSE; end
gp.Bests.rmse=menores;
gp.Bests.num_reg=gp.geracoes.reg_menoresRMSE;
gp.Bests.reg_signif=gp.geracoes.reg_signif_menoresRMSE;
gp.Bests.propreg_signif=gp.geracoes.propreg_signif_menoresRMSE;

%% Teste

gp.runcontrol.dataset='teste';
[gp]=evalfitness(gp);

% gp.geracoes.teste.r2ajustado_teste=r2ajustado_teste;
% gp.geracoes.teste.r2_teste=r2_teste;
% gp.geracoes.teste.rmse_teste=rmse_teste;

