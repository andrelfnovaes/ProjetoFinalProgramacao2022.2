function gp=rungp(config_file)
%RUNGP Runs GPTIPS using the specified parameter file.
%
%   GP=RUNGP('yourfile') or GP=RUNGP(@yourfile) runs GPTIPS using the
%   parameters contained in the file yourfile.m and returns the results in
%   the GP data structure.  Post run analysis commands may then be
%   run on this structure, e.g. SUMMARY(GP) or RUNTREE(GP,'BEST').
%
%   Demos:
%   Use 'gpdemo1' at the commmand line for a demonstration of simple symbolic
%   regression (not multigene).
%
%   For demos involving multigene symbolic regression see gpdemo2,
%   gpdemo3 and gpdemo4.
%
%    ----------------------------------------------------------------------
%
%    Copyright (C) 2010  Dominic Searson
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   ----------------------------------------------------------------------
%
%   See also: GPDEMO1, GPDEMO2, GPDEMO3, GPDEMO4
%
%   v1.0 10th January, 2010.

% if no config supplied then report error
if nargin<1
    disp('Cannot execute. To run GPTIPS a configuration m file must be ')
    disp('specified, e.g. gp=rungp(''quartic_poly'') ')
    return
end

% generate prepared data structure with some default parameter values
gp=gpdefaults();

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

% main generation loop
for count= 1:gp.runcontrol.num_gen % for 1 : #gera��es
    
    gp.geracoes.num_geracoes=count;

    if count==1;

        % generate the initial population
        gp=initbuild(gp);

        % calculate fitnesses of population members
        gp=evalfitness(gp);
        
        if gp.benchmarks.achou==1
            return
        end
        
%         gp=fitness_correto(gp);

        % update run statistics
%         gp=updatestats(gp)

        % call user defined function (if defined)
%         gp=gp_userfcn(gp)

        % display current stats on screen
%         displaystats(gp)

    else % Crossover, Muta��o e GAP: Variantes
         % Press�o Lexogr�fica: Variantes
        
        num_inicial=6;
        % 1 UP
        if count>=1 && count<=round(0.30*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.05;
            
            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial;
        
        % 2 DOWN    
        elseif count>round(0.30*gp.runcontrol.num_gen) && count<=round(0.40*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.30;
            
            gp.geracoes.pressao.booleano=false;
        
        % 3 UP
        elseif count>round(0.40*gp.runcontrol.num_gen) && count<=round(0.70*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.10;
            
            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial+2;
            
        % 4 DOWN
        elseif count>round(0.70*gp.runcontrol.num_gen) && count<=round(0.80*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.35;
            
            gp.geracoes.pressao.booleano=false;
            
        % 5 UP
        elseif count>round(0.80*gp.runcontrol.num_gen) && count<=round(0.90*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.15;
            
            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial+4;
            
        % 6 DOWN
        elseif count>round(0.90*gp.runcontrol.num_gen) && count<=round(1.00*gp.runcontrol.num_gen)
            count
            gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
            gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
            gp.selection.elite_fraction=0.40;
            
            gp.geracoes.pressao.booleano=true;
            gp.geracoes.pressao.num_classes=num_inicial+6;
            
        end

        % use crossover, mutation etc. to generate a new population
        gp=popbuild(gp);

        % calculate fitnesses of population members
        gp=evalfitness(gp);
        
        if gp.benchmarks.achou==1
            return
        end
        
%         gp=fitness_correto(gp);

        % update run statistics
%         gp=updatestats(gp)

        % call user defined function
%         gp=gp_userfcn(gp)

        % display current stats on screen
%         displaystats(gp)

    end

    %save gp structure % IMPORTANTE! Vem pra c� depois de gp_userfcn
    if  gp.runcontrol.savefreq && ~mod(gp.state.count-1,gp.runcontrol.savefreq)
        save gptips_tmp gp % ONDE EST� ISSO?
    end

    % break out of generation loop if termination required
    if gp.state.terminate
        disp('Fitness criterion met. Terminating run.')
        break
    end

end  %end generation loop

%% TREINO: Gr�ficos

% % Melhor_RMSE: RMSE do Melhor indiv�duo / gera��o
% f=figure('name','RMSE do Melhor indiv�duo / gera��o','numbertitle','off','visible','off');
% plot(min(gp.geracoes.rmse),'k');
% ylabel('RMSE');
% xlabel('Gera��es');
% title(['RMSE do Melhor indiv�duo / gera��o']);
% saveas(f,'Melhor_RMSE.png');

% �ndices dos Melhores indiv�duos (em todas as gera��es)
menores=min(gp.geracoes.rmse);
r2ajustado_menoresRMSE=[];
reg_menoresRMSE=[];
ind_menoresRMSE=[];
for i=1:length(menores)
    ind_menor=find(gp.geracoes.rmse(:,i)==menores(i));
    if length(ind_menor)==1
        r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)];
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
    else
        ind_menor=ind_menor(1);
        r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)];
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
    end
end

gp.geracoes.r2ajustado_menoresRMSE=r2ajustado_menoresRMSE;
gp.geracoes.reg_menoresRMSE=reg_menoresRMSE;
gp.geracoes.ind_menoresRMSE=ind_menoresRMSE;

gp.Bests.R2ajustado=gp.geracoes.r2ajustado_menoresRMSE;
gp.Bests.rmse=menores;
gp.Bests.num_reg=gp.geracoes.reg_menoresRMSE;

% % Melhor_R2ajustado   : R2ajustado do Melhor indiv�duo / gera��o
% f=figure('name','R2ajustado do Melhor indiv�duo / gera��o','numbertitle','off','visible','off');
% plot(r2ajustado_menoresRMSE,'g');
% hold on;
% plot(max(gp.geracoes.r2ajustado),'r');
% ylabel('R2ajustado');
% xlabel('Gera��es');
% title(['R2ajustado do Melhor indiv�duo / gera��o']);
% hold off;
% saveas(f,'Melhor_R2ajustado.png');
% 
% % Melhor_numreg : #reg E #reg_signif do Melhor indiv�duo / gera��o
% f=figure('name','#reg E #reg_signif do Melhor indiv�duo / gera��o','numbertitle','off','visible','off');
% plot(reg_menoresRMSE,'m');
% ylabel('#reg');
% xlabel('Gera��es');
% title(['#reg do Melhor indiv�duo / gera��o']);
% saveas(f,'Melhor_numreg.png');
% 
% % Melhor_prop_numreg_signif : prop_numreg_signif do Melhor indiv�duo / gera��o
% f=figure('name','prop_numreg_signif do Melhor indiv�duo / gera��o','numbertitle','off','visible','off');
% plot(reg_menoresRMSE,'m');
% ylabel('#reg');
% xlabel('Gera��es');
% title(['#reg do Melhor indiv�duo / gera��o']);
% saveas(f,'Melhor_numreg.png');
% 
% 
% 
% 
% 
% 
% 
% % M�dia_RMSE      : M�DIA dos indiv�duos / gera��o
% f=figure('name','RMSE: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','off');
% plot(mean(gp.geracoes.rmse),'r');
% ylabel('RMSE: M�dia dos indiv�duos');
% xlabel('Gera��es');
% title(['RMSE: M�dia dos indiv�duos / gera��o']);
% saveas(f,'M�dia_RMSE.png');
% 
% % M�dia_R2ajustado: M�DIA dos indiv�duos / gera��o
% f=figure('name','R2ajustado: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','off');
% plot(mean(gp.geracoes.r2ajustado),'k');
% ylabel('R2ajustado: M�dia dos indiv�duos');
% xlabel('Gera��es');
% title(['R2ajustado: M�dia dos indiv�duos / gera��o']);
% saveas(f,'M�dia_R2ajustado.png');
% 
% % M�dia_numreg    : M�DIA dos indiv�duos / gera��o
% f=figure('name','#Regressores: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','off');
% plot(mean(gp.geracoes.num_reg),'g');
% ylabel('#Regressores: M�dia dos indiv�duos');
% xlabel('Gera��es');
% title(['#Regressores: M�dia dos indiv�duos / gera��o']);
% saveas(f,'M�dia_numreg.png');

%% TESTE (Dataset)

[gp]=evalfitness_teste(gp);

% gp.geracoes.teste.r2ajustado_teste=r2ajustado_teste;
% gp.geracoes.teste.r2_teste=r2_teste;
% gp.geracoes.teste.rmse_teste=rmse_teste;

%% IGM
% gp.geracoes.testes=gp.geracoes.teste1+gp.geracoes.teste2
% gp.geracoes.IGM=find(gp.geracoes.testes==2);
% gp.geracoes.IGMrmse=gp.geracoes.rmse(gp.geracoes.IGM)
% gp.geracoes.IGMmedia=mean(gp.geracoes.IGMrmse)

% % finalise the run
% gp=gpfinalise(gp);
