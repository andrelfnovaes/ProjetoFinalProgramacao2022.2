function [gp]=evalfitness_teste(gp)
%EVALFITNESS GPTIPS function to call the user specified fitness function.
%
%    [GP]=EVALFITNESS(GP) evaluates the fitnesses of individuals stored 
%    in the GP structure and updates various other fields of GP accordingly.
%
%    (c) Dominic Searson 2009
%
%    v1.0
%
%    See also TREE2EVALSTR

% Loop through population and calculate fitnesses

% for i=1:gp.runcontrol.pop_size
    
%     gp.geracoes.pop_size=i;

    % update state to reflect the index of the individual that is about to
    % be evaluated
    
    %% Indivíduo com menor Fitness
    i=gp.geracoes.ind_menoresRMSE(end);
    
    %%
    gp.state.current_individual=i;
   
    %First preprocess the cell array of string expressions into a form that
    %Matlab can evaluate, EVALSTR = individ pra cálculo de fitness!
    evalstr=tree2evalstr(gp.pop{i},gp);
    %evalstr = gp.pop{i}

    %store number of nodes (sum total for all genes)
%     gp.fitness.numnodes(i,1)=getnumnodes(gp.pop{i});
%     gp.geracoes.numero_nos(i,gp.geracoes.num_geracoes)=gp.fitness.numnodes(i,1);



    % Evaluate gp individual using fitness function
    % (the TRY CATCH is to assign a poor fitness value
    % to trees that violate Matlab's
    % daft 'Nesting of {, [, and ( cannot exceed a depth of 32.' error.
% try

    
    temp_evalstr=evalstr; % no fim, evalstr=temp_evalstr
%     save temp_evalstr
    [a,b,expr_sym,c]=gppretty(gp,i,evalstr); % [expr_sym,cell_expr_sym]
%     gp.geracoes.individuos{i,gp.geracoes.num_geracoes}=expr_sym;
    
    [eval1,num_reg]=regressores(expr_sym);
%     gp.geracoes.num_reg(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=num_reg;
    
%     save eval1
    evalstr=eval1;
    
    % Trocas: 'x', '^', '*'
    pat='x(\d+)';
    evalstr=regexprep(evalstr,pat,'gp.userdata.xtest(:,$1)');
    pat='^';
    evalstr=strrep(evalstr,pat,'.^');
    pat='*';
    evalstr=strrep(evalstr,pat,'.*');
    
%     [gp,r2ajustado_teste,r2_teste,rmse_teste]=feval(gp.fitness.fitfun_teste,evalstr,gp);
    [gp]=regressmulti_fitfun_teste(evalstr,gp);
    
%     gp.fitness.values(i)=fitness;
%     gp.geracoes.fitness(i,gp.geracoes.num_geracoes)=fitness;

% catch
%     if ~strncmpi(lasterr,'Nesting of {',12)
%         error(lasterr)
%     end
%     if gp.fitness.minimisation
%         gp.fitness.values(i)=Inf
%     else
%         gp.fitness.values(i)=-Inf
%     end
% end

%     evalstr=temp_evalstr

% end



