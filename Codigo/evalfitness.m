function [gp]=evalfitness(gp)
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

if strcmp(gp.runcontrol.dataset,'treino')==1

    for i=1:gp.runcontrol.pop_size

        gp.geracoes.pop_size=i;

        gp.state.current_individual=i;

        evalstr=tree2evalstr(gp.pop{i},gp);

        gp.fitness.numnodes(i,1)=getnumnodes(gp.pop{i});
        gp.geracoes.numero_nos(i,gp.geracoes.num_geracoes)=gp.fitness.numnodes(i,1);

        temp_evalstr=evalstr; % no fim, evalstr=temp_evalstr
        [~,~,expr_sym,~]=gppretty(gp,i,evalstr); % [expr_sym,cell_expr_sym]

        [eval1,num_reg]=regressores(expr_sym);
        gp.geracoes.num_reg(gp.geracoes.pop_size,gp.geracoes.num_geracoes)=num_reg;
        gp.geracoes.individuos{i,gp.geracoes.num_geracoes}=eval1;

        evalstr=eval1;

        % Trocas: 'x', '^', '*'
        pat='x(\d+)';
        evalstr=regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)');
        pat='^';
        evalstr=strrep(evalstr,pat,'.^');
        pat='*';
        evalstr=strrep(evalstr,pat,'.*');
        % gp.geracoes.evalstr{gp.geracoes.pop_size,gp.geracoes.num_geracoes}=evalstr;

        % Adiocionei (14/05/14)
        % pat='/';
        % evalstr=strrep(evalstr,pat,'./');
        
        [fitness,gp]=feval(gp.fitness.fitfun,evalstr,gp);
        gp.fitness.values(i)=fitness;
        gp.geracoes.fitness(i,gp.geracoes.num_geracoes)=fitness;

        %% Datasets Benchmarks (não uso mais)
        
        if gp.benchmarks.opcao==1

            if gp.benchmarks.estat_signif==1
                individuo=eval1(gp.geracoes.index_reg_signif{gp.geracoes.pop_size,gp.geracoes.num_geracoes});
                betas_temp=gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes};
                beta_0=betas_temp(1,1);
                betas_temp=betas_temp(2:end,1);
                betas_signif=betas_temp(gp.geracoes.index_reg_signif{gp.geracoes.pop_size,gp.geracoes.num_geracoes});
            else
                individuo=eval1;
            end

            matriz=[];
            for ii=1:length(gp.benchmarks.elementos)
                matriz=[matriz;strcmp(gp.benchmarks.elementos{1,ii},individuo)];
                if isempty(individuo)==1
                   break
                end
                if sum(matriz(ii,:))==0
                    break
                end
            end

            tamanho_matriz=size(matriz);
            if (tamanho_matriz(1)==tamanho_matriz(2))&&((sum(sum(matriz)))==length(gp.benchmarks.elementos))

                if gp.benchmarks.estat_signif==1

                    indiv=individuo';
                    elem=gp.benchmarks.elementos'; % vertical
                    coef=gp.benchmarks.coef; % vertical

                    [~,elem_ordem]=sort(elem); % '~' = elem_sort
                    coef_ordem=coef(elem_ordem);

                    % Cria limites p/ significantes
                    coef_ordem_UP=coef_ordem+0.05;
                    coef_ordem_DOWN=coef_ordem-0.05;

                    [~,indiv_ordem]=sort(indiv); % '~' = indiv_sort
                    betas_signif_ordem=betas_signif(indiv_ordem);

                    % Teste p/ significantes
                    down_signif=(coef_ordem_DOWN<=betas_signif_ordem);
                    up_signif=(betas_signif_ordem<=coef_ordem_UP);
                    if down_signif==up_signif
                       teste_signif=true;
                    end

                    % Cria limites p/ INsignificantes
                    betas_INsignif=betas_temp(gp.geracoes.index_reg_INsignif{gp.geracoes.pop_size,gp.geracoes.num_geracoes});
                    [linhas_beta,~]=size(betas_INsignif);
                    down_2=zeros(linhas_beta,1)-0.01;
                    up_2=zeros(linhas_beta,1)+0.01;
                    if isempty(betas_INsignif)
                        betas_INsignif=1;
                        down_2=1;
                        up_2=1;
                    end

                    % Teste p/ INsignificantes
                    down_INsignif=(down_2<=betas_INsignif);
                    up_INsignif=(betas_INsignif<=up_2);
                    if down_INsignif==up_INsignif
                       teste_INsignif=true;
                    end

                    if (teste_signif==1) && (teste_INsignif==1)
                        gp.benchmarks.achou=1;
                        gp.benchmarks.results.individuo=eval1';
                        gp.benchmarks.results.beta=gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes};                                     
                        return
                    end

                    return

                else
                    gp.benchmarks.achou=1;
                    ind=eval1;
                    ind=ind';
                    gp.benchmarks.results.individuo=ind;
                    gp.benchmarks.results.beta=gp.geracoes.betas{gp.geracoes.pop_size,gp.geracoes.num_geracoes};
                    return
                end

                % return

            end

        end

    end

elseif strcmp(gp.runcontrol.dataset,'teste')==1
    
    % Indivíduo com menor fitness
    i=gp.geracoes.ind_menoresRMSE(end);
    
    gp.state.current_individual=i;
    evalstr=tree2evalstr(gp.pop{i},gp);
    temp_evalstr=evalstr;
    [a,b,expr_sym,c]=gppretty(gp,i,evalstr);
    [eval1,num_reg]=regressores(expr_sym);
    evalstr=eval1;
    
    pat='x(\d+)';
    evalstr=regexprep(evalstr,pat,'gp.userdata.xtest(:,$1)');
    pat='^';
    evalstr=strrep(evalstr,pat,'.^');
    pat='*';
    evalstr=strrep(evalstr,pat,'.*');
    
    % [gp]=regressmulti_fitfun(evalstr,gp);
    [~,gp]=feval(gp.fitness.fitfun,evalstr,gp);
    
end

