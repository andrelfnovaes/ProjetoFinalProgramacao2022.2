function [ind]=selection(gp)
%SELECTION GPTIPS function to probabilistically select an individual from 
%the current population based on fitness. Currently, only tournament selection 
%is implemented.
%   
%   [IND]=SELECTION(GP) returns the index IND in GP.POP of the selected 
%   individual. 
%
%   Remarks:
%   For a tournament of size N:
%
%   1) N individuals are randomly selected from the population WITH
%   RESELECTION ALLOWED.
%   2) The population index of the best individual in the tournament is
%   returned.
%   3) If one or more individuals have the same fitness then one of these
%   is selected randomly unless GP.SELECTION.TOURNAMENT.LEX_PRESSURE is 
%   set to true. IN THIS CASE, the one with the fewest NODES is selected. 
%
%   (if there is more than one individual with the best fitness and the 
%   same number of nodes then one of these individuals is randomly 
%   selected.)
%
%   (c) Dominic Searson 2009
%
%   1.0
%
%   See also: POPBUILD, INITBUILD
 
    % Escolhe índices dos indivíduos que participarão do Torneio
    tour_ind=ceil(rand(gp.selection.tournament.size,1)*gp.runcontrol.pop_size)
       
    %% Pressão Lexográfica 1: #classes fixo
    
%     n_classes=round((max(gp.fitness.values)-min(gp.fitness.values))/(std(gp.fitness.values)/2))
%
%     fitness_real=gp.fitness.values
%     
%     ind_porclasse=3
%     fitness_auxiliar=sort(gp.fitness.values)
%     
%     if mod(gp.runcontrol.pop_size,ind_porclasse)==0
%         n_classes=gp.runcontrol.pop_size/ind_porclasse
%     else
%         n_classes=fix(gp.runcontrol.pop_size/ind_porclasse)+1
%     end
%     
%     indices_fitness_limite=[]
%     for ii=1:(n_classes-1)
%         indices_fitness_limite=[indices_fitness_limite ind_porclasse*ii]
%     end
%     
%     fitness_limite=fitness_auxiliar(indices_fitness_limite)

    %% Pressão Lexográfica 2: Histograma
    
    if gp.geracoes.pressao.booleano==true
   
        fitness_real=gp.fitness.values;

        % quantidade_porbarra=hist(gp.fitness.values,round(gp.runcontrol.pop_size/gp.geracoes.fatordePL));
        quantidade_porbarra=hist(gp.fitness.values,gp.geracoes.pressao.num_classes);
        % hist(gp.fitness.values)

        [fitness_auxiliar,index]=sort(gp.fitness.values,1);

        % [vals, idx] = sort(B,2);
        % A(:,idx)

        soma=0;
        quantidade_porbarra_acum=[];
        for ii=1:length(quantidade_porbarra)
            soma=quantidade_porbarra(ii)+soma;
            quantidade_porbarra_acum=[quantidade_porbarra_acum soma];
        end
        quantidade_porbarra_acum=unique(quantidade_porbarra_acum);

        fitness_limite=fitness_auxiliar(quantidade_porbarra_acum);

        novo_fitness=[];
        fitness_limite2=sort(fitness_limite,'descend');
        for jj=1:length(fitness_limite2)
            replace=find(fitness_auxiliar<=fitness_limite2(jj));
            novo_fitness(replace)=fitness_limite2(jj);
        end
        novo_fitness=novo_fitness';

        final_fitness=[];
        pp=1;
        for kk=1:length(index)
            final_fitness(index(pp),1)=novo_fitness(kk);
            pp=pp+1;
        end

        gp.fitness.values=final_fitness;

    end
            
    %%
        
    % Retorna fitness desses indivíduos
    tour_fitness=gp.fitness.values(tour_ind)
    
    % Max ou Min?
    if gp.fitness.minimisation
        bestfitness=min(tour_fitness) % Diz o menor fitness
    else
        bestfitness=max(tour_fitness)
    end
    
    
    % Matriz só com os índices dos melhores (referência é o vetor do
    % Torneio!)
    bestfitness_tour_ind=find(tour_fitness==bestfitness)
    
    % Número de melhores
    number_of_best=numel(bestfitness_tour_ind) % #de bests
    
        
        %Use plain lexicographic parsimony pressure a la Sean Luke
        %and Livui Panait (Lexicographic Parsimony Pressure, GECCO 2002, pp. 829-836).
        %According to them, this works best when limiting depth of trees as well.
         if number_of_best>1 && gp.selection.tournament.lex_pressure
            
            %each individual may consist of more than one gene so add up
            %the number of nodes of each gene
            lowest_tnodes=inf
            for i=1:number_of_best
               
               tnodes=0
               % Acumula os nós
               for j=1:numel(gp.pop{tour_ind(bestfitness_tour_ind(i))}) % de 1 ao #genes do indiv.
                   tnodes=tnodes+getnumnodes(gp.pop{tour_ind(bestfitness_tour_ind(i))}{j}) % pega o número de cada gene
               end % COMO FUNCIONA o getnumnodes na linha acima?
                   % tnodes: acumula o #nós, gene a gene, por indivíduo
               
               % Entra se o novo_melhor é MENOR(#nós) do que o
               % antigo_melhor
               if tnodes<lowest_tnodes
                 smallest_indiv_ind=bestfitness_tour_ind(i)
                 lowest_tnodes=tnodes
               end
                   
            end
            bestfitness_tour_ind=smallest_indiv_ind
            
        else    %otherwise randomly pick one
        bestfitness_tour_ind=bestfitness_tour_ind(ceil(rand*number_of_best))
        end
    
    
 ind=tour_ind(bestfitness_tour_ind) % ind é a saída = indíce, em pop, do indivíduo selecionado
 
 if gp.geracoes.pressao.booleano==true
    gp.fitness.values=fitness_real
 end
 
