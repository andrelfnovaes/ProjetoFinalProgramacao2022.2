function [gp]=initbuild(gp)
%INITBUILD GPTIPS function to generate an initial population of GP 
%individuals.
%
%   [GP]=INITBUILD(GP) creates an initial population using the parameters 
%   in the structure GP. Various changes to fields of GP are made.
%
%   Remarks:
%   Each individual in the population can contain 1 or more separate trees
%   (if more than 1, then each tree is referred to as a gene of the 
%   individual). This is set by the user in the config file in the field 
%   GP.GENES.MULTIGENE.
%
%   Each individual is a cell array. You can use cell addressing to 
%   retrieve the individual genes. E.g. GP.POP{3}{2} will return the 2nd 
%   gene in the third individual.
%
%   Trees are (by default) constructed using a probabilistic version
%   of Koza's ramped 1/2 and 1/2 method.
%   E.g. if maximum tree depth is 5 and population is 100
%   then, on average, 20 will be generated at each depth (1/2 using
%   'grow' and 1/2 using 'full').
%
%   Multiple copies of genes in an individual are disallowed when
%   individuals are created. There is, however, no such restriction 
%   imposed on future generations.
%
%   (c) Dominic Searson 2009
%
%   v1.0
%
%   See also POPBUILD, TREEGEN




% Extract temp variables from the gp structure
pop_size=gp.runcontrol.pop_size;
max_nodes=gp.treedef.max_nodes;
init_genes=gp.genes.max_genes;

%override any gene settings if using single gene gp
if ~gp.genes.multigene
    init_genes=1;
end


%initialise vars ; "GP.POP{3}{2} will return the 2nd gene in the third individual."
gp.pop=cell(pop_size,1); % estrutura com o tamanho da população
num_genes=1;

%initialise gen build counter
gp.state.count=1;

%building process
for i=1:pop_size %loop through population
    
    %randomly pick num of genes in individual
    if init_genes>1
        num_genes=ceil(rand*init_genes); % num_genes é RAND!
    end

    individ=cell(1,(num_genes)); % construct empty individual

    for z=1:num_genes %loop through genes in each individual and generate a tree for each


        %Generate gene z and check that genes 1..z-1 are different
        loop_counter=0;
        while 1
            loop_counter=loop_counter+1;

            %generate a trial tree for gene z
            temp=treegen(gp);
            %temp = ['x' sprintf('%d',ceil(rand*gp.nodes.inputs.num_inp))]
            numnodes=getnumnodes(temp); % retorna #nós de temp=tree_string

            if numnodes<=max_nodes 

                copy_detected=false;

                if z>1 %check previous genes for copies

                    for j=1:z-1
                        if strcmp(temp,individ{1,j})
                            copy_detected=true; % = 1
                            break
                        end
                    end

                end

                if ~copy_detected % ~copy_detected=0 -> retorna 1
                    break
                end

            end %max nodes check

            % display a warning if having difficulty building trees due to
            % user constraints
            if loop_counter>3
                disp('initbuild: iterating tree build loop because of constraints. Consider setting gp.treedef.max_nodes higher')
            end

        end %while loop

        % ok so far. Assign tree as gene in current individual
        individ{1,z}=temp; % ATENÇÃO em individ, criado por 'cell'

    end %gene loop

    %all genes done. write new individual to population cell array
    gp.pop{i,1}=individ; % ATENÇÃO em gp.pop, criado por 'cell'
end



