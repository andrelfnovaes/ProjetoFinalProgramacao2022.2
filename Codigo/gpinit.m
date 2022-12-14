function gp=gpinit(gp)
%GPINIT GPTIPS function to initialise the GP run.
%
%   GP=GPINIT(GP)  
%
%   (c) Dominic Searson 2009
%
%   v1.0
%
%   See also GPCHECK

%call subfunction to process function nodes before run
gp=procfuncnodes(gp);

% Throw an error if there are no inputs, p_ERC=0 AND there are no arity zero
% functions active... ESTRANHO!
if  gp.nodes.inputs.num_inp==0 && gp.nodes.const.p_ERC==0 && ...
        isempty(find(gp.nodes.functions.arity(logical(gp.nodes.functions.active))==0)) 
    error('No terminals (inputs, constants or zero arity functions) have been defined for this run!') 
end



%next, initialise some state and tracker variables
gp.state.count=[];
gp.state.best.fitness=[];
gp.state.best.individual=[];
gp.state.run_completed=false;
gp.state.current_individual=[];
gp.state.std_devfitness=[];
gp.state.terminate=false;
gp.fitness.returnvalues=cell(gp.runcontrol.pop_size,1);

gp.geracoes.individuos=cell(1);
gp.geracoes.fitness=zeros(gp.runcontrol.pop_size,1);
gp.geracoes.betas=cell(gp.runcontrol.pop_size,1);
gp.geracoes.numero_nos=zeros(gp.runcontrol.pop_size,1);
gp.geracoes.teste1=zeros(gp.runcontrol.pop_size,1);
gp.geracoes.teste2=zeros(gp.runcontrol.pop_size,1);
gp.geracoes.teste3=zeros(gp.runcontrol.pop_size,1);

%process mutation probabilities vector
gp.operators.mutation.cumsum_mutate_par=cumsum(gp.operators.mutation.mutate_par);


%init. history variables
gp.results.history.bestfitness=zeros(gp.runcontrol.num_gen,1);
gp.results.history.meanfitness=zeros(gp.runcontrol.num_gen,1);
gp.results.history.std_devfitness=zeros(gp.runcontrol.num_gen,1);


%best of run fields
gp.results.best.fitness=[];
gp.results.best.individual=[];
gp.results.best.returnvalues=[];
gp.results.best.foundatgen=[];


%assign field holding fitnesses to gp structure
gp.fitness.values=zeros(gp.runcontrol.pop_size,1);

gp.fitness.numnodes=zeros(gp.runcontrol.pop_size,1);

if ~gp.runcontrol.quiet % "~vetor" = se vetor(i)=0, retorna 1.
    %disp run info
    fns=[];
    for i=1:length(gp.nodes.functions.active_name_UC)
        fns=[fns ' ' gp.nodes.functions.active_name_UC{i}];
    end

    
    if gp.selection.tournament.lex_pressure
        lex_inf='True';
    else
        lex_inf='False';
    end
    
    disp(' ') % pula 1 linha
    disp(' ')
    
    disp('-------------------------------------------------------------------------')
    disp('GPTIPS: Genetic programming toolbox')
    disp('Copyright (C) 2010 Dominic Searson')
    disp(' ')
    disp('This program is free software: you can redistribute it and/or modify')
    disp('it under the terms of the GNU General Public License as published by')
    disp('the Free Software Foundation, either version 3 of the License, or')
    disp('(at your option) any later version.')
    disp(' ')
    disp('This program is distributed in the hope that it will be useful,')
    disp('but WITHOUT ANY WARRANTY; without even the implied warranty of')
    disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the')
    disp('GNU General Public License for more details: http://www.gnu.org/licenses')
    disp('-------------------------------------------------------------------------- ')
    disp(' ')
    disp('Run parameters')
    disp('--------------')
    disp(['Population size:         ' int2str(gp.runcontrol.pop_size)])
    disp(['Number of generations:   ' int2str(gp.runcontrol.num_gen)])
    disp(['Tournament size:         ' int2str(gp.selection.tournament.size)])
    disp(['Lexicographic selection: ' lex_inf])
    disp(['Max tree depth:          ' int2str(gp.treedef.max_depth)])
    disp(['Max nodes per tree:      ' int2str(gp.treedef.max_nodes)])
    disp(['Using function set:     ' fns])
    disp(['Number of inputs:        ' int2str(gp.nodes.inputs.num_inp)])

    if gp.genes.multigene
    disp(['Max genes:               ' int2str(gp.genes.max_genes)])
    end

    if ~gp.nodes.const.p_ERC % IMPORTANTE: "gp.nodes.const.p_ERC" responde pelo uso de ctes!
        disp('Using no constants')
    else
    disp(['Constants range:         [' num2str(gp.nodes.const.range) ']'])
    end

    disp(['Using fitness function:  ' func2str(gp.fitness.fitfun) '.m'])
    disp(' ')
end

%log run start time
gp.info.start_time=datestr(now,0);





function gp=procfuncnodes(gp)
% PROCFUNCNODES GPTIPS function to to process required function node information prior to a run.



%loop through function nodes and generate arity list
for i=1:length(gp.nodes.functions.name)

    arity=nargin(gp.nodes.functions.name{i}); % #par?metros de cada opera??o. Ex: arity(times)=2
    
    
    
    %some functions have a variable number of input arguments (e.g. rand)
    %In this case generate an error message and exit
    if arity==-1
        error(['The function ' gp.nodes.functions.name{i} ' may not be used as a tree node because it has a variable number of arguments.'])
    end
    
    gp.nodes.functions.arity(i)=arity;
    
end

if ~isfield(gp.nodes.functions,'active') || isempty(gp.nodes.functions.active)
   gp.nodes.functions.active=ones(1,length(gp.nodes.functions.name));
end


gp.nodes.functions.active=logical(gp.nodes.functions.active);

%check max number of allowed functions not exceeded
gp.nodes.functions.num_active=numel(find(gp.nodes.functions.active));
if gp.nodes.functions.num_active>22 % IMPORTANTE: # m?ximo de opera??es ? 22.
    error('Maximum number of active functions allowed is 22')
end


% Generate single char Active Function Identifiers (AFId)(a->z excluding x,e,i,j) to stand in for function names
% whilst processing expressions. Exclusions are because x is reserved for input nodes, e is used
% for expressing numbers in standard form by Matlab and, by default, i and
% j represent sqrt(-1) in Matlab.
charnum=96; skip=0;
for i=1:gp.nodes.functions.num_active
    while true      %e                          %i                    %j                         %x
        if (charnum+i+skip)==101 || (charnum+i+skip)==105 || (charnum+i+skip)==106 || (charnum+i+skip)==120
            skip=skip+1;
        else
            break
        end

    end
    afid(i)=char(charnum+i+skip); % afid = abcdf, por exemplo (? uma string corrida)
end


% Extract upper case active function names for later use
gp.nodes.functions.afid=afid;
temp=cell(gp.nodes.functions.num_active,1); % temp = [ [] [] ... [] ]length(afid)X1
[temp{:}]=deal(gp.nodes.functions.name{gp.nodes.functions.active}); % temp recebe ['times' ... ]
[gp.nodes.functions.active_name_UC]=upper(temp); % upper -> letras Mai?sculas




%Generate index locators for arity>0 and arity=0 active functions (The treegen function needs
%this info later for identifying which functions are terminal and which are internal nodes)
active_ar=(gp.nodes.functions.arity(gp.nodes.functions.active)); % = s? arity's ativados
fun_argt0=active_ar>0; % mesmo array anterior: repleto de 1s
fun_areq0=~fun_argt0;  % mesmo array anterior: repleto de 0s

gp.nodes.functions.afid_argt0=gp.nodes.functions.afid(fun_argt0); %functions with arity>0
gp.nodes.functions.afid_areq0=gp.nodes.functions.afid(fun_areq0); %functions with arity=0
gp.nodes.functions.arity_argt0=active_ar(fun_argt0);