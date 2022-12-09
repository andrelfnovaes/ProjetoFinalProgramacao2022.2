function tree_string=treegen(gp,fixed_depth)
%TREEGEN GPTIPS function to generate a new tree expression in prefix 
%notation.
%
%   [TREE_STRING]=TREEGEN(GP) generates a symbolic expression TREE_STRING
%   using the parameters stored in the structure GP.
%
%   [TREE_STRING]=TREEGEN(GP,DEPTH) generates an expression TREE_STRING of 
%   depth DEPTH but otherwise using the parameters stored in the GP 
%   structure.
%
%   IMPORTANTE: Trees are nominally built to a depth specified in the field 
%   GP.TREEDEF.MAX_DEPTH according to the method specified in the field 
%   GP.TREEDEF.BUILD_METHOD. However, if the optional additional input 
%   argument DEPTH is specified (which must be an integer > 0) then a full 
%   tree of this depth will be built irrespective of the settings in the GP
%   structure.
%
%   Remarks:
%   In GPTIPS, a tree of depth 1 contains a single node. 
%
%   (c) Dominic Searson 2009
%
%   v1.0
%
%   See also EXTRACT


%Check input arguments for depth overide parameter
if nargin<2
    fixed_depth=0;
end


%Extract parameters from gp structure
max_depth=gp.treedef.max_depth  ;       % Will use user max depth unless overridden by input argument
build_method=gp.treedef.build_method;   % 1=full 2=grow 3=ramped 1/2 and 1/2
p_ERC=gp.nodes.const.p_ERC;             % Terminal config. [0=no constants, 0.5=half constants half inputs 1=no inputs]

num_inp=gp.nodes.inputs.num_inp;
range=gp.nodes.const.range;
rangesize=range(2)-range(1);

use_mformat=gp.nodes.const.use_matlab_format; %Set to 1 if using matlab defaults for formatting numeric strings (slower, but may be necessary for some problems)
fi =gp.nodes.const.format;


afid_argt0=gp.nodes.functions.afid_argt0; %functions with arity>0
afid_areq0=gp.nodes.functions.afid_areq0; %functions with arity=0;


arity_argt0=gp.nodes.functions.arity_argt0; % args(a)=2, args(b)=2, ...


fun_lengthargt0=numel(afid_argt0); % numel('abcdf') = 5
fun_lengthareq0=numel(afid_areq0); % numel(") = numel([]) = 0



%if a fixed depth tree is required, use 'full' build method
if fixed_depth % fixed_depth = 0 ---> não pede uma tree de fixed depth
    max_depth=fixed_depth;
    build_method=1;
end


% If using ramped 1/2 and 1/2 then RANDOMLY choose max_depth and
% build_method
if build_method==3
    max_depth=ceil(rand*gp.treedef.max_depth);
    build_method=ceil(rand*2);  %set either to 'full' or 'grow' for duration of function
end



%Initial string structure (nodes/subtrees to be built are marked with the $ character)
tree_string='($)';


%recurse through branches = recursão através dos ramos
while true % ATENÇÃO! ATENÇÃO! ATENÇÃO! : pulei vários blocos aqui (cada caso é um caso)

    %Find leftmost node token $ and store string position
    node_token_ind=findstr('$',tree_string); % findstr('$','($)') = 2, pois tem o "("
    
    %breakout if no-more nodes to fill; AQUI É QUE TERMINA!
    if isempty(node_token_ind)
        break
    end

    node_token=node_token_ind(1); % node_token = 1a posição do '$' (ou do '@')

    %REPLACE this node token with 'ACTIVE' node token, @ : ($) ---> (@)
    tree_string(node_token)='@';

    %Count brackets from beginning of string to @, to work out depth of @
    left_seg=tree_string(1:node_token);

    %count number of open brackets
    num_open=numel(findstr(left_seg,'(')); % ordem?...Ok! ; numel = #elementos (# de '(') no array
    num_closed=numel(findstr(left_seg,')'));
    depth=num_open-num_closed;




    % Choose type of node to insert based on depth and building method.
    % If root node, then pick a non-terminal node (unless otherwise
    % specified).
    if depth==1
        node_type=1; % NÓS: 1 = INTERNAL, 2 = TERMINAL, 3 = EITHER
            if max_depth==1 %if max_depth is 1 must always pick a terminal
            node_type=2;
        end

        %if less than max_depth and
        %method is 'full' then only pick a function with arity>0
    elseif build_method==1 && depth<max_depth
        node_type=1;
    elseif depth==max_depth % if depth is max_depth, then just pick terminals
        node_type=2;
    else %pick either with equal prob.
        node_type=ceil(rand*2);
    end



    if node_type==1 % i.e a function with arity>0 (IMPORTANTE)
        fun_choice=ceil(rand*fun_lengthargt0); % fun_lengthargt0 NÃO MUDA! É RAND AQUI!
        fun_name=afid_argt0(fun_choice); % EX: fun_choice = 4 ---> fun_name = 'd' (de 'abcdf')
        num_fun_arg=arity_argt0(fun_choice); % EX: " ---> num_fun_arg = 1 (de [2 2 2 1 1])

        %create appropriate replacement string e.g. a($,$) for 2 argument function
        fun_rep_str=[fun_name '($']; % cria estrutura
        if num_fun_arg>1
            for j=2:num_fun_arg
                fun_rep_str=[fun_rep_str ',$'];
            end
            fun_rep_str=[fun_rep_str ')'];
        else % for single argument functions
            fun_rep_str=[fun_name '($)']; % diferente da linha 134 (fechou o caixão aqui)
        end

        %now replace active node token @ with replacement string
        tree_string=strrep(tree_string,'@',fun_rep_str);

    elseif node_type==2 %pick a terminal (or an arity 0 function, if active)

        %choose a function or input with equal probability
        if fun_lengthareq0 && (rand>=0.5 || (num_inp==0 && p_ERC==0))
            fun_choice=ceil(rand*fun_lengthareq0);
            fun_name=afid_areq0(fun_choice);


            %create appropriate replacement string for 0 arity function
            fun_rep_str=[fun_name '()'];

            % CRIA A TREE!!! (now replace active node token @ with replacement string
            tree_string=strrep(tree_string,'@',fun_rep_str);

        else %choose an input (if it exists) or
            %check if ERCs are switched on; if so, pick a ERC node/arithmetic tree token(ï¿½) with p_ERC prob.
            r_num=rand;
            if r_num>=p_ERC
                r_num=1; %then use an ORDINARY TERMINAL, not a constant
            else
                r_num=2;
            end


            if r_num==1 %then use an ordinary terminal not a constant
                inp_choice=ceil(rand*num_inp);

                tree_string=strrep(tree_string,'@',['x' sprintf('%d',inp_choice)]); % resultado: @ ---> '(x139)'

            elseif r_num==2 %use an ERC token (? character)
                tree_string=strrep(tree_string,'@','?');
            end
        end

        if depth>max_depth
            error('Tree depth exceeded')
        end
    end

    %end of while loop
end




% constant processing
if p_ERC

    % insert ERCs for ? tokens
    c_ind=findstr(tree_string,'?');

    numc=numel(c_ind); % Se c_ind = [], então numel(c_ind) = 0

    %Generate reqd. number of constants in the allowed range: GEROU CTE!
    const=rand(1,numc)*rangesize + range(1);

    for k=1:numc


        % 3) Convert the constant to square bracketed string form:
        if use_mformat
            arg=['[' num2str(const(k)) ']'];
        else

            if const(k)==fix(const(k)) % EX: const = 2.1412 ; fix(2.1412) = 2
                arg=['[' sprintf('%.0f',const(k)) ']'];
            else
                arg=['[' sprintf(fi,const(k)) ']']; % retorna: [2.141225] ; Obs: 6 casas decimais
            end
        end
        %replace ? with const
        [main_tree]=extract(c_ind(1),tree_string); % main_tree ~= tree_string ('?'->'$')
        tree_string=strrep(main_tree,'$',arg); % troca '$' pela cte
        c_ind=findstr(tree_string,'?'); % c_ind = []
    end

end

%strip off outside brackets
tree_string=tree_string(2:end-1);
