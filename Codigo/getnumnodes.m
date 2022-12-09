function numnodes=getnumnodes(expr)
%GETNUMNODES GPTIPS function to return the number of nodes in a GP 
%expression or the total node count for a cell array of expressions.
%
%   NUMNODES=GETNUMNODES(EXPR) finds the number of nodes in the GP
%   expression EXPR or cell array of expressions EXPR.
%
%   (c) Dominic Searson 2009
%
%    v1.0
%
%    See also GETDEPTH

if isa(expr,'char') % Se expr é tipo char, true.

    numnodes=getnn(expr); % função lá embaixo
    return

elseif iscell(expr) % vindo de evalfitness, veio por aqui
                    % gp.pop{i} não é char, 2o GPTIPS.
    numexpr=length(expr); % numexpr = #genes do indivíduo

    if numexpr<1
        error('Cell array must contain at least one valid symbolic expression')
    else
        numnodes=0;
        for i=1:numexpr % numexpr = #genes do indivíduo
            numnodes=numnodes+getnn(expr{i});
        end
    end
else
    error('Illegal argument')
end

    function numnodes=getnn(expr)
    % Sub-function to get numnodes from a single symbolic string

    %NUMBER OF NODES = number of open brackets + number of inputs + number of
    %constants ; #nós = #parênteses_esquerdos + #inputs + #ctes
    open_br=findstr(expr,'(');
    open_sq_br=findstr(expr,'[');
    inps=findstr(expr,'x');

    num_open=numel(open_br);
    num_const=numel(open_sq_br);
    num_inps=numel(inps);

    numnodes=num_open+num_const+num_inps;
    % #nós = #parênteses_esquerdos + #inputs + #ctes