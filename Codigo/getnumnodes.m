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

if isa(expr,'char') % Se expr � tipo char, true.

    numnodes=getnn(expr); % fun��o l� embaixo
    return

elseif iscell(expr) % vindo de evalfitness, veio por aqui
                    % gp.pop{i} n�o � char, 2o GPTIPS.
    numexpr=length(expr); % numexpr = #genes do indiv�duo

    if numexpr<1
        error('Cell array must contain at least one valid symbolic expression')
    else
        numnodes=0;
        for i=1:numexpr % numexpr = #genes do indiv�duo
            numnodes=numnodes+getnn(expr{i});
        end
    end
else
    error('Illegal argument')
end

    function numnodes=getnn(expr)
    % Sub-function to get numnodes from a single symbolic string

    %NUMBER OF NODES = number of open brackets + number of inputs + number of
    %constants ; #n�s = #par�nteses_esquerdos + #inputs + #ctes
    open_br=findstr(expr,'(');
    open_sq_br=findstr(expr,'[');
    inps=findstr(expr,'x');

    num_open=numel(open_br);
    num_const=numel(open_sq_br);
    num_inps=numel(inps);

    numnodes=num_open+num_const+num_inps;
    % #n�s = #par�nteses_esquerdos + #inputs + #ctes