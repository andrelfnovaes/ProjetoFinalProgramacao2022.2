function [gp]=gpdemo1_config(gp);
%GPDEMO1_CONFIG example configuration file demonstrating simple symbolic 
%regression.
%
%   The simple quartic polynomial (y=x+x^2+x^3+x^4) from John Koza's 1992 
%   Genetic Programming book is used.
%
%   [GP]=GPDEMO1_CONFIG(GP) returns the user specified parameter structure 
%   GP for the quartic polynomial problem.
%   
%   Example:
%   GP=GPTIPS('GPDEMO1_CONFIG') performs a GPTIPS run using this 
%   configuration file and returns the results in a structure called GP.
%
%   v1.0
%
%   (C) Dominic Searson 2009
% 
%   See also QUARTIC_FITFUN, GPDEMO1



% Main run control parameters
% ---------------------------
gp.runcontrol.pop_size=50			% Population size
gp.runcontrol.num_gen=100			% Number of generations to run for including generation zero
                                    % (i.e. if set to 100, it'll finish after generation 99).
gp.runcontrol.verbose=5          % Set to n to display run information to screen every n generations




% Selection method options
% -------------------------

gp.selection.tournament.size=2       
gp.selection.elite_fraction=0.02     



% Fitness function and optimisation specification 
% ------------------------------------------------
gp.fitness.fitfun=@quartic_fitfun           % Function handle of the user's fitness function (filename with no .m extension).
gp.fitness.minimisation=true                % True if to minimise the fitness function (if false it is maximised).
gp.fitness.terminate=true                   % True to terminate run early if fitness threshold met
gp.fitness.terminate_value=1e-3

% User data specification (sets up quartic polynomial data)  
% ----------------------------------------------------------------
x=linspace(-1,1,20)' %generate 20 data points in the range [-1 1]
gp.userdata.x=x
gp.userdata.y=x+x.^2+x.^3+x.^4 %generate y



% Input configuration
% --------------------
% This sets the number of inputs 
gp.nodes.inputs.num_inp=1 		         



% Constants
% ---------

% When building a tree this is
% the probability with which a constant will
% be selected instead of a terminal.
% [1=all ERCs, 0.5=1/2 ERCs 1/2 inputs, 0=no ERCs]
gp.nodes.const.p_ERC=0		   % quartic example doesn't need constants            



% Tree build options
% -------------------
             
% Maximum depth of trees 
gp.treedef.max_depth=12
 	              
% Maximum depth of sub-trees created by mutation operator
gp.treedef.max_mutate_depth=7



% Define function nodes
% ---------------------

gp.nodes.functions.name{1}='times'
gp.nodes.functions.name{2}='minus'
gp.nodes.functions.name{3}='plus'
gp.nodes.functions.name{4}='rdivide'



