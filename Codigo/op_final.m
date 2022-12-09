function [fitness_popfinal]=op_final(gp)

% DATASETS REAIS: Calcula o fitness da população final, SEM os regressores insignificantes

if gp.benchmarks.opcao==0
    fitness_popfinal=zeros(gp.runcontrol.pop_size,1);
    for i=1:gp.runcontrol.pop_size
        if isempty(gp.geracoes.evalstr_signif{i,gp.runcontrol.num_gen})==0
            [fitness_popfinal(i,1),gp]=feval(gp.fitness.fitfun,gp.geracoes.evalstr_signif{i,gp.runcontrol.num_gen},gp); 
        end
    end

% DATASETS BENCHMARKS: Identifica se achou a regressão, SEM os regressores insignificantes

end

end