function [indice_vencedor]=selection(gp)

% 'Escolhe Índices' dos indivíduos que participarão do Torneio
% Obs: pode haver repetição de indivíduos
indices_torneio=ceil(rand(gp.selection.tournament.size,1)*gp.runcontrol.pop_size);

%% Pressão Lexográfica (Histograma): Entrada
    
if gp.geracoes.pressao.booleano==true

    fitness_real=gp.fitness.values;

    % quantidade_porbarra=hist(gp.fitness.values,round(gp.runcontrol.pop_size/gp.geracoes.fatordePL));
    quantidade_porbarra=hist(gp.fitness.values,gp.geracoes.pressao.num_classes);
    % hist(gp.fitness.values)

    [fitness_auxiliar,index]=sort(gp.fitness.values,1);

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

% 'Retorna Fitness' dos indivíduos que participarão do Torneio
fitness_torneio=gp.fitness.values(indices_torneio);

% Max ou Min?
if gp.fitness.minimisation
    bestfitness=min(fitness_torneio);
else
    bestfitness=max(fitness_torneio);
end

posicoes_temp=find(bestfitness==fitness_torneio);
posicoes_torneio=indices_torneio(posicoes_temp);

% % NÚMERO de regressores significantes (quantidade)
numreg_signif_torneio=gp.geracoes.numreg_signif(posicoes_torneio);
[~,indice_temp]=max(numreg_signif_torneio);
indice_vencedor=posicoes_torneio(indice_temp);

% % PROPORÇÃO de regressores significantes
% prop_numreg_signif_torneio=gp.geracoes.prop_numreg_signif(posicoes_torneio);
% [~,indice_temp]=max(prop_numreg_signif_torneio);
% indice_vencedor=posicoes_torneio(indice_temp);

% % MÉDIA de Significância
% media_signif_torneio=gp.geracoes.media_signif(posicoes_torneio);
% [~,indice_temp]=max(media_signif_torneio);
% indice_vencedor=posicoes_torneio(indice_temp);

%% Pressão Lexográfica (Histograma): Saída
if gp.geracoes.pressao.booleano==true
    gp.fitness.values=fitness_real;
end

end