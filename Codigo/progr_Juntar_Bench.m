num_rodadas=50;

for i=1:num_rodadas
    
    load(['gp' num2str(i) '.mat']) 
    
    Bench.tempos(i,:)=gp.geracoes.tempo
    Bench.acertou(i,:)=gp.benchmarks.achou
    Bench.geracao(i,:)=gp.geracoes.num_geracoes
    
end

% Tempo
Bench.tempo.tempo_medio=mean(Bench.tempos)
Bench.tempo.tempo_total=sum(Bench.tempos)

% Geração
Bench.taxa_acerto=sum(Bench.acertou)./num_rodadas

% Qdo acerta
indices_acertos=find(Bench.acertou==1)
Bench.qdo_acerta.rodadas_acertivas=indices_acertos

Bench.qdo_acerta.tempo=Bench.tempos(indices_acertos)
Bench.qdo_acerta.tempo_medio=mean(Bench.qdo_acerta.tempo)
Bench.qdo_acerta.tempo_total=sum(Bench.qdo_acerta.tempo)

Bench.qdo_acerta.numger=Bench.geracao(indices_acertos)
Bench.qdo_acerta.numger_medio=mean(Bench.qdo_acerta.numger)

save Bench