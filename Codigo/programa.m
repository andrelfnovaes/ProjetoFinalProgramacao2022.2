for i=1:5
    
    gpdemo4
    
    Bests.R2ajustado(i,:)=gp.geracoes.r2ajustado_menoresRMSE;
    Bests.rmse(i,:)=min(gp.geracoes.rmse);
    Bests.num_reg(i,:)=gp.geracoes.reg_menoresRMSE;
    
    filename=['gp' num2str(i) '.mat' ];
    
    save(filename);
    
end

save Bests

%% MM: R2ajustado, RMSE, #Regressores

% MM_R2ajustado: M�DIA dos Melhores / gera��o
f=figure('name','R2ajustado: M�dia dos Melhores / gera��o','numbertitle','off','visible','off');
plot(mean(Bests.R2ajustado),'k');
ylabel('R2ajustado: M�dia dos Melhores');
xlabel('Gera��es');
title(['R2ajustado: M�dia dos Melhores / gera��o']);
saveas(f,'MM_R2ajustado.png');

% MM_rmse: M�DIA dos Melhores / gera��o
f=figure('name','RMSE: M�dia dos Melhores / gera��o','numbertitle','off','visible','off');
plot(mean(Bests.rmse),'r');
ylabel('RMSE: M�dia dos Melhores');
xlabel('Gera��es');
title(['RMSE: M�dia dos Melhores / gera��o']);
saveas(f,'MM_rmse.png');

% MM_numreg: M�DIA dos Melhores / gera��o
f=figure('name','#Regressores: M�dia dos Melhores / gera��o','numbertitle','off','visible','off');
plot(mean(Bests.num_reg),'g');
ylabel('#Regressores: M�dia dos Melhores');
xlabel('Gera��es');
title(['#Regressores: M�dia dos Melhores / gera��o']);
saveas(f,'MM_numreg.png');

