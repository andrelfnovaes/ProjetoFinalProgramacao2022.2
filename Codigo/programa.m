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

% MM_R2ajustado: MÉDIA dos Melhores / geração
f=figure('name','R2ajustado: Média dos Melhores / geração','numbertitle','off','visible','off');
plot(mean(Bests.R2ajustado),'k');
ylabel('R2ajustado: Média dos Melhores');
xlabel('Gerações');
title(['R2ajustado: Média dos Melhores / geração']);
saveas(f,'MM_R2ajustado.png');

% MM_rmse: MÉDIA dos Melhores / geração
f=figure('name','RMSE: Média dos Melhores / geração','numbertitle','off','visible','off');
plot(mean(Bests.rmse),'r');
ylabel('RMSE: Média dos Melhores');
xlabel('Gerações');
title(['RMSE: Média dos Melhores / geração']);
saveas(f,'MM_rmse.png');

% MM_numreg: MÉDIA dos Melhores / geração
f=figure('name','#Regressores: Média dos Melhores / geração','numbertitle','off','visible','off');
plot(mean(Bests.num_reg),'g');
ylabel('#Regressores: Média dos Melhores');
xlabel('Gerações');
title(['#Regressores: Média dos Melhores / geração']);
saveas(f,'MM_numreg.png');

