% Gráficos

load gp

%% MELHOR Indivíduo

% RMSE: Melhor indivíduo / geração
f=figure('name','RMSE: Melhor indivíduo / geração','numbertitle','off','visible','on');
plot(min(gp.geracoes.rmse),'-ko','LineWidth',2,'MarkerSize',5);
legend('RMSE: Melhor indivíduo');
ylabel('RMSE');
xlabel('Gerações');
title('RMSE: Melhor indivíduo / geração');
saveas(f,'Melhor_RMSE.png');

% Índices dos Melhores indivíduos (em todas as gerações)
menores=min(gp.geracoes.rmse);
r2ajustado_menoresRMSE=[];
reg_menoresRMSE=[];
ind_menoresRMSE=[];
reg_signif_menoresRMSE=[];
propreg_signif_menoresRMSE=[];

for i=1:length(menores)
    ind_menor=find(gp.geracoes.rmse(:,i)==menores(i));
    if length(ind_menor)==1
        r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)];
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
        reg_signif_menoresRMSE=[reg_signif_menoresRMSE gp.geracoes.numreg_signif(ind_menor,i)];
        propreg_signif_menoresRMSE=[propreg_signif_menoresRMSE gp.geracoes.prop_numreg_signif(ind_menor,i)];
    else
        ind_menor=ind_menor(1);
        r2ajustado_menoresRMSE=[r2ajustado_menoresRMSE gp.geracoes.r2ajustado(ind_menor,i)];
        reg_menoresRMSE=[reg_menoresRMSE gp.geracoes.num_reg(ind_menor,i)];
        ind_menoresRMSE=[ind_menoresRMSE ind_menor];
        reg_signif_menoresRMSE=[reg_signif_menoresRMSE gp.geracoes.numreg_signif(ind_menor,i)];
        propreg_signif_menoresRMSE=[propreg_signif_menoresRMSE gp.geracoes.prop_numreg_signif(ind_menor,i)];
    end
end

gp.geracoes.r2ajustado_menoresRMSE=r2ajustado_menoresRMSE;
gp.geracoes.reg_menoresRMSE=reg_menoresRMSE;
gp.geracoes.ind_menoresRMSE=ind_menoresRMSE;
gp.geracoes.reg_signif_menoresRMSE=reg_signif_menoresRMSE;
gp.geracoes.propreg_signif_menoresRMSE=propreg_signif_menoresRMSE;

gp.Bests.R2ajustado=gp.geracoes.r2ajustado_menoresRMSE;
gp.Bests.rmse=menores;
gp.Bests.num_reg=gp.geracoes.reg_menoresRMSE;
gp.Bests.reg_signif=gp.geracoes.reg_signif_menoresRMSE;
gp.Bests.propreg_signif=gp.geracoes.propreg_signif_menoresRMSE;

% R2ajustado: Melhor indivíduo (/geração) e R2ajustado:Comparativo
f=figure('name','R2ajustado: Melhor indivíduo (/geração) e R2ajustado:Comparativo','numbertitle','off','visible','on');
plot(r2ajustado_menoresRMSE,'-go','LineWidth',2,'MarkerSize',5);
legend('R2ajustado: Melhor indivíduo');
hold on;
plot(max(gp.geracoes.r2ajustado),':go','LineWidth',3,'MarkerSize',5);
legend('R2ajustado: Comparativo');
ylabel('R2ajustado');
xlabel('Gerações');
title('R2ajustado: Melhor indivíduo (/geração) e R2ajustado:Comparativo');
hold off;
saveas(f,'MelhorR2ajustado.png');

% #reg E #reg.signif do Melhor indivíduo / geração
f=figure('name','#reg E #reg.signif: Melhor indivíduo / geração','numbertitle','off','visible','on');
plot(reg_menoresRMSE,'-ro','LineWidth',2,'MarkerSize',5);
legend('#reg: Melhor indivíduo');
hold on;
plot(reg_signif_menoresRMSE,':ro','LineWidth',3,'MarkerSize',5);
legend('#reg.signif: Melhor indivíduo');
ylabel('#reg E #reg.signif');
xlabel('Gerações');
title('#reg E #reg.signif: Melhor indivíduo / geração');
hold off;
saveas(f,'Melhor#regE#reg.signif.png');

% Melhor prop.#reg.signif : prop.#reg.signif do Melhor indivíduo / geração
f=figure('name','prop.#reg.signif: Melhor indivíduo / geração','numbertitle','off','visible','on');
plot(propreg_signif_menoresRMSE,'-bo','LineWidth',2,'MarkerSize',5);
legend('prop.#reg.signif: Melhor indivíduo');
ylabel('prop.#reg.signif');
xlabel('Gerações');
title('prop.#reg.signif: Melhor indivíduo / geração');
saveas(f,'Melhorprop.#reg.signif.png');

%% MÉDIA dos Indivíduos

% Média RMSE      : MÉDIA dos indivíduos / geração
f=figure('name','RMSE: Média dos indivíduos / geração','numbertitle','off','visible','on');
plot(mean(gp.geracoes.rmse),'-k','LineWidth',2);
legend('RMSE: Média');
ylabel('RMSE: Média dos indivíduos');
xlabel('Gerações');
title('RMSE: Média dos indivíduos / geração');
saveas(f,'MédiaRMSE.png');

% Média R2ajustado: MÉDIA dos indivíduos / geração
f=figure('name','R2ajustado: Média dos indivíduos / geração','numbertitle','off','visible','on');
plot(mean(gp.geracoes.r2ajustado),'-g','LineWidth',2);
legend('R2ajustado: Média');
ylabel('R2ajustado: Média dos indivíduos');
xlabel('Gerações');
title('R2ajustado: Média dos indivíduos / geração');
saveas(f,'MédiaR2ajustado.png');

% Média #reg e Média #reg.signif: MÉDIA dos indivíduos / geração
f=figure('name','#reg E #reg.signif: Média dos indivíduos / geração','numbertitle','off','visible','on');
plot(mean(gp.geracoes.num_reg),'-r','LineWidth',2);
legend('#reg: Média');
hold on;
plot(mean(gp.geracoes.numreg_signif),':r','LineWidth',3)
legend('#reg.signif: Média');
ylabel('#reg E #reg.signif: Média dos indivíduos');
xlabel('Gerações');
title('#reg E #reg.signif: Média dos indivíduos / geração');
hold off;
saveas(f,'Média#regE#reg.signif.png');

% Média prop.#reg.signif: MÉDIA dos indivíduos / geração
f=figure('name','prop.#reg.signif: Média dos indivíduos / geração','numbertitle','off','visible','on');
plot(mean(gp.geracoes.prop_numreg_signif),'b','LineWidth',2);
legend('prop.#reg.signif: Média');
ylabel('prop.#reg.signif: Média dos indivíduos');
xlabel('Gerações');
title('prop.#reg.signif: Média dos indivíduos / geração');
saveas(f,'Médiaprop.numreg.signif.png');

