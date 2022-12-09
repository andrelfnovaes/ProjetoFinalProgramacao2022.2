% Gr�ficos

load gp

%% MELHOR Indiv�duo

% RMSE: Melhor indiv�duo / gera��o
f=figure('name','RMSE: Melhor indiv�duo / gera��o','numbertitle','off','visible','on');
plot(min(gp.geracoes.rmse),'-ko','LineWidth',2,'MarkerSize',5);
legend('RMSE: Melhor indiv�duo');
ylabel('RMSE');
xlabel('Gera��es');
title('RMSE: Melhor indiv�duo / gera��o');
saveas(f,'Melhor_RMSE.png');

% �ndices dos Melhores indiv�duos (em todas as gera��es)
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

% R2ajustado: Melhor indiv�duo (/gera��o) e R2ajustado:Comparativo
f=figure('name','R2ajustado: Melhor indiv�duo (/gera��o) e R2ajustado:Comparativo','numbertitle','off','visible','on');
plot(r2ajustado_menoresRMSE,'-go','LineWidth',2,'MarkerSize',5);
legend('R2ajustado: Melhor indiv�duo');
hold on;
plot(max(gp.geracoes.r2ajustado),':go','LineWidth',3,'MarkerSize',5);
legend('R2ajustado: Comparativo');
ylabel('R2ajustado');
xlabel('Gera��es');
title('R2ajustado: Melhor indiv�duo (/gera��o) e R2ajustado:Comparativo');
hold off;
saveas(f,'MelhorR2ajustado.png');

% #reg E #reg.signif do Melhor indiv�duo / gera��o
f=figure('name','#reg E #reg.signif: Melhor indiv�duo / gera��o','numbertitle','off','visible','on');
plot(reg_menoresRMSE,'-ro','LineWidth',2,'MarkerSize',5);
legend('#reg: Melhor indiv�duo');
hold on;
plot(reg_signif_menoresRMSE,':ro','LineWidth',3,'MarkerSize',5);
legend('#reg.signif: Melhor indiv�duo');
ylabel('#reg E #reg.signif');
xlabel('Gera��es');
title('#reg E #reg.signif: Melhor indiv�duo / gera��o');
hold off;
saveas(f,'Melhor#regE#reg.signif.png');

% Melhor prop.#reg.signif : prop.#reg.signif do Melhor indiv�duo / gera��o
f=figure('name','prop.#reg.signif: Melhor indiv�duo / gera��o','numbertitle','off','visible','on');
plot(propreg_signif_menoresRMSE,'-bo','LineWidth',2,'MarkerSize',5);
legend('prop.#reg.signif: Melhor indiv�duo');
ylabel('prop.#reg.signif');
xlabel('Gera��es');
title('prop.#reg.signif: Melhor indiv�duo / gera��o');
saveas(f,'Melhorprop.#reg.signif.png');

%% M�DIA dos Indiv�duos

% M�dia RMSE      : M�DIA dos indiv�duos / gera��o
f=figure('name','RMSE: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','on');
plot(mean(gp.geracoes.rmse),'-k','LineWidth',2);
legend('RMSE: M�dia');
ylabel('RMSE: M�dia dos indiv�duos');
xlabel('Gera��es');
title('RMSE: M�dia dos indiv�duos / gera��o');
saveas(f,'M�diaRMSE.png');

% M�dia R2ajustado: M�DIA dos indiv�duos / gera��o
f=figure('name','R2ajustado: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','on');
plot(mean(gp.geracoes.r2ajustado),'-g','LineWidth',2);
legend('R2ajustado: M�dia');
ylabel('R2ajustado: M�dia dos indiv�duos');
xlabel('Gera��es');
title('R2ajustado: M�dia dos indiv�duos / gera��o');
saveas(f,'M�diaR2ajustado.png');

% M�dia #reg e M�dia #reg.signif: M�DIA dos indiv�duos / gera��o
f=figure('name','#reg E #reg.signif: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','on');
plot(mean(gp.geracoes.num_reg),'-r','LineWidth',2);
legend('#reg: M�dia');
hold on;
plot(mean(gp.geracoes.numreg_signif),':r','LineWidth',3)
legend('#reg.signif: M�dia');
ylabel('#reg E #reg.signif: M�dia dos indiv�duos');
xlabel('Gera��es');
title('#reg E #reg.signif: M�dia dos indiv�duos / gera��o');
hold off;
saveas(f,'M�dia#regE#reg.signif.png');

% M�dia prop.#reg.signif: M�DIA dos indiv�duos / gera��o
f=figure('name','prop.#reg.signif: M�dia dos indiv�duos / gera��o','numbertitle','off','visible','on');
plot(mean(gp.geracoes.prop_numreg_signif),'b','LineWidth',2);
legend('prop.#reg.signif: M�dia');
ylabel('prop.#reg.signif: M�dia dos indiv�duos');
xlabel('Gera��es');
title('prop.#reg.signif: M�dia dos indiv�duos / gera��o');
saveas(f,'M�diaprop.numreg.signif.png');

