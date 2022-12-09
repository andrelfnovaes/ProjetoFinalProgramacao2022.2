% Taxas Adaptativas

gp.runcontrol.num_gen=50;
gp.operators.crossover.p_cross=0.50;
gp.operators.mutation.p_mutate=0.50;

cross=zeros(1,length(gp.runcontrol.num_gen));
cross(1,1)=0.5;

muta=zeros(1,length(gp.runcontrol.num_gen));
muta(1,1)=0.5;

% for count=1:gp.runcontrol.num_gen
% 
%     if count>1 && count<=round(0.40*gp.runcontrol.num_gen)
%         gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.50)/(0.40*gp.runcontrol.num_gen))
%         cross(1,count)=gp.operators.crossover.p_cross;
%         
%         gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.50)/(0.40*gp.runcontrol.num_gen))
%         muta(1,count)=gp.operators.mutation.p_mutate;
% 
%     elseif count>round(0.40*gp.runcontrol.num_gen) && count<=round(0.56*gp.runcontrol.num_gen)
%         gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.60)/((0.56-0.40)*gp.runcontrol.num_gen))
%         cross(1,count)=gp.operators.crossover.p_cross;
%         
%         gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.60)/((0.56-0.40)*gp.runcontrol.num_gen))
%         muta(1,count)=gp.operators.mutation.p_mutate;
% 
%     elseif count>round(0.56*gp.runcontrol.num_gen) && count<=round(0.88*gp.runcontrol.num_gen)
%         gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.60)/((0.88-0.56)*gp.runcontrol.num_gen))
%         cross(1,count)=gp.operators.crossover.p_cross;
%         
%         gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.60)/((0.88-0.56)*gp.runcontrol.num_gen))
%         muta(1,count)=gp.operators.mutation.p_mutate;
% 
%     elseif count>round(0.88*gp.runcontrol.num_gen) && count<=round(1.00*gp.runcontrol.num_gen)
%         gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.60)/((1.00-0.88)*gp.runcontrol.num_gen))
%         cross(1,count)=gp.operators.crossover.p_cross;
%         
%         gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.60)/((1.00-0.88)*gp.runcontrol.num_gen))
%         muta(1,count)=gp.operators.mutation.p_mutate;
% 
%     end
%     
% end

for count=1:gp.runcontrol.num_gen
    
    if count>=1 && count<=round(0.30*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.50)/(0.30*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    elseif count>round(0.30*gp.runcontrol.num_gen) && count<=round(0.40*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.10)/((0.40-0.30)*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    elseif count>round(0.40*gp.runcontrol.num_gen) && count<=round(0.70*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.80-0.10)/((0.70-0.40)*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    elseif count>round(0.70*gp.runcontrol.num_gen) && count<=round(0.80*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.80-0.20)/((0.80-0.70)*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    elseif count>round(0.80*gp.runcontrol.num_gen) && count<=round(0.90*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.50-0.20)/((0.90-0.80)*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    elseif count>round(0.90*gp.runcontrol.num_gen) && count<=round(1.00*gp.runcontrol.num_gen)
        gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
        cross(1,count)=gp.operators.crossover.p_cross;
        
        gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.50-0.40)/((1.00-0.90)*gp.runcontrol.num_gen));
        muta(1,count)=gp.operators.mutation.p_mutate;

    end

    
end

f=figure('name','Crossover e Mutação','numbertitle','off','visible','on');
plot(cross,'k');
hold on;
plot(muta,'b');
ylabel('Crossover e Mutação');
xlabel('Gerações');
legend('Crossover','Mutação');
title(['Crossover e Mutação / geração']);
hold off;

%% Crossover e Mutação Adaptativos (ANTIGO)

% if count>=1 && count<=round(0.40*gp.runcontrol.num_gen)
%     count
%     gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.50)/(0.40*gp.runcontrol.num_gen));
%     gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.50)/(0.40*gp.runcontrol.num_gen));
% 
% elseif count>round(0.40*gp.runcontrol.num_gen) && count<=round(0.56*gp.runcontrol.num_gen)
%     count
%     gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.60)/((0.56-0.40)*gp.runcontrol.num_gen));
%     gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.60)/((0.56-0.40)*gp.runcontrol.num_gen));
% 
% elseif count>round(0.56*gp.runcontrol.num_gen) && count<=round(0.88*gp.runcontrol.num_gen)
%     count
%     gp.operators.crossover.p_cross=gp.operators.crossover.p_cross+((0.90-0.60)/((0.88-0.56)*gp.runcontrol.num_gen));
%     gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate-((0.90-0.60)/((0.88-0.56)*gp.runcontrol.num_gen));
% 
% elseif count>round(0.88*gp.runcontrol.num_gen) && count<=round(1.00*gp.runcontrol.num_gen)
%     count
%     gp.operators.crossover.p_cross=gp.operators.crossover.p_cross-((0.90-0.60)/((1.00-0.88)*gp.runcontrol.num_gen));
%     gp.operators.mutation.p_mutate=gp.operators.mutation.p_mutate+((0.90-0.60)/((1.00-0.88)*gp.runcontrol.num_gen));
% 
% end
