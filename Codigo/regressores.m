function [genes,num_reg]=regressores(individuo)

% load ind
% individuo=ind1
% individuo=real2

genes=cell(1);
j=1;
length(individuo)

%% FUNÇÃO 1: Criando vetor de regressores (=genes)

for i=1:length(individuo)
    if individuo(i)=='x'
        temp=[];
        i;
        k=i;
        while individuo(k)~=' '
            temp=[temp individuo(k)];
            k=k+1;
            if k > length(individuo)
                break
            end
        end
        genes{j}=temp;
        j=j+1;
    end
end

%% FUNÇÃO 2: Excluindo regressores ilegais

tamanho_genes=0;
for j=1:length(genes)
    
    temp=genes{j};
    
    if ~isempty(temp)
        tamanho_genes=tamanho_genes+1;
    end
    
    numero_x=length(find(temp=='x'));
    if numero_x>1
        for k=(j+1):(j+numero_x-1)
            genes{:,k}=[];
        end
    end
    
end

%% FUNÇÃO 3: ReSize de genes, tirando matrizes vazias

genes_temp=cell(1,tamanho_genes);

k=1;
for j=1:length(genes)
    temp=genes{j};
    if ~isempty(temp)
       genes_temp{k}=temp;
       k=k+1;
    end
end
genes=genes_temp;
num_reg=length(genes);

%% FUNÇÃO 5: Preparando genes SÓ para eval
% (Substitui a FUNÇÃO 4: Transformando notação In-fixada em Pré-fixada)

% % Troca o x
% pat='x(\d+)'
% genes=regexprep(genes,pat,'gp.userdata.xtrain(:,$1)')
% 
% % Troca o '^'
% pat='^'
% genes=strrep(genes,pat,'.^')
% 
% % Troca o '*'
% pat='*'
% genes=strrep(genes,pat,'.*')


end