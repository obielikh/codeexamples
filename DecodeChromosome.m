function x = DecodeChromosome(chromosome,numberOfVariables,variableRange)

    nGenes = size(chromosome,2);
    nPart = fix(nGenes/numberOfVariables);

    for i = 1:numberOfVariables
        x(i) = 0.0;
        for j = 1:nPart
            x(i) = x(i) + chromosome(j+(i-1)*nPart)*2^(-j);
        end
        x(i) = -variableRange + 2*variableRange*x(i)/(1 - 2^(-nPart));
    
    end
    
end