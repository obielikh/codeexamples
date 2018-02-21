function [ updatedPopulation ] =  InsertBestIndividual (population, ...
    bestIndividualIndex, nCopies)

updatedPopulation = population;

for i = 1:nCopies
    
    updatedPopulation(i,:) = population(bestIndividualIndex,:);

end
end

