function [ updatedPopulation ] =  InsertBestIndividual (population, ...
    bestIndividualIndex, nCopies)

updatedPopulation = population;

for j = 1:nCopies
    
    updatedPopulation(j,:) = updatedPopulation(bestIndividualIndex,:);

end
end

