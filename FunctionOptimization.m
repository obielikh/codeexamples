%Genetic algorithm with tournament selection, crossover, mutation, elitism

clear all

%initialization
populationSize = 100;
numberOfGenes = 50;
crossoverProbability = 0.9;
mutationProbability = 0.035;
tournamentSelectionParameter = 0.75;
tournamentSize=2;
variableRange=5.0;
numberOfGenerations = 200;
fitness = zeros(populationSize,1);
numberOfVariables = 2;
nCopies = 1;
numberOfRuns = 20;


for iRuns = 1:numberOfRuns
    
    population = InitializePopulation(populationSize, numberOfGenes);
    
    for iGeneration = 1:numberOfGenerations
        
        maximumFitness = 0.0;
        xBest = zeros (1,2);
        bestIndividualIndex = 0;
        
        %finding the individual with the highest fitness
        for i = 1:populationSize
            chromosome = population(i,:);
            x = DecodeChromosome(chromosome, numberOfVariables, variableRange);
            decodedPopulation(i,:)=x;
            fitness(i) = EvaluateIndividual(x);
            if (fitness(i) > maximumFitness)
                maximumFitness =fitness(i);
                bestIndividualIndex=i;
                xBest=x;
            end
        end
             
        tempPopulation = population;
        
        %tournament selection
        for i=1:2:populationSize
            i1 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
            i2 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
            chromosome1 = population(i1,:);
            chromosome2 = population(i2,:);
       
            %crossover
            r=rand;
            if (r<crossoverProbability)
                newChromosomePair = Cross (chromosome1,chromosome2);
                tempPopulation(i,:) = newChromosomePair(1,:);
                tempPopulation(i+1,:) = newChromosomePair(2,:);
            else
                tempPopulation(i,:) = chromosome1;
                tempPopulation(i+1,:) = chromosome2;
            end
        end % Loop over population
        
        %mutation
        for i = 1:populationSize
            originalChromosome = tempPopulation(i,:);
            mutatedChromosome = Mutate(originalChromosome,mutationProbability);
            tempPopulation(i,:)=mutatedChromosome;
        end

        %elitism
        population = InsertBestIndividual (population, bestIndividualIndex, nCopies);
        
        for j=1:nCopies
            tempPopulation(j,:) = population(j,:);
        end
        population = tempPopulation;
        
        
        avgxBest(iGeneration,:) = xBest;
        avgminVal(iGeneration) = 1/maximumFitness;  
    end %Loop over generations 
    
end

averagexBest = mean (avgxBest);
averageMinimumValue = mean (avgminVal);

