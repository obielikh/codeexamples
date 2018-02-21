function iSelected = TournamentSelect(fitness,tournamentSelectionParameter,...
    tournamentSize)

%tournamentSelectionParameter=0.75;
%tournamentSize=3;
populationSize = size(fitness,1);
    
    for i = 1:tournamentSize
        iTmp(i) = 1 + fix(rand*populationSize);
    end
    
   
    iSelected=0;
    
    while length(iTmp)>1
        r = rand;
        if (r < tournamentSelectionParameter)
            iSelected = max(iTmp);
            break;
        else
            [m,maxPos] = max (iTmp);
            iTmp(maxPos) = [];
            
        end
    end
    
    if length (iTmp)==1
        iSelected = iTmp(1);
    end
 
end

