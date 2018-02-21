%Selection  drift in an isolated population

clear all
close all

%N = 10^3;
s = 10^(-3);
ring = zeros(6,1);
iters=1;
%
for N=[2*10^2 5*10^2 10^3 5*10^3 10^4 10^5]
    for numberOfExperiments=1:10^3
        %         if rem(numberOfExperiments,100)==0
        %             numberOfExperiments
        %         end
        individuals=zeros (N,2);
        individuals(randi(N),1)=1;
        fAA=0;
        fAa=1;
        faa=N-1;
        time=0;
        while (fAA~=N) && (faa~=N)
            time=time+1;
            %average fitness
            favg=(1-s)*faa+(1-(s/2))*fAa+fAA;
            %proportion of gametes in gamete pool
            pa=((1-(s/2))*fAa+2*(1-s)*faa)/favg;
            pA=((1-(s/2))*fAa+2*fAA)/favg;
            %normalization
            p=pA+pa;
            pa=pa/p;
            pA=pA/p;
            %nullify frequencies
            fAA=0;
            fAa=0;
            faa=0;
            %creating new generation, calculating new frequencies
            for i=1:N
                for j=1:2
                    if rand>pa
                        individuals(i,j)=1;
                    else
                        individuals(i,j)=0;
                    end
                end
                if (individuals(i,1)==1) && (individuals(i,2)==1)
                    fAA=fAA+1;
                elseif (individuals(i,1)==1) && (individuals(i,2)==0)
                    fAa=fAa+1;
                elseif (individuals(i,1)==0) && (individuals(i,2)==1)
                    fAa=fAa+1;
                else
                    faa=faa+1;
                end
            end
        end
        if fAA==N
            ring(iters)=ring(iters)+1;
        end
    end
    psim(iters)=ring(iters)/numberOfExperiments;
    iters=iters+1;
end

N=[2*10^2 5*10^2 10^3 5*10^3 10^4 10^5];

%my simulation
loglog (N,psim,'*g')
hold on

%Kimura
%p0=1/(2*N);
for i=1:length(N)
    p0=1/(2*N(i));
    pfix(i) = (1-exp(-2*s*N(i)*p0))/(1-exp(-2*s*N(i)));
end
loglog (N,pfix,'*b')
hold on

%Haldane
loglog (N,s,'*r')

title('Fixation probability of allele A for different values of N');
h = legend('Our simulation','$Kimura','$Haldane');
ylabel('Probability')
xlabel('Value of N')
ylim([0.5*10^(-3),0.5*10^(-1)])