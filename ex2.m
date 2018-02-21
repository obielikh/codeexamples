%Classification using Kohonen’s algorithm

clear all
close all

%reading data
a=load ('wine.data.txt');
inputPattern = a(:,2:14);
N=length(inputPattern);
inputPattern=zscore(inputPattern);

%ordering phase variables
torder=1000;
t=[1:1:torder];
tau=300;
sigma0=30;
sigmao=sigma0*exp(-t/tau);
eta0=0.1;
etao=eta0*exp(-t/tau);

%convergence phase variables
sigmac=0.9;
etac=0.01;
tconv=20000;

%random initialization of weights
for i=1:20
    for j=1:20
        for k=1:13
            w(i,j,k)=2*rand-1;
        end
    end
end

%ordering phase
for time=t
    
    
    
    %picking the pattern
    pat=randi(N,1);
    pickedPattern=(inputPattern(pat,:));
    b=zeros(20,20);
    %finding the winner
    for i=1:20
        for j=1:20
            for k=1:13
                b(i,j)=b(i,j)+(pickedPattern(k)-w(i,j,k))^2;
            end
            b(i,j)=sqrt(b(i,j));
        end
    end
    bmin=min(b(:));
    [i0,j0]=find (b==bmin);
    %Kohonen step
    for i=1:20
        for j=1:20
            L(i,j)=exp((-(i-i0)^2-(j-j0)^2)/(2*(sigmao(time))^2));
            for k=1:13
                deltaw(i,j,k)=etao(time)*L(i,j)*(pickedPattern(k)-w(i,j,k));
            end
        end
    end
    w=w+deltaw;
    
end

%convergence phase
t=[1:1:tconv];
for time=t
    
    
    
    %picking the pattern
    pat=randi(N,1);
    pickedPattern=(inputPattern(pat,:));
    b=zeros(20,20);
    %finding the winner
    for i=1:20
        for j=1:20
            for k=1:13
                b(i,j)=b(i,j)+(pickedPattern(k)-w(i,j,k))^2;
            end
            b(i,j)=sqrt(b(i,j));
        end
    end
    bmin=min(b(:));
    [i0,j0]=find (b==bmin);
    %Kohonen step
    for i=1:20
        for j=1:20
            L(i,j)=exp((-(i-i0)^2-(j-j0)^2)/(2*(sigmac)^2));
            for k=1:13
                deltaw(i,j,k)=etac*L(i,j)*(pickedPattern(k)-w(i,j,k));
            end
        end
    end
    w=w+deltaw;
    
end
%colouring
colours=a(:,1);
colouredWeights=zeros(20,20);
for z=1:N
    pickedPattern=(inputPattern(z,:));
    for i=1:20
        for j=1:20
            for k=1:13
                b(i,j)=b(i,j)+(pickedPattern(k)-w(i,j,k))^2;
            end
            b(i,j)=sqrt(b(i,j));
        end
    end
    bmin=min(b(:));
    [i0,j0]=find (b==bmin);
    colouredWeights(i0,j0)=colours(z);
end
%plotting
for i=1:20
    for j=1:20
        switch colouredWeights(i,j)
            case 0
                plot(i,j,'k.','MarkerSize',20)
                hold on
            case 1
                plot(i,j,'r.','MarkerSize',20)
                hold on
            case 2
                plot(i,j,'g.','MarkerSize',20)
                hold on
            case 3
                plot(i,j,'b.','MarkerSize',20)
                hold on
        end
    end
end

