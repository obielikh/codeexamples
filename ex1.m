%One-dimensional Kohonen network

clear all
close all
%generating patterns in triangle
N = 1000;                    % # Number of points
V = [-0.5, 0; 0, 0.866; 0.5, 0];  % # Triangle vertices, pairs of (x, y)
t = sqrt(rand(N, 1));
s = rand(N, 1);
P = (1 - t) * V(1, :) + bsxfun(@times, ((1 - s) * V(2, :) + s * V(3, :)), t);
for i = 1:size(P)
    P(i,1) = P(i,1) + 0.5;
end
%scatter(P(:, 1), P(:, 2), '.')
hold on

%ordering phase variables
torder=1000;
t=[1:1:torder];
tau=200;
sigma0=100;
sigmao=sigma0*exp(-t/tau);
eta0=0.1;
etao=eta0*exp(-t/tau);

%convergence phase variables
sigmac=0.9;
etac=0.01;
tconv=40000;

%random initialization of weights
%w=randi(100,100,2)./100;
w=randi(20,100,2)./100+0.4;
%scatter(w(:, 1), w(:, 2), '*')

%ordering phase
for time=t
    
    
    
    %picking the pattern
    pat=randi(N,1);
    pickedPattern=(P(pat,:));
    
    %finding the winner
    for i=[1:100]
        b(i)=sqrt((pickedPattern(1)-w(i,1))^2+(pickedPattern(2)-w(i,2))^2);
        [bmin,i0]=min(b);
    end
    %Kohonen step
    for i=[1:100]
        L(i)=exp(-(i-i0)^2/(2*(sigmao(time))^2));
        deltaw(i,1)=etao(time)*L(i)*(pickedPattern(1)-w(i,1));
        deltaw(i,2)=etao(time)*L(i)*(pickedPattern(2)-w(i,2));
    end
    w=w+deltaw;
    
end

%plotting
subplot(1,2,1)
scatter(w(:, 1), w(:, 2), '.')
hold on
line(w(:, 1), w(:, 2))
hold on
line([0,1], [0, 0], 'Color', 'r');
hold on
line([1,0.5], [0, 0.866], 'Color', 'r');
hold on
line([0.5,0], [0.866, 0], 'Color', 'r');
t=[1:1:tconv];
hold off

%convergence phase
for time=t
    
    %picking the pattern
    pat=randi(N,1);
    pickedPattern=(P(pat,:));
    
    %finding the winner
    for i=[1:100]
        b(i)=sqrt((pickedPattern(1)-w(i,1))^2+(pickedPattern(2)-w(i,2))^2);
        [bmin,i0]=min(b);
    end
    
    %Kohonen step
    for i=[1:100]
        L(i)=exp(-(i-i0)^2/(2*sigmac^2));
        deltaw(i,1)=etac*L(i)*(pickedPattern(1)-w(i,1));
        deltaw(i,2)=etac*L(i)*(pickedPattern(2)-w(i,2));
    end
    w=w+deltaw;
end

%plotting
subplot(1,2,2)
scatter(w(:, 1), w(:, 2), '.')
hold on
line(w(:, 1), w(:, 2))
hold on
line([0,1], [0, 0], 'Color', 'r');
hold on
line([1,0.5], [0, 0.866], 'Color', 'r');
hold on
line([0.5,0], [0.866, 0], 'Color', 'r');
