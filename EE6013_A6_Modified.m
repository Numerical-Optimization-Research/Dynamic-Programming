%% EE6013_A6.m


%% Clean Up

clear all;
close all;
clc


%% Variables
global a  b c d e kT Tref kSoC socRef kt socMax socMin
a = 5.7905;
b = -6.8292;
c = 3.3209;
d = 5.3696*10^-1;
e = 6.1638e-2;
kt = 1.49e-6;
kT = 0.0693;
Tref = 0.25;
kSoC = 1.04;
socRef = 0.5;
socMax = 0.98;
socMin = 0.2;

initial = 0.5;
N = 1000;
delta = 0.1;
Ceol = 800000;

%% Initialize

% x = zeros(1,N);
x(1) = socMax - socMin;     % DoD
x(2) = 1;     % Temp
x(3) = 0.85;    % SoC
x(4) = 1;    % time
u = -0.01:0.001:0.01;


%% Main

% First run works out to figure out the costs at each stage of the equation


u = nextU(x);
outcost = recursive(x,N);


% It doesn't make sense to have a zero state as the answer
% although it would be the minimum it is not practical.
index(1) = minimumX_nonzero(xMin(1,:));
index(2) = minimumX_nonzero(xMin(2,:));
index(3) = minimumX_nonzero(xMin(3,:));
index(4) = minimumX_nonzero(xMin(4,:));
indexFinal = min(index(:));

xfinal(1) = xMin(1,indexFinal);
xfinal(2) = xMin(2,indexFinal);
xfinal(3) = xMin(3,indexFinal);
xfinal(4) = xMin(4,indexFinal);

plotDOD(xfinal);
plotCost(outcost,u,N);

fd_final = evalFd(xfinal(:));
deg_dot = -kt*exp(kSoC*(xfinal(3) - socRef))*exp(kT*(x(2) - Tref)*(Tref/x(2)))*exp(-fd_final);
totalCost = -Ceol*deg_dot;


%% Functions


function plotCost(outcost,u,N)

    figure('name','Plot of ');
    plot(u,outcost(1,:),'LineWidth',1.5);
    hold on
    for i = 2:N
        plot(u,outcost(N,:),'LineWidth',1.5);
    end
    hold off
    grid on;
%     xlim([0 1])
%     ylim([0 1])
%     title('Discretized Battery Power');
    xlabel('Control Action u');
    ylabel('Recursive Function Value');
end


function plotDOD(xfinal)

    x = -1:0.001:1;
    func = zeros(4,length(x));
    choicefunc = zeros(1,4);
    for i = 1:4
        for j = 1:length(x)
            func(i,j) = eval(x(j),i);
        end
        choicefunc(i) = eval(xfinal(i),i);
        figure('name','Plot of Functions');
        
        plot(x,func(i,:),xfinal(i),choicefunc(i),'b*','LineWidth',1.5);
        switch i
            case 1
               title('Depth of Discharge Stress Factor'); 
            case 2
               title('Temperature Stress Factor');
            case 3
               title('State of Charge Stress Factor');
            case 4
               title('Time Stress Factor');
        end
        grid on;
    end
    
end


function x = minimumX_nonzero(xMin)
    x = 10000;
    for i = 1:length(xMin)
        if (xMin(i) < x)&&(xMin(i) >= 0)
            x = i;
        end
    end
    
end


function xmin = enforceBC(xmin)

    for i = 1:length(xmin)
        if xmin(i) < 0 
            xmin(i) = 0;
        elseif xmin(i) > 1
            xmin(i) = 1;
        end
    end

end


function xout = updateX(umin,x)

    xout = zeros(1,length(umin)+1);
    for i = 1:length(umin)+1
        if i == 1
            xout(i) = x;
        else
            xout(i) = xout(i-1) + umin(i-1);
        end
    end

end


function minimum = findMin(outcost,N,u)

    u_min = zeros(1,N-1);
    for i = 1:N-1

%         minU = min(outcost(i,:));
        u_min(i) = u(find(outcost(i,:) > 0,1,'first'));

    end
    minimum = u_min;
    
end


function cost = recursive(x,N,u)

cost = zeros(1,N);
    for i = N:-1:1
        
        if i == N
            cost(i) = eval(x,3)*eval(x,2)*eval(x,4);
        else
            cost(i) = cost(i+1) + u(i);
        end
        
    end
    
end


function u = nextU(x)
    
    u = eval(x,3)*eval(x,2)*eval(x,4);

end


function funcVal = eval(x,i)
global a b c d e kt Tref kSoC socRef kT 

    switch i
        case 1
            funcVal = a*x(1)^4 + b*x(1)^3 + c*x(1)^2 + d*x(1) + e;
        case 2
            funcVal = exp(kT*(x(2) - Tref)*(Tref/x(2)));
        case 3
            funcVal = exp(kSoC*(x(3) - socRef));
        case 4
            funcVal = kt*x(4);
    end
end


function fd = evalFd(x)
global a b c d e kt Tref kSoC socRef kT

    funcVal(1) = a*x(1)^4 + b*x(1)^3 + c*x(1)^2 + d*x(1) + e;
    funcVal(2) = exp(kT*(x(2) - Tref)*(Tref/x(2)));
    funcVal(3) = exp(kSoC*(x(3) - socRef));
    funcVal(4) = kt*x(4);
    fd = funcVal(4)*funcVal(3)*funcVal(2) + funcVal(1)*funcVal(3)*funcVal(2);

end


%% Edits
%
%
%
%
%