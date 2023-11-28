%% EE6013_A6.m


%% Clean Up

clear all;
close all;
clc


%% Variables
global a  b c d e kT Tref kSoC socRef kt socMax socMin dodMin dodMax ...
    minVal
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
dodMax = 0.85;

dodMin = 0.65;
minVal = 0.25;

initial = 0.5;
N = 1000;
delta = 0.1;
Ceol = 800000;

%% Initialize

% x = zeros(1,N);
x(1) = socMax - socMin;     % DoD
x(2) = 1;     % Temp
x(3) = socMin;    % SoC
x(4) = 1;    % time
u = -0.1:0.001:0.1;


%% Main

% First run works out to figure out the costs at each stage of the equation

outcost = recursive(x(2),u,N,2);
uMin(2,:) = findMin(outcost,N,u);
xMin(2,:) = updateX(uMin(2,:),x(2),0);
% xMin(2,:) = enforceBC(xMin(2,:));


% SOC 
outcost = recursive(x(3),u,N,3);
uMin(3,:) = findMin(outcost,N,u);
xMin(3,:) = updateX(uMin(3,:),x(3),1);
index(3) = minimumX_nonzero(xMin(3,:));
socMinNew = xMin(3,index(3));
% xMin(3,:) = enforceBC(xMin(3,:));


% DOD depends on SOC so calculate SOC first?
x(1) = socMax - socMinNew;
outcost = recursive(x(1),u,N,1);
uMin(1,:) = findMin(outcost,N,u);
xMin(1,:) = updateX(uMin(1,:),x(1),2);
% xMin(1,:) = enforceBC(xMin(1,:));


outcost = recursive(x(4),u,N,4);
uMin(4,:) = findMin(outcost,N,u);
xMin(4,:) = updateX(uMin(4,:),x(4),0);
% xMin(4,:) = enforceBC(xMin(4,:));


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
totalCost = -Ceol*deg_dot*60*24;
fprintf('Final Cost at the optimal values for 1 day is:\n');
fprintf("$ %6.2f\n\n",totalCost);
% Since the equation already includes seconds I think, it makes sense to
% include a conversion to days, so 60min/hour 24hour/day 1day
% Otherwise the cost will be a decimal.

fd  = -2:0.001:4;
e_fd = exp(-fd);
figure('name','fd versus exp(fd)');
plot(fd,e_fd,fd,fd,'LineWidth',1.5);
xlabel('Linearly Increasing fd');
ylabel('exp(-fd)');
grid on
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
global dodMin dodMax socMin socMax minVal
    x = -1:0.001:1;
    func = zeros(4,length(x));
    choicefunc = zeros(1,4);
    
    figure('name','Plot of Functions and their Optimum');
    hold on
    for i = 1:4
        for j = 1:length(x)
            func(i,j) = eval(x(j),i);
        end
        choicefunc(i) = eval(xfinal(i),i);
        
        subplot(2,2,i)
        switch i
            case 1
               plot(x,func(i,:),xfinal(i),choicefunc(i),'ks',dodMin,eval(dodMin,1),'r+',dodMax,eval(dodMax,1),'r+','LineWidth',1.5);
               
               title('Depth of Discharge Stress Factor'); 
               ylim([0 2])
            case 2
               plot(x,func(i,:),xfinal(i),choicefunc(i),'ks',minVal,eval(minVal,2),'r+',1,eval(1,2),'r+','LineWidth',1.5);
               title('Temperature Stress Factor');
               ylim([0 5])
            case 3
               plot(x,func(i,:),xfinal(i),choicefunc(i),'ks',socMin,eval(socMin,3),'r+',socMax,eval(socMax,3),'r+','LineWidth',1.5);
               title('State of Charge Stress Factor');
               ylim([0 2])
               xlim([-1 1])
            case 4
               plot(x,func(i,:),xfinal(i),choicefunc(i),'ks',minVal,eval(minVal,4),'r+',1,eval(1,4),'r+','LineWidth',1.5);
               title('Time Stress Factor');
               ylim([-0.5e-6 2e-6])
        end
        grid on;
    end
    hold off
    
end


function x = minimumX_nonzero(xMin)
    x = 10000;
    for i = 1:length(xMin)
        if (xMin(i) < x)&&(xMin(i) >= 0)
            x = i;
        end
    end
    
end


function xmin = enforceBCSOC(xmin)
global socMin socMax
    for i = 1:length(xmin)
        if xmin(i) < socMin 
            xmin(i) = socMin;
        elseif xmin(i) > socMax
            xmin(i) = socMax;
        end
    end

end

function xmin = enforceBCMin(xmin)
global minVal
    for i = 1:length(xmin)
        if xmin(i) < minVal 
            xmin(i) = minVal;
        elseif xmin(i) > 1
            xmin(i) = 1;
        end
    end

end


function xmin = enforceBCDOD(xmin)
global dodMin dodMax
    for i = 1:length(xmin)
        if xmin(i) < dodMin 
            xmin(i) = dodMin;
        elseif xmin(i) > dodMax
            xmin(i) = dodMax;
        end
    end

end


function xout = updateX(umin,x,soc)

    xout = zeros(1,length(umin));
    for i = 1:length(umin)
        if i == 1
            xout(i) = x;
        else
            xout(i) = xout(i-1) + umin(i-1);
        end
        
        if soc == 1
            xout(i) = enforceBCSOC(xout(i));
        elseif soc == 2
            xout(i) = enforceBCDOD(xout(i));
        else
            xout(i) = enforceBCMin(xout(i));
        end
    end
    
    
end


function minimum = findMin(outcost,N,u)

    u_min = zeros(1,N);
    for i = 1:N

%         minU = min(outcost(i,:));
        u_min(i) = u(find(outcost(i,:) > 0,1,'first'));

    end
    minimum = u_min;
    
end


function cost = recursive(x,u,N,j)

cost = zeros(N,length(u));
    for i = N:-1:1
        
        if i == N
            cost(i,:) = eval(x,j)*ones(1,length(u));
        else
            cost(i,:) = cost(i+1,:) + u.^3;
        end
        
    end
    
end


function funcVal = eval(x,i)
global a b c d e kt Tref kSoC socRef kT 

    switch i
        case 1
            funcVal = a*x^4 + b*x^3 + c*x^2 + d*x + e;
        case 2
            funcVal = exp(kT*(x - Tref)*(Tref/x));
        case 3
            funcVal = exp(kSoC*(x - socRef));
        case 4
            funcVal = kt*x;
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