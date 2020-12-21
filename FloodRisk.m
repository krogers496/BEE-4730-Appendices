clear
clc

load DataName.mat
load DataNum.mat
%% Initialize Constants
Num = Data2;
Name = Data1;

L = Num(:,15);            %[ft] longest flow path
Sl = Num(:,12);           %[unitless] average slope along the longest flow path
A = Num(:,10)*1000^2;     %[m^2]                   
Q_III = Num(:,18);
CN = Num(:,7);            % Curve Number
N = length(L);


%% CN Synthetic Triangular Hydrograph Method
tc = zeros(N,1);     %[hrs] Kirpich (1940)
S = zeros(N,1);      %[in] watershed storage
Ia = zeros(N,1);     %[in] initial abstraction
Q_tri = zeros(N,1);  %[m] runoff depth


for i = 1:N
    S(i) = (1000/CN(i)-10)*1.42;
    Ia(i) = 0.05*S(i);
    tc(i) = (0.0078 * L(i)^(0.77) * Sl(i)^(-0.385)) / 60; 
    Q_tri(i) = (Q_III(i)*2.937*tc(i)*3600) / (2*A(i));     %(P(i)-Ia)^2 / (P(i)-Ia+S)*0.0254; 
    syms x
    eqn = Q_tri(i)-(x-Ia(i))^2 / (x-Ia(i)+S(i)) *0.0254 == 0;
    y = solve(eqn,x,'Real',true);
    if y(1)>0
        P(i)=y(1);
    else
        P(i)=y(2);
    end
end

