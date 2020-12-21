clear
clc

load DataName.mat
load DataNum.mat
%% Initialize Constants
Num = Data2;
Name = Data1;

L = Num(:,15);            %[ft] longest flow path
Sl = Num(:,12);           %[unitless] average slope along the longest flow path
P = Num(:,17);            %[in] rainfall depth (tc value in hrs at 50 yr event form NRCC)
A = Num(:,10)*1000^2;     %[m^2]                   
CN = Num(:,7);            % Curve Number
N = length(L);


%% CN Synthetic Triangular Hydrograph Method
tc = zeros(N,1);     %[hrs] Kirpich (1940)
S = zeros(N,1);      %[in] watershed storage
Ia = zeros(N,1);     %[in] initial abstraction
Q_tri = zeros(N,1);  %[m] runoff depth
qp_tri = zeros(N,1); %[m^3/s] 
x_tri = zeros(N,3);
y_tri = zeros(N,3);


for i = 1:N
    tc(i) = (0.0078 * L(i)^(0.77) * Sl(i)^(-0.385)) / 60; 
    S(i) = (1000/CN(i)-10)*1.42;
    Ia(i) = 0.05*S(i);
    Q_tri(i) = (P(i)-Ia(i))^2 / (P(i)-Ia(i)+S(i))*0.0254; 
    qp_tri(i) = 2*Q_tri(i)*A(i) / (2.937*tc(i)*3600);
    
end

tp_tri = 1.1*tc;                 %[hrs] peak time, Kirpich
tr_tri = 1.67*tp_tri;            %[hrs] recession time, Kirpich

for j = 1:N
    x_tri(j,2) = tp_tri(j);
    x_tri(j,3) = tp_tri(j)+tr_tri(j);
    y_tri(j,2) = qp_tri(j);
    plot(x_tri(j,:), y_tri(j,:))
    hold on
end

hold off

legend(Name(1,1), Name(2), Name(3), Name(4), Name(5), Name(6), Name(7),...
    Name(8), Name(9), Name(10), Name(11), Name(12), Name(13),Name(14),...
    Name(15), Name(16), Name(17), Name(18), Name(19), Name(20), Name(21))
xlabel('Time [hrs]')
ylabel('Storm Runoff [m^3/s]')
title('Synthetic Triangular Hydrograph of 21 Different Stream Crossings for a 50 year Return Period')

