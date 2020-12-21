%%
clear
clc

load culvertdata.mat
load culvertshape.mat
%%Initialize Constants
Num = culvertdata1;
shape = CapstoneStreamCrossingsS1;

Width = Num(:,1)*0.3048; %[m] culvert width 
Height = Num(:,2)*0.3048; %[m] culvert height
HW = Num(:,3)*0.3048; %[m] above top of culvert to road surface
S = Num(:,5);   %slope of culvert
L = Num(:,6)*0.3048;   %[m] length of culvert
Ke = Num(:,7);  % Ke value for each culvert
n = Num(:,8);   % Manning's n for each culvert 
c = Num(:,9);   % Tabulated constant for FHA equations
Y = Num(:,10);  % Tabulated constant for FHA equations
g = 9.81;     %gravitational acceleration constant (m/s^2)
N = length(Width); 

%Shape dependent calculations
a = zeros(N,1);   %[m^2] area of culvert
WP = zeros(N,1);  %[m] wetted perimeter of culvert


for o = 1:N
   if shape(o)=='Box'
       a(o) = Width(o)*Height(o);
       WP(o) = 2*Width(o)+2*Height(o);
   else
       a(o)=pi()*Width(o)/2*Height(o)/2;
       WP(o) = 2*pi()*sqrt(((Width(o)/2)^2+(Height(o)/2)^2)/2);
   end
end
 
%% Calculation of Type II, Type III and FHA
%Type II
H = zeros(N,1);     
Kc = zeros(N,1);   %[m^-1]
q_II = zeros(N,1); 

H(1) = HW(1)+0.4*Height(1);
Kc(1) = 2*g*n(1)^2/((a(1)/WP(1))^(4/3));
q_II(1) = (a(1)*(sqrt(2*g*H(1))))/(sqrt(1+Ke(1)+Kc(1)*L(1)));

for i = 2:N
    H(i) = L(i)*S(i)+(HW(i))-(0.6*Height(i));
    Kc(i) = 2*g*n(i)^2/((a(i)/WP(i))^(4/3)); 
    q_II(i) = (a(i)*(sqrt(2*g*H(i))))/(sqrt(1+Ke(i)+Kc(i)*L(i)));
end 


%Type III
q_III = zeros(N,1);
H_III = zeros(N,1); 
C = 0.611;  %Vena contracta, constant

for j = 1:N
    H_III(j) = HW(j)+Height(j)/2;
    q_III(j) = a(j)*C*(sqrt(2*g*H_III(j)));
end

%FHA
q_FHA = zeros(N,1);
ks = -0.5; %slope adjestment constant 
H_FHA = zeros(N,1);
for k = 1:N
    H_FHA(k) = HW(k)+Height(k);
    q_FHA(k) = a(k)*(sqrt(Height(k)))*sqrt(((H_FHA(k)/(Height(k)))-Y(k)-(ks*S(k)))/c(k));
end 


