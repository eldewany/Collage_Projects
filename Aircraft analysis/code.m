clc;
clear all;
close all;

%% Given

% Airplane Specifications

Cd_0   = 0.0204; 
Cl_0   = 0.244;  % for wing
C_L_max= 1.5; 
Cl_alpha=0.08129;
e      = 0.8  ;
AR     = 6.4    ;
S      = 250;                             % [ft^2]
W      = 6490.8;                            % [1bs ]
We     = 4056.7;                            % [1bs ]
eta_pr = 0.82; 
C_f    = 0.77/(550*3600);   % Specific fuel consumption
P_BHP  = 702;
P_av   = 702*550;                         % [ft.1b/s]
V_st   = 61*1.688;                        % [ft/sec ]

% Take-Off Properties

h_obst=50;  % [ft] transition and climb distance  

% Landing Properties (roscam)

delta_Cd_0_L=0.068;
e_L=0.72;

% Cruise Properties

rho_cr = 0.001265;                  % [slug/ft^3]
V_cr   =250*1.688;                        % [   fps   ]

% Sea Level Properties 

rho_0 =23.77*10^(-4);                     % [slug/ft^3]
P_0   = 14.696*12^2;                      % [1bs/ft^3 ]
T_0   = 518.67;                           % [    R    ]
g   = 32;                               % [  ft/s^2 ]
R     = 1716;
a     = -3.567*10^(-3);                   % [   R/ft  ]
mu    = 0.04;
% Standard Atmosphere Relations

h     = linspace(0,30000,100);            % [ft]
T_h   = T_0+a*h;
P_h   = P_0*(T_h/T_0).^(-g/(a*R));
rho_h = rho_0*(T_h/T_0).^(-1-g/(a*R));

% Velocity Range 

V = linspace(0,1200,100);  %%%                  % [fps]

% Lift coefficient vector
C_L=linspace(0,3,30);
% Landing distance inputs
n_values=[1 3 6];

clmax_l=3.2;
mu_r=0.4;
h_1=6.9; %wing height [ft]
b=40; %wing span [ft]
N=5; %[sec] [Time increment for the free roll (assumed]
th_a=3; % degree (assumed)


%%  Lift Coefficient versus alpha [Cl-alpha]

% Calculations

alpha = -5:0.1:20;
Cl=Cl_0+Cl_alpha.*alpha;   % Cl alpha for thin airofils = 0.11

% Plotting

figure
plot(alpha,Cl,'LineWidth', 2)
grid on
xlabel('\alpha (deg)')
ylabel('Cl')
title('Cl Vs. \alpha');

%%  Drag Coefficient versus alpha [Cd-alpha]

% Calculations
K=1/(pi*e*AR);
Cd=Cd_0+K.*Cl.^2;
xlabel('\alpha (deg)')
ylabel('Cd')

% Plotting
figure 
grid on;
hold on;
plot(alpha,Cd,'LineWidth', 2);
xlabel('\alpha');
ylabel('C_D ');
title(' C_D Vs \alpha ') ;

%%   Drag Coefficient versus Lift Coefficient [Cd-CL]
% Plotting
figure
plot(Cd,Cl,'LineWidth',2)
title('C_l - C_D Polar')
xlabel('C_D')
ylabel('Cl')
grid on

%% Cl^(1/2)/C_d  Versus V 

% Calculation 

Cl_12_Cd     = V.^3*sqrt(2*W/(rho_cr*S))./(Cd_0*V.^4+K*(2*W/(rho_cr*S))^2);
Cl_12_Cd_max = max(Cl_12_Cd);
           i = find(Cl_12_Cd==Cl_12_Cd_max);
           
% Plotting

figure ;
grid on;
hold on;
plot(V,Cl_12_Cd,'LineWidth', 2);
plot(V(i),Cl_12_Cd_max,'o','MarkerSize',10)
xlabel('V_\infty [ft/s]');
ylabel('C_L^1^/^2/C_d ');
title(' C_L^1^/^2/C_d Vs  V_\infty ') ;
legend('C_L^1^/^2/C_d','C_L^1^/^2/C_d _m_a_x');

%%  Cl/C_d  Versus V 

% Calculations 

Cl_Cd=W./((0.5*rho_0.*V.^2*S*Cd_0)+(K*S*(W/S)^2./(0.5*rho_0.*V.^2)));

Cl_Cd_max = max(Cl_Cd);
           r = find(Cl_Cd==Cl_Cd_max);
           
% Plotting

figure ;
grid on;
hold on;
plot(V,Cl_Cd,'LineWidth', 2);
plot(V(r),Cl_Cd_max,'o','MarkerSize',10)
xlabel('V_\infty [ft/s]');
ylabel('C_L / C_d ');
title(' C_L / C_d Vs  V_\infty ') ;
legend('C_L / C_d','C_L / C_d _m_a_x');


%% Cl^(3/2)/C_d  Versus V 

% Calculations 

Cl_32_Cd=sqrt((2*W./(rho_cr.*V.^2.*S)).^3)...
    ./(Cd_0+K.*(2*W./(rho_cr.*V.^2.*S)).^2);
Cl_32_Cd_max = max(Cl_32_Cd);
           m = find(Cl_32_Cd==Cl_32_Cd_max);
           
% Plotting

figure ;
grid on;
hold on;
plot(V,Cl_32_Cd,'LineWidth', 2);
plot(V(m),Cl_32_Cd_max,'o','MarkerSize',10)
xlabel('V_\infty [ft/s]');
ylabel('C_L^3^/^2/C_d ');
title(' C_L^3^/^2/C_d Vs  V_\infty ') ;
legend('C_L^3^/^2/C_d','C_L^3^/^2/C_d _m_a_x');
%%   [ Cl12C_d_V ] & [ ClC_d_V ] & [ Cl32C_d_V ]

% Cl^(1/2)/C_d  Versus V [ Cl12C_d_V ]

figure
hold on
 for n=n_values
cl_12_cd=V.^3*sqrt(2*n*W/(rho_0*S))./(Cd_0*V.^4+K*(2*n*W/(rho_0*S))^2);
plot(V,cl_12_cd,'LineWidth',2)
end
grid on 
title('C_L^1^/^2/C_d - V_\infty at different n')
xlabel('V_\infty [ft/sec]')
ylabel('C_L^1^/^2/C_d')
axis tight
legend('n=1','n=3','n=6')

% Cl/C_d  Versus V [ ClC_d_V ]
figure
hold on
q=0.5*rho_0*V.^2;
 for n=n_values
CL_CD=n*W./((q*S*Cd_0)+(K*S*(n*W/S)^2./q));
plot(V,CL_CD,'LineWidth',2)
 end
grid on 
title('C_L/C_d - V_\infty at different n')
xlabel('V_\infty [ft/sec]')
ylabel('C_L/C_d')
axis tight
legend('n=1','n=3','n=6')

% Cl^(3/2)/C_d  Versus V [ Cl32C_d_V ]

figure
hold on
    for n=n_values
    cl32cd=sqrt((2*n*W./(rho_0.*V.^2.*S)).^3)./(Cd_0+K.*(2*n*W./(rho_0.*V.^2.*S)).^2);
    plot(V,cl32cd,'LineWidth',2)
    end
grid on 
title('C_L^3^/^2/C_d - V_\infty at different n')
xlabel('V_\infty [ft/sec]')
ylabel('C_L^3^/^2/C_d')
axis tight
legend('n=1','n=3','n=6')

%%  Thrust versus Velocity at different Altitudes [ D,T_V(h) ]


% Density values at different different Altitudes(5000,10000,15000,20000)

rho_vector = [2.0482*10^(-3),1.7556*10^(-3),1.5455*10^(-3),0.001265];

% Calculation & Plotting

figure
grid on;
hold on; 
  for n = 1:length(rho_vector)
      
      % Calculation
      
      A1   = 0.5*rho_h(n)*S*Cd_0;
      B1   = 2*K*S*((W/S)^2)/rho_vector(n);
      Tr   = A1.*V.^2+B1./V.^2;  
      T_av = P_av./V;
      
      % Plotting
      
      plot(V,Tr,'LineWidth', 2);

  end
  
    plot(V,T_av,'LineWidth', 2,'LineStyle','--');
    xlabel('V_\infty [ft/s]');
    ylabel('T_R  [1b] ');
    title(' T_R  Vs   V_\infty ');
    set(gca,'Ylim',[0,12000]);
    legend('T_R at h1 = 5000 ft','T_R at h2 = 10000 ft',...
          'T_R at h3 = 15000 ft','T_R at h4 = 20000 ft','T_a_v ');
      
      
%%   Power versus Velocity at different Altitudes [ P-V(h) ]

P1=(0.5*rho_cr.*V.^2*Cd_0.*V*S)/550;   % Power to overcome parasite drag 
P2=((2*K*S*(W/S)^2)./(rho_cr*V))/550;   % Power to overcome induced drag
P=P1+P2;
figure
hold on
plot(V,P1,'--','LineWidth', 2)
plot(V,P2,'--','LineWidth', 2)
plot(V,P,'LineWidth', 2)
plot(V,ones(length(V))*400,'--','LineWidth', 2)
xlabel('Velocity (ft/sec)')
ylabel('Power (HP)')
xlim([0 400])
ylim([0 1000])
grid on
title('Power Vs. Velocity');
legend('P due to Parasite drag','P due to Induced drag','Total Power'...
    ,'Available Power')


%%  Max. Velocity & Min. velcity ( Vmax , Vmin)

% Calculation 

er=@(V_max) V_max-sqrt(((P_av/V_max/S)+sqrt((P_av/V_max/S)^2....
            -4*Cd_0*K*(W/S)^2))/(rho_cr*Cd_0));
V_max=fsolve(er,600);

er=@(V_min) V_min-sqrt(((P_av/V_min/S)-sqrt((P_av/V_min/S)^2....
                     -4*Cd_0*K*(W/S)^2))/(rho_cr*Cd_0));
V_min=fsolve(er,60);

disp(['V_max = ',num2str(V_max),' [ft/sec]']);
disp(['V_min = ',num2str(V_min),' [ft/sec]']);
disp('V_min < V_stall ')
disp('V_min=V_stall')
disp(['V_max = ',num2str(V_max),' [ft/sec]']);
disp(['V_min = ',num2str(V_st),' [ft/sec]']);
%%  V min. power

V_MP = sqrt(2*W*sqrt(1/(3*pi*e*AR*Cd_0))/(rho_cr*S));
disp(['V min. power = ',num2str(V_MP)]);


%% Landing distance

% Inputs
V_f=1.23*V_st; % for n= 1.2
V_td=1.15*V_st;
RR=V_f^2/(0.2*g);
V_inf=0.7*V_td;
h_f=RR*(1-cosd(th_a));
S_a=(50-h_f)/tand(th_a);
S_f=RR*sind(th_a);
T_r=0.4*0.5*rho_0*V_inf^2*S*(Cd_0+delta_Cd_0_L+clmax_l^2*K);
J_t=T_r/W+mu_r;
G=16*(h_1/b)^2/(1+16*(h_1/b)^2);
J_a=rho_0/(2*W/S)*(Cd_0+delta_Cd_0_L+(K+G/(pi*AR*e_L))*clmax_l^2-mu_r*clmax_l);
S_g=N*V_td+1/(2*g*J_a)*log(1+(J_a/J_t)*V_td^2);
disp(' ' )
disp(['Approach distance(S_a) = ',num2str(S_a),' [ft]'])
disp(' ')
disp(['Flare distance(S_f) = ',num2str(S_f),' [ft]'])
disp(' ')
disp(['Ground Roll distance(S_g) = ',num2str(S_g),' [ft]'])
disp(' ')
disp(['Total landing distance = ',num2str(S_a+S_f+S_g),' [ft]'])

%%  Height versus Velocity at q = constant [ h_V (qconst) ]

% Different values for q (q = constant at each vaue) 

q = linspace(50,200,4);

% Calculation & Plotting

figure ;

for n= 1:length(q)
    
% Calculation

v=sqrt(2*q(n)./rho_h);

%  Plotting

plot(v,h,'LineWidth', 2)
hold on;

end

grid on;
xlabel('V_\infty [ft/s]');
ylabel('h  [ft] ');
title(' h  Vs   V_\infty (q = const) ');
legend('q1 = 50','q2 = 100','q3 = 150','q4 = 200');
%%  Height veresus Velocity at E = constant   [ h-V)@(E=const ]

Vm=[0:2:800]';
He=[1000:5000:21000];
hh=He-Vm.^2/(2*g);
figure
plot(Vm,hh,'LineWidth',2)
xlabel('Velocity (ft/sec)')
ylabel('Heigh (Ft)')
legend('He= 1000','He= 6000','He= 11000','He= 16000','He= 21000')
title('Flight envelope , Height Vs. Velocity at constant specific energy');
grid on
%%    Height veresus Velocity at P_s = constant   [ h-V)@(P_s=const ]

%Calculations

H=linspace(0,35000,100);
v=linspace(0,350,100);
aa=-6.5*10^-3;  %[K/m]
[v,H]=meshgrid(v,H);
oc=0:5:40;
T=288.16+aa.*H.*0.305 ;%(1 ft = 0.305 m)   
Rho=rho_0.*(T./288.16).^-(9.81/287/aa+1);
P_req=0.5.*Rho.*S.*Cd_0.*v.^3+2*K*W^2./(Rho.*S.*v);
P_s=(P_av*(Rho./rho_0)-P_req)./W;

% Plotting

figure ;
grid on;
hold on;
contour(v,H,P_s,oc,'LineWidth',2)
title('Altitude-V_\infty')
xlabel('V_\infty [ft/sec]')
ylabel('Altitude [ft]')

%%  Rho(h)
figure
h1=linspace(0,7000,1000); % [m]
[myT, mya, myP, myrho] = atmosisa(h1);
plot(h1*3.37,myrho*0.00194,'LineWidth',2)
grid on 
title('rho - Hight')
xlabel('Hight [ft]')
ylabel('Rho [slug/ft^3]')
axis tight
%%  max ceiling

% Calculation

R_C_max = eta_pr*P_av*(rho_h./rho_0)/W - sqrt(2*sqrt(K/(3*Cd_0))...
          *(W/S)./rho_h).*(1.155*sqrt(4*K*Cd_0));
      
%  Absolute ceiling    
      
rho_abs = fzero(@(rho_abs)eta_pr*P_av*(rho_abs/rho_0)...
                          /W - sqrt(2*sqrt(K/(3*Cd_0))*(W/S)/rho_abs)*...
                           (1.155*sqrt(4*K*Cd_0)), 0.00135);
                       
T_abs   = T_0*(rho_abs/rho_0).^(1/(-1-g/(a*R)));
h_abs   =(T_abs-T_0)/a;
disp(['Absolute ceiling  = ',num2str(h_abs)]);

%  Service ceiling    

rho_serv = fzero(@(rho_serv)  eta_pr*P_av*(rho_serv/rho_0)/W...
                            - sqrt(2*sqrt(K/(3*Cd_0))*(W/S)/rho_serv)...
                              *(1.155*sqrt(4*K*Cd_0))-1.67, 0.00135);
                          
T_serv  = T_0*(rho_serv/rho_0).^(1/(-1-g/(a*R)));
h_serv  =(T_serv-T_0)/a;
disp(['Service ceiling  = ',num2str(h_serv)]);

% Plotting

figure ;
grid on;
hold on;
plot(R_C_max,h,'LineWidth', 2);
line([0,0] , [h(1),h(end)],'LineStyle', '--','LineWidth', 1) ;
line([1.67,1.67] , [h(1),h(end)],'LineStyle', '-.','LineWidth', 1);
legend('R/C_m_a_x','Absolute Ceiling','Service Ceiling')
xlabel('R/C_m_a_x  [ft/s]');
ylabel('h  [ft] ');
title(' R/C_m_a_x  Vs   T ');
%%  Flight envelop 

syms V
kn=0.707;  % index
h=[0:60:12000]';
P_A0=400*550;
P=((0.5*rho_0*V^3*Cd_0*S)+((2*K*S*(W/S)^2)./(rho_0*V)))/550;
n=size(h,1);
Vv=zeros(n,2);
rho_i=0.0024*exp(-6.4714e-5*h);
    for i=1:n  
        sol=real(vpasolve((P_A0*(rho_i(i)/rho_0)^kn)-((5.0220*rho_i(i)*V^3)+(1.9717e4/(rho_i(i)*V)))));
        Vv(i,1)=sol(1);
        Vv(i,2)=sol(2);
    end
    ind=find(Vv(:,1)==Vv(:,2));
figure
Vv(ind(1),2)=Vv(ind(1),2)*-1;
Vv(ind(1),1)=Vv(ind(1),2);
plot(Vv(1:min(ind),:),h(1:min(ind)),'LineWidth',2)
xlabel('Velocity (ft/sec)')
ylabel('Heigh (Ft)')
title('Flight envelope , Height Vs. Velocity');
grid on
%% Rate of climb & Rate of Decent versus V  [ROC/ROD-V]

% Rate of climb 

% Calculation 

v_1=linspace(29.1109,V_max,1000);
D=0.5*rho_0.*v_1.^2*S.*(Cd_0+K.*(2*W./(rho_0.*S.*v_1.^2)).^2);
ROC=(P_av-D.*v_1)./W;

% Plotting 

figure
subplot(2,1,1)
plot(v_1,ROC,'LineWidth',2)
xlabel('V_\infty [ft/sec]')
ylabel('ROC [fpm]')
title('ROC - V_\infty')
grid on
subplot(2,1,2)

% Rate of Descent 

% Calculation 

ROD=D.*v_1./W;

% Plotting 

plot(v_1,ROD,'LineWidth',2)
set(gca, 'YDir','reverse')
xlabel('V_\infty [ft/sec]')
ylabel('ROD [fpm]')
title('ROD - V_\infty')
grid on 

%%  Max climb/decent angle - V
%max. climb angle - V
figure
V = linspace(20,800,100);                  % [fps]

V_theta_max=4*K*(W/S)/(rho_0*eta_pr*(P_av/W));
Theta=eta_pr*P_av./V/W -0.5*rho_0*V.^2*(W/S)^(-1)*Cd_0 -(W/S)*2*K./(rho_0*V.^2);
plot(V,180/pi*asin(Theta),'LineWidth',2)
grid on 
title('max climb angle - V')
xlabel('V_\infty [ft/sec]')
ylabel('Angle [deg]')
axis tight
disp(['V_theta_max  = ',num2str(V_theta_max)]);

%max. decent angle - V

%%   approach angle 

% Calculation
C_L_Landing=2.1; % assumed
C_d_landing=Cd_0+delta_Cd_0_L+K*C_L_Landing^2;
theta_a=asind(1/(C_L_Landing/C_d_landing)-P_av/(W*V_cr));

disp(['Approach Angle  = ',num2str(theta_a)]);
%%   V_n diagram


W_S= W/S;
Cl_max=3.2;
Cl_max_n=-1.2; 
n_max_p=3.6;
n_max_n=-1.4;
V_star=round(sqrt((2*n_max_p)/(rho_0*Cl_max)*W_S));
V=[0:(V_star/20):V_star]';
q=0.5*rho_0*V.^2;
n_p=q*Cl_max/W_S;
n_b=0.5*rho_0*V_star^2*Cl_max/W_S;
V_star_n=round(sqrt((2*n_max_n)/(rho_0*Cl_max_n)*W_S));
V_n=[0:(V_star_n/20):V_star_n]';
q=0.5*rho_0*V_n.^2;
n_n=q*Cl_max_n/W_S;
figure
hold on
plot(V,n_p,'LineWidth',2)
plot([V_star,V_max],[n_max_p,n_max_p],'LineWidth',2)
plot(V_n,n_n,'LineWidth',2)
plot([V_star_n,V_max],[n_max_n,n_max_n],'LineWidth',2)
plot([V_max,V_max],[n_max_p,n_max_n],'LineWidth',2)
xlabel('Velocity (ft/sec)')
ylabel('n')
xlim([0 800])
legend('Loading factor as V Increase','Max +ve Load factor at Cl max','Max negative load factor as V increase','Max -ve load factor at -ve Cl Max')
title('V-N Maneuver Graph');
grid on
hold off
%% Vertical , Horizontal turn rate - V
%vertical turn rate - V
figure
hold on
for n=[1 2 3]
Omega_Vertical=g*(n-1)./V;
plot(V,Omega_Vertical,'LineWidth',2)
end
grid on 
title('vertical turn rate - V')
xlabel('V_\infty [ft/sec]')
ylabel('turn rate []')
axis tight
legend('n=1','n=2','n=3')
%Horizontal turn rate - V
figure
hold on
for n=[1 2 3]
Omega_Horizontal=g*(n^2-1)^0.5./V;
plot(V,Omega_Horizontal,'LineWidth',2)
end
grid on 
title('Horizontal turn rate - V')
xlabel('V_\infty [ft/sec]')
ylabel('turn rate []')
axis tight
legend('n=1','n=2','n=3')
%%   Atmosphere Properties versus Altitude ( Atmosphere(h) )

% Calculations

T_h   = T_0+a*h;
P_h   = P_0*(T_h/T_0).^(-g/(a*R));
rho_h = rho_0*(T_h/T_0).^(-1-g/(a*R));

% Plotting

figure ;
plot(T_h,h,'LineWidth', 2);
grid on;
xlabel('T  [R]');
ylabel('h  [ft] ');
title(' h  Vs   T ');

figure ;
plot(rho_h,h,'LineWidth', 2);
grid on;
xlabel('\rho  [slug/ft\^3]');
ylabel('h  [ft] ');
title(' \rho  Vs   T ');

figure ;
plot(P_h,h,'LineWidth', 2);
grid on;
xlabel('P  [1b/ft/^2]');
ylabel('h  [ft] ');
title(' h  Vs  P  ');




