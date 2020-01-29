clc
clear all
% %  Assumptions : 
% % 1- The diaphram is at the middle of the tube
% % 2- Some of the equations are specific only for air (gamma= 1.4)
p4=input('Please enter P4: \n');
p1=input('Please enter P1: \n');
gamma1=input('Please enter gamma for the first gas: \n');
gamma2=input('Please enter gamma for the second gas: \n');
t1=input('Please enter T1: \n');
t4=input('Please enter T4: \n');
l=input('Please enter Tube length: \n');
mw1=input('Please enter Molecular weight of gas in region 1: \n');
mw4=input('Please enter Molecular weight of gas in region 4: \n');
rho1=input('Please enter density of gas in region 1: \n');
rho4=input('Please enter density of gas in region 4: \n');
r1=8314/mw1;
r4=8314/mw4;
a1=sqrt(gamma1*r1*t1);
a4=sqrt(gamma1*r4*t4);
inp=p4/p1;
fun1=@(x) (inp-x*(1-(((gamma1-1)*(a1/a4)*(x-1))/(sqrt(2*gamma1*((2*gamma1)+((gamma1+1)*(x-1)))))))^(-2*gamma1/(gamma1-1)));
p2p1=fsolve(fun1,1);
p2=p2p1*p1;
t2t1=p2p1*((((gamma1+1)/(gamma1-1))+p2p1)/(1+((gamma1+1)/(gamma1-1))*p2p1));
t2=t2t1*t1;
a2=sqrt(gamma1*r1*t2);
rho2rho1=(1+((gamma1+1)/(gamma1-1))*p2p1)/(((gamma1+1)/(gamma1-1))+p2p1);
rho2=rho2rho1*rho1;
% be calculated
Ws=a1*sqrt(((gamma1+1)/(2*gamma1))*(p2p1-1)+1);
Ms=Ws/a1;
Up=(a1/gamma1)*(p2p1-1)*sqrt((2*gamma1/(gamma1+1))/(p2p1+((gamma1-1)/(gamma1+1))));
fun2=@(x) ((x/(x^2-1))-(Ms/(Ms^2-1))*sqrt(1+((2*((gamma2-1)/(gamma2+1)))*(Ms^2-1)*(gamma2+((Ms^2+1)/Ms^2)))));
Mr=fsolve(fun2,3);
Wr=(Mr*a2)-Up;
% % For gamma = 1.4 , N.S Relations, after reflection
p5p2=(7*Mr^2-1)/6;
rho5rho2=(6*Mr^2)/(Mr^2+5);
t5t2=p5p2/rho5rho2;
% %  before expansion waves
p3p4=p2p1*(p1/p4);
rho3rho4=p3p4^(1/gamma2);
rho3=rho3rho4*rho4;
t3t4=p3p4^((gamma1-1)/gamma2);
t3=t3t4*t4;
a3=sqrt(gamma2*r4*t3);
% time for the sw to reach the the tube end
T1=(l/2)/Ws;
disp(['The time limits of the simulation is 0 and ',num2str(T1)])
dp=T1*Up;
dr=(l/2)-dp;
fun3= @(x) ((x/Up)-((l-x)/Wr));
T2=fsolve(fun3,1)/Up;
t=input('Please enter the time desired: \n');
fun4= @(x) ((2/(gamma2+1))*(a4+(x/t)));
x_3=fsolve(fun4,l/2);
d1=l/2-(Ws*t);
l_1=linspace(Ws*t,l/2,10);
d2=d1-(Up*t);
l_2=linspace(Up*t,Ws*t,10);
d3=d2+x_3;
l_31=linspace(0,Up*t,10);
l_32=linspace(0,x_3,10);
l_4=linspace(-l/2,x_3,10);
p_1=p1*ones(length(l_1));
p_2=p2*ones(length(l_2));
p_31=p2*ones(length(l_31));
u=linspace(Up,0,10);
p_32=p4*((1-((gamma2-1)/2).*(u./a4)).^((2*gamma2)/(gamma2-1)));
t_32=t4*((1-((gamma2-1)/2).*(u./a4)).^2);
rho_32=rho4*((1-((gamma2-1)/2).*(u./a4)).^(1/(gamma2-1)));
p_4=p4*ones(length(l_4));
% % plotting
figure
plot(l_1,p_1,'r')
ylim([0 10])
xlim([-l/2 l/2])
hold on
plot(l_2,p_2,'g')
plot(l_31,p_31,'b')
plot(l_32,p_32,'b')
plot(l_4,p_4,'k')
ylim([0 inf])
xlim([-l/2 l/2])
xlabel('X [m]')
ylabel('P(x) ')
figure
plot(l_1,0,'r')
hold on
plot(l_2,ones(length(l_2))*Up,'g')
plot(l_31,ones(length(l_31))*Up,'k')
plot(l_32,u,'k')
xlabel('X [m]')
ylabel('U(x) [m/s]')
xlim([-l/2 l/2])
figure
plot(l_1,ones(length(l_2))*t1,'r')
hold on
plot(l_2,ones(length(l_2))*t2,'g')
plot(l_31,ones(length(l_31))*t3,'b')
plot(l_32,t_32,'b')
plot(l_4,ones(length(l_31))*t4,'k')
xlabel('X [m]')
ylabel('T(x) [Kelvin]')
xlim([-l/2 l/2])
figure
plot(l_1,rho1*ones(length(l_1)),'r')
ylim([0 10])
xlim([-l/2 l/2])
hold on
plot(l_2,rho2*ones(length(l_2)),'g')
plot(l_31,rho3*ones(length(l_31)),'b')
plot(l_32,rho_32,'b')
plot(l_4,rho4*ones(length(l_4)),'k')
ylim([0 inf])
xlim([-l/2 l/2])
xlabel('X [m]')
ylabel('Density(x) ')


