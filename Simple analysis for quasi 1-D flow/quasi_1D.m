clc
clear all
x=linspace(0,1,20);
PbP0=[1,0.98,0.97,0.96085,0.95,0.9,0.8,0.7,0.6,0.5,0.435,0.4,0.3,0.2,0.1,0.04,0];
fun1=@(x) (11.9774.*x.^2-15.4774.*x+6);
A=fun1(x);
plot(x,sqrt(A/pi))
hold on
plot(x,-sqrt(A/pi))
hold off
[Ath,i]=min(A);
Aex=A(end);
fun2=@(x) ((Aex/Ath)-(1/x)*((1+0.2*x^2)/1.2)^3);
Msub=fsolve(fun2,0.1);
Msup=fsolve(fun2,3);
fun3= @(x) ((1+0.2*x^2)^-3.5);
PexP0_sub=fun3(Msub);
PexP0_sup=fun3(Msup);
PexP0_NS=fun3(1);
% % Plots
figure
Min=linspace(0,1,length(PbP0));
plot(PbP0,Min)
Pplot=figure;
xlabel('x')
ylabel('P')
hold on
Mplot=figure;
xlabel('x')
ylabel('M')
hold on
Xplot=figure;
hold on
for i=1:length(PbP0)
    
 if PexP0_sub<PbP0(i)
    Me=fsolve(fun3,0.1); 
    fun4=@(x) ((Aex/x)-(1/Me)*((1+0.2*Me^2)/1.2)^3);
    As=fsolve(fun4,1);
        for i=1:length(A)
            fun5{i}=@(x) ((A(i)./As)-(1./x).*((1+0.2.*x^2)./1.2).^3);
        end
        for k=1:length(A)
            M(k)=fsolve(fun5(k),0.1);
            PP0(k)=fun3(M(k));
        end
elseif PexP0_NS<PbP0(i)<PexP0_sub
    Me=fsolve(fun3,1);
    As=Ath;
        for i=1:length(A)
            fun5{i}=@(x) ((A(i)./As)-(1./x).*((1+0.2.*x^2)./1.2).^3);
        end
        for k=1:length(A)
            M(k)=fsolve(fun5(k),0.1);
            PP0(k)=fun3(M(k));
        end  
     
 end
    figure(Pplot)
    plot(x,PP0)
    figure(Mplot)
    plot(x,M)
end
