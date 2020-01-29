clc
clear
close all
% Mohamed Hassan Hosny
% Sec:2 B.N:12

%% inputs
% Input materials 
disp('Note that all dimension are SI units')
N_r=input('Number of Elements = ');
disp('Please enter number of materials')
N_m=input('Number of materials = ');
E=zeros(N_m,1);

for i=1:N_m
    E(i)=input(['E ',num2str(i),'=']);
end
%Input property 
disp('Please enter number of Properties')
N_p=input('Number of Areas = ');
Ar=zeros(N_p,1);
for i=1:N_p
    Ar(i)=input(['Area ',num2str(i),'= ']);
end
t=zeros(N_r,1);
s=zeros(N_r,1);

if N_m==1 && N_p==1
    t=ones(N_r,1);
    s=ones(N_r,1);
else
    for i=1:N_r
    
    t(i)=input([' Element ',num2str(i),' material num = ']);
    s(i)=input([' Element ',num2str(i), ' property num = ']);
    
end
end
M_info=[ t s ];

N_n=input('Number of nodes = ');
x_co=zeros(N_n,1);
y_co=zeros(N_n,1);

for i=1:N_n
disp(['for Node num (',num2str(i),')'])
x_co(i)=input(' X coordinates = ');
y_co(i)=input(' Y coordinates = ');
end

N_co=[x_co , y_co];
f_ma=zeros(N_r,1);
l_ma=zeros(N_r,1);

for i=1:N_r

    disp('Please enter the elements start and end nodes')
    disp(['for element num (',num2str(i),')'])
    f_ma(i)=input('First Node of the element number = ');
    l_ma(i)=input('First Node of the element number = ');
end 

Conect=[f_ma,l_ma];

v=zeros(N_n,1);
u=zeros(N_n,1);

disp('Please enter DOF')
for i=1:N_n
    disp('Note that if DOF is Free Please enter Nan')
    u(i)=input(['horizontal displacement for node ',num2str(i),' = ']);
    v(i)=input(['vertical displacement for node ',num2str(i),' = ']);
 
end

constr=[u,v];

F_x=zeros(N_n,1);
F_y=zeros(N_n,1);
for i=1:N_n
    F_x(i)=input(['Fx on node ',num2str(i),' = ']);
    F_y(i)=input(['Fy on node ',num2str(i),' = ']);
end

Px=zeros(N_r,1);
Py=zeros(N_r,1);
M=zeros(N_r,1);
wx=zeros(N_r,1);
wy=zeros(N_r,1);

disp('1- Forces on elemnts')
disp('2- No Forces on elements')
h=input('Please chose on of the above : ');
if h==1
for i=1:N_r
    disp(['For element num (',num2str(i),')'])
    Px(i)=input('Force in x direction = ');
    Py(i)=input('Force in y direction = ');
    M(i)=input(' element Moment load = ');
    wx(i)=input('Element distributed  load in x direction = ');
    wy(i)=input('Element distributed  load in y direction = ');
end
end



num_node=1:1:N_n;
num_ele=1:1:N_r ;
N_load=[num_node' F_x F_y];
M_load=[num_ele' M Px Py wx wy];

%% Stiffness Matrix
D_j(:,1)=reshape(constr',[2*N_n,1]);
R_L=1-isnan(D_j);
S_j=sparse(2*N_n,2*N_n);

for i=1:N_r
    N1=Conect(i,1);
    N2=Conect(i,2);
    dx=N_co(N2,1)-N_co(N1,1);
    dy=N_co(N2,2)-N_co(N1,2);
    L_m=sqrt(dx^2+dy^2);
    cx=dx/L_m;
    cy=dy/L_m;
    E_m=E(M_info(i,1));
    A_me=Ar(M_info(i,2));
    Dof=[2*N1-1 ,2*N1 ,2*N2-1 , 2*N2];
    
    S_md(:,:,i)=E_m*A_me/L_m* [  cx^2    cx*cy   -cx^2    -cx*cy 
                                 cx*cy    cy^2   -cx*cy   -cy^2  
                                -cx^2    -cx*cy   cx^2    cx*cy  
                                -cx*cy   -cy^2   cx*cy     cy^2  ];
    
    S_j(Dof,Dof)=S_j(Dof,Dof)+S_md(:,:,i);
    
    RT(:,:,i)=[cx  cy   0   0
               -cy cx   0   0 
                0   0   cx  cy
                0   0  -cy  cx ];
    L(i)=L_m;
end    
%% Forces
A_E=zeros(2*N_n,1);
A_ml=zeros(N_r,4);
A_mlb=A_ml;
for i=1:size(M_load,1)
    
    id=M_load(i,1);
    M_m=M_load(i,2);
    P_xm=M_load(i,3);
    P_ym=M_load(i,4);
    w_xm=M_load(i,5);
    w_ym=M_load(i,6);
    N1=Conect(id,1);
    N2=Conect(id,2);
    L_m=M_info(id);
    Dof=[2*N1-1 , 2*N1 , 2*N2-1 , 2*N2];
         
    A_ml(id,1)=-P_xm/2-w_xm*L_m/2;
    A_ml(id,2)=M_m/L_m-P_ym/2-w_ym*L_m/2;
     
    A_ml(id,3)=A_ml(i,1);
    A_ml(id,4)=-M_m/L_m-P_ym/2-w_ym*L_m/2;
    
    A_mlb(id,:)=transpose(RT(:,:,id)*transpose(A_ml(id,:)));
     
    A_E(Dof,1)=A_E(Dof,1)-A_mlb(id,:)';
end

%% Nodal load 

A=zeros(2*N_n,1);

for i=1:size(N_load,1)
    
    id=N_load(i,1);
    A(2*id-1:2*id,1)=N_load(i,2:3);
end 
A_c=A+A_E;

%% Renumbering DOF
ind=[find(R_L==0);find(R_L==1)];
Re(ind,1)=(1:2*N_n)';
ndpos=length(find(R_L==0));

%Renumbering and Partioning

S_jre(Re,Re)=S_j;
S=S_jre(1:ndpos,1:ndpos);
S_dr=S_jre(1:ndpos,ndpos+1:2*N_n);
S_rd=transpose(S_dr);
S_rr=S_jre(ndpos+1:2*N_n,ndpos+1:2*N_n);

A_cre(Re,1)=A_c;
A_d=A_cre(1:ndpos);
A_rl=-A_cre(ndpos+1:2*N_n);
D_jre(Re,1)=D_j;
D_r=D_jre(ndpos+1:2*N_n);

%% Solution 
D=S^-1*(A_d-S_dr*D_r);
D_jre=[D;D_r];
D_j=D_jre(Re);
A_rd=S_rd*D + S_rr*D_r;
A_R=A_rl+A_rd;

%% Post-Processing
DSF=2;
figure
hold on

for i=1:N_r
    N1=Conect(i,1);
    N2=Conect(i,2);
    De_fx=D_j([2*N1-1,2*N2-1],1);
    De_fy=D_j([2*N1,2*N2],1) ;
    
    X_vec=N_co([N1,N2],1);
    Y_vec=N_co([N1,N2],2);
    plot(X_vec,Y_vec,'--')  %Original Structure
    plot(X_vec+DSF*De_fx,Y_vec+DSF*De_fy)
end


xbar1=(0:0.1:0.5)';
xbar2=xbar1+0.5;
Load=zeros(N_r,5);
id=M_load(:,1);
Load(id,:)=M_load(:,2:6);

for i=1:N_r
    
    N1=Conect(i,1);
    N2=Conect(i,2);
    Dof=[2*N1-1 , 2*N1 , 2*N2-1 , 2*N2];
    Dno=D_j(Dof,1);   Lm=L(i);
    A_m(i,:)=A_ml(i,:)+(RT(:,:,i)*S_md(:,:,i)*Dno)';
    
    BM1=-A_m(i,2)*xbar1*L_m-0.5*Load(i,5)*xbar1.^2*L_m^2;
    BM2=-A_m(i,2)*xbar2*L_m-0.5*Load(i,5)*xbar2.^2*L_m^2-Load(i,3)*(xbar2-0.5)*L_m+Load(i,1);
    SF1=-A_m(i,2)-Load(i,5)*xbar1*L_m;
    SF2=-A_m(i,2)-Load(i,5)*xbar2*L_m-Load(i,3);
    TF1=-A_m(i,1)-Load(i,4)*xbar1*L_m;
    TF2=-A_m(i,1)-Load(i,4)*xbar2*L_m-Load(i,2);
    
    figure
    subplot(3,1,1)
    plot(xbar1,TF1,xbar2,TF2)
    xlabel('x_m')
    title(['Thrust diagram of member',num2str(i)])
    subplot(3,1,2)
    plot(xbar1,SF1,xbar2,SF2)
    xlabel('x_m')
    title(['Shear diagram of member',num2str(i)])
    subplot(3,1,3)
    plot(xbar1,BM1,xbar2,BM2)
    xlabel('x_m')
    title(['Bending diagram of member',num2str(i)])
end


