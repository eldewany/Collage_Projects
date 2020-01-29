%% Initialization of Variables
clear all
clc
% Initialization of variables
E=zeros;
L=zeros;
I=zeros;
Constr=zeros;
Nodal=zeros;
Mem_Loads=zeros;
Dj=zeros;
%% Input type selection
in=input('Please choose 1 for Vectorized input or 2 for step by step input: ');
if in==1 
%% Vectorized set of inputs
E=[200e9;70e9;70e9];
L=[2;2;4];
I=[8e-4;12e-4;12e-4];
Constr=[1;1;0;0;0;0;1;0];
Nodal=[0;0;-20e3;0;-10e3;0;0;0];
Mem_Loads=[-40e3 0 0;0 -60e3 0;-20e3 30e3 -5e3];
Dj=[1;1;0;0;0;0;-0.01;0];
Nm=3;
Nn=Nm+1;   
%% 
else 
%% Step by step Input
% Input the no. of elements
Nm=input('Please enter no. of memebers: ');
% No. Of nodes
Nn=Nm+1;
% Gathering member information
for i=1:Nm    
E(i,1)=input(['Please enter material of member no. ',num2str(i),':']);
L(i,1)=input(['Please enter Length of member no. ',num2str(i),':']);
I(i,1)=input(['Please enter the MOI of member no. ',num2str(i),':']);
end
% Constructing memeber information matrix
% Forming restrained DOF List
fprintf('Please put 1 for yes and 0 for no \n')
for i=1:Nn
    Constr(2*i-1,1)=input(['Is Translation DOF restrained for node ',num2str(i),':']);
    Constr(2*i,1)=input(['Is Rotation DOF restrained for node ',num2str(i),':']);
end
% Constructing the nodal loads vector
for i=1:Nn
    Nodal(2*i-1,1)=input(['Please input the Shear nodal force at node no. ',num2str(i),':']);
    Nodal(2*i,1)=input(['Please input the Moment nodal force at node no. ',num2str(i),':']);
end
% Constructing the member loads vector
for i=1:Nm
    Mem_Loads(i,1)=input(['Please input the Shear member force on member no. ',num2str(i),':']);
    Mem_Loads(i,2)=input(['Please input the Moment member force on member no. ',num2str(i),':']);
    Mem_Loads(i,3)=input(['Please input the Distributed force on member no. ',num2str(i),':']);
end
fprintf("Please Enter the settlment value if available , if it's restrained Please put 1 , If it's free put 0 \n \n \n")
for i=1:Nn
    Dj((2*i-1),1)=input(['Please input the Axial settlment on node no. ',num2str(i),':']);
    Dj(2*i,1)=input(['Please input the rotational settlment in radian on node no. ',num2str(i),':']);
end
end
%% Processing
% Preparations before solving
Sm=zeros(4,4,Nm);
AE=zeros(2*Nn,1);
M_info=[E , I , L];
RL=Constr;
temp1=find(Dj==0);
temp=find(Dj==1);
Dj(temp1)=NaN;
Dj(temp)=0;
%% Forming the model

% Constructing arbitrary stiffness matrix
%   Zeros is better than sparse in initializing the Sj matrix variable as it's faster in general
    Sj=zeros(2*Nn,2*Nn); 
for i=1:Nm
    DOF=(2*i-1):2*(i+1);
    Em=M_info(i,1);
    Im=M_info(i,2);
    Lm=M_info(i,3);
    Sm(:,:,i)=(Em*Im/Lm^3)*[ 12   6*Lm     -12   6*Lm;
                           6*Lm   4*Lm^2  -6*Lm  2*Lm^2;
                            -12  -6*Lm      12  -6*Lm;
                           6*Lm   2*Lm^2  -6*Lm  4*Lm^2];
    Sj(DOF,DOF)=Sj(DOF,DOF)+Sm(:,:,i);
end
% Constructing arbitrary loads vector
    Aml=zeros;
for i=1:Nm
    DOF=(2*i-1):2*(i+1);
    Pm=Mem_Loads(i,1);
    Mm=Mem_Loads(i,2);
    Wm=Mem_Loads(i,3);
    Lm=M_info(i,3);
    Aml(i,1) =(1.5*Mm/Lm)-(Pm/2)-(Wm*Lm/2);
    Aml(i,2) =(Mm/4)-(Pm*Lm/8)-(Wm*Lm^2/12);
    Aml(i,3) =(-1.5*Mm/Lm)-(Pm/2)-(Wm*Lm/2);
    Aml(i,4) =(Mm/4)+(Pm*Lm/8)+(Wm*Lm^2/12);
    AE(DOF,1)=AE(DOF,1)-Aml(i,:)';
end
% Forming combined loads vector
Ac=Nodal+AE;
% Finding indices of zero and non zero elements in Restrained List
ind=[find(RL==0);find(RL==1)];
% Renumbering Sj,Ac,Dj
Sj_re=Sj(ind,ind);
Ac_re=Ac(ind);
Dj_re=Dj(ind);
% Partitoning Matrices
pos_dof=length(find(RL==0));
S=Sj_re(1:pos_dof,1:pos_dof);
Sdr=Sj_re(1:pos_dof,(1+pos_dof):2*Nn);
Srd=Sj_re((1+pos_dof):2*Nn,1:pos_dof);
Srr=Sj_re((1+pos_dof):2*Nn,(1+pos_dof):2*Nn);
% Partitioninig Ac,Dj
Ad=Ac_re(1:pos_dof);
Arl=-Ac_re((1+pos_dof));
Dr=Dj_re((1+pos_dof):2*Nn);
%% Solving 
% % % % % % Solving % % % % % %
D=inv(S)*(Ad-(Sdr*Dr));
Dj_re=[D;Dr];
Dj(ind)=Dj_re;
Ard=(Srd*D)+(Srr*Dr);
Ar=Arl+Ard;
%% Plotting
% Plotting the graphs
xbar=(0:0.1:1)';
xstart=0;
xbar1=(0:0.1:0.5)';
xbar2=(0.5:0.1:1)';
% Beam Deflection Diagram
def=figure; 
title('Deflection Diagram','interpreter','latex');
xlabel('Beam Length','interpreter','latex');
ylabel('Deflection (m)','interpreter','latex')
hold on
% Figure for Bending Moment Diagram
BMD=figure; 
title('Bending Moment Diagram','interpreter','latex');
xlabel('Beam Length','interpreter','latex');
ylabel('B.M (N.m)','interpreter','latex')
hold on
% Figure for Shear Force Diagram
SFD=figure; 
title('Shear Force Diagram','interpreter','latex');
xlabel('Beam Length','interpreter','latex');
ylabel('S.F (N)','interpreter','latex')
hold on
load=zeros(Nm,3);
id=1:Nm;
load(id,:)=Mem_Loads(:,:);
Am=zeros(Nm,4);
for i=1:Nm
 Lm=M_info(i,3);
 DOF=(2*i-1):2*(i+1);
 % Member End Deflection
 Dm=Dj(DOF,1); 
 % Member End Action
 Am(i,:)=Aml(i,:)+(Sm(:,:,i)*Dm)'; 
 % Deflection
 phi(:,1)=1-3*xbar.^2+2*xbar.^3;
 phi(:,2)=Lm*(xbar-2*xbar.^2+xbar.^3);
 phi(:,3)=3*xbar.^2-2*xbar.^3;
 phi(:,4)=Lm*(-xbar.^2+xbar.^3);
 vm=phi*Dm;
 figure(def)
 plot(xstart+xbar*Lm,vm) 
 % Bending Moment
 MB1=-Am(i,2)+Am(i,1)*xbar1*Lm+0.5*load(i,3)*xbar1.^2*Lm^2;
 MB2=-Am(i,2)+Am(i,1)*xbar2*Lm+0.5*load(i,3)*xbar2.^2*Lm^2+load(i,1)*...
 (xbar2-.5)*Lm-load(i,2);
 figure(BMD)
 plot(xstart+xbar1*Lm,MB1,xstart+xbar2*Lm,MB2)
 % Shear force
 V1=-Am(i,1)-load(i,3)*xbar1*Lm;
 V2=-Am(i,1)-load(i,3)*xbar2*Lm-load(i,1);
 figure(SFD)
 plot(xstart+xbar1*Lm,V1,xstart+xbar2*Lm,V2)
 xstart=xstart+Lm;
end


