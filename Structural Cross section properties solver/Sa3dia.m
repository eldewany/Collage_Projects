clear
clc
import Tshape
import Cshape
import Lshape
import Ishape
import Zshape

disp('Welcome to Sa3dia Calculator ')
disp('(1)Solving Fixed Shapes Only')
disp('(2)Solving Custom Shapes only')
disp('(3)Solving wing proplem with Fixed Shapes only')
disp('(4)Solving wing proplem with Custom Shapes only')
disp('(5)Solving wing proplem with both Shapes')
UU=input('Choose your prefered option : ');
x=0;
if UU==1||UU==3||UU==5
    x=input('please enter number of Fixed Shapes : ');
end
A_tot1=zeros(1,x);

if UU==2||UU==4||UU==5
    
op=input('Please enter number of Custom Shapes : ');

 c=x+1;
    A_tot1=zeros(1,(op+x));

end
if UU==1||UU==3||UU==5

Iy=zeros(1,x);
Iz=zeros(1,x);
   
for n=1:x
disp('(1) T Shape.')
disp('(2) C Shape.')
disp('(3) I Shape.')
disp('(4) Z Shape.')
disp('(5) L Shape.')
for ii=1:x
    disp(['For Shape(',num2str(ii),')'])
    disp(['   Iy=',num2str(Iy(ii)/10^6),'e6 mm^4        ,                   Iz=',num2str(Iz(ii)/10^6),'e6 mm^4'])
    disp(['A_total=',num2str(A_tot1(ii)),' mm^2'])
end
p=input('pleas enter your Shape number:');
if p==1
disp('                                                                        ')    
disp('                                          ^Z                             ')
disp('                                          |                             ')
disp('                                   _______|_______                                    ')
disp('                                  |_______|_______|  <-----P.2                           ')
disp('                                        | | |                        ')
disp('                         P.1----->      | | |                         ')
disp('                                     ---|-|-|---------->Y                                     ')
disp('                                        | | |              ')
disp('                                        | | |              ')
disp('                                        |_|_|                               ')
disp('                                                                        ')

 
a=Tshape;
prompt1=' width of part 1 (b1) = ';
prompt2='width of part 2 (b2) = ';
b=[input(prompt1),input(prompt2)];
prompt1='hight of part 1 (h1) = ';
prompt2='hight of part 2 (h2) = ';
h=[input(prompt1),input(prompt2)];
a.b_t=b;
a.h_t=h;
A_tot1(n)=a_t(a);
Iy(n)=iy(a);
Iz(n)=iz(a);
Iyz_z=0;
if Iy(n)>Iz(n)
    alpha1=0;
    I1=Iy(n);
    I2=Iz(n);
else
    alpha1=90;
    I1=Iz(n);
    I2=Iy(n);
end
clc
elseif p==2
    
disp('                                       ^Z                                 ')
disp('                                       |                                 ')
disp('                               ________|_____                                         ')
disp('                              |________|_____|  <-----P.2                                        ')
disp('                              | |      |                                   ')
disp('                              | |      |                           ')
disp('                              | |      |        ')
disp('                 P.1 ----->   | |  ----|------------>Y                                              ')
disp('                              | |      |                                   ')
disp('                              | |      |                                   ')
disp('                              |_|____________                                         ')
disp('                              |______________|                                          ')
disp('                                                                        ')

     
a=Cshape; 
prompt1=' width of part 1 (b1) = ';
prompt2='width of part 2 (b2) = ';
b=[input(prompt1),input(prompt2)];
prompt1='hight of part 1 (h1) = ';
prompt2='hight of part 2 (h2) = ';
h=[input(prompt1),input(prompt2)];
a.b_c=b;
a.h_c=h;
A_tot1(n)=a_c(a);
Iy(n)=iy(a);
Iz(n)=iz(a);
Iyz_z=0;
if Iy(n)>Iz(n)
    alpha1=0;
    I1=Iy(n);
    I2=Iz(n);
    
else
    alpha1=90;
    I1=Iz(n);
    I2=Iy(n);
end
clc
elseif p==3
    
disp('                        ^Z                                                ')
disp('                        |                                                ')
disp('                 _______|_______                                                       ')
disp('                |_______|_______| <------ P.2                                           ')
disp('                      | | |                                              ')
disp('          P.1 ----->  | | |                                              ')
disp('                   ---|-|-|----------->Y                                                     ')
disp('                      | | |                                              ')
disp('                 _____|___|_____                                                      ')
disp('                |_______________|                                                       ')
disp('                                                                        ')
disp('                                                                        ')
disp('                                                                        ')
disp('                                                                        ')


    
    a=Ishape;
prompt1=' width of part 1 (b1) = ';
prompt2='width of part 2 (b2) = ';
b=[input(prompt1),input(prompt2)];
prompt1='hight of part 1 (h1) = ';
prompt2='hight of part 2 (h2) = ';
h=[input(prompt1),input(prompt2)];
a.b_i=b;
a.h_i=h;
A_tot1(n)=a_i(a);
Iy(n)=iy(a);
Iz(n)=iz(a);
Iyz_z=0;
if Iy(n)>Iz(n)
    alpha1=0;
    I1=Iy(n);
    I2=Iz(n);
else
    alpha1=90;
    I1=Iz(n);
    I2=Iy(n);
end
clc
    
elseif p==4
    
disp('                       ^Z                                                 ')
disp('                       |                                                 ')
disp('               ________|_                                                           ')
disp('              |________|_|  <------P.2                                                      ')
disp('                     | | |                                                 ')
disp('       P.1 ----->    | | |               ')
disp('                  ---|-|-|--------->Y                                                      ')
disp('                     |   |                                                 ')
disp('                     |___|______                                                   ')
disp('                     |__________|                                                   ')
disp('                                                                        ')
disp('                                                                        ')
disp('                                                                        ')

     
     a=Zshape;
prompt1=' width of part 1 (b1) = ';
prompt2='width of part 2 (b2) = ';
b=[input(prompt1),input(prompt2)];
prompt1='hight of part 1 (h1) = ';
prompt2='hight of part 2 (h2) = ';
h=[input(prompt1),input(prompt2)];
a.b_z=b;
a.h_z=h;
A_tot1(n)=a_z(a);
Iy(n)=iy(a);
Iz(n)=iz(a);
Iyz_z=Iyz(a);
[I_34,alpha2]=cart2pol(((Iy(n)-Iz(n))./2),(-Iyz_z));
alpha1=alpha2/2;
I1=((Iy(n)+Iz(n))./2)+I_34;
I2=((Iy(n)+Iz(n))./2)-I_34;


    clc
   
elseif p==5
    
disp('                ^Z                  ')
disp('                |                                                        ')
disp('               _|_                                                          ')
disp('              | | |                                                        ')
disp('    P.1---->  | | |                                                        ')
disp('              | | |     ')
disp('            --|-|-|-------->Y                                                            ')
disp('              | | |                                                        ')
disp('              |___|________                                                        ')
disp('              |____________| <----P.2                                                         ')
disp('                                                                        ')
disp('                                                                        ')

 
 a=Lshape;
prompt1=' width of part 1 (b1) = ';
prompt2='width of part 2 (b2) = ';
b=[input(prompt1),input(prompt2)];
prompt1='hight of part 1 (h1) = ';
prompt2='hight of part 2 (h2) = ';
h=[input(prompt1),input(prompt2)];
a.b_l=b;
a.h_l=h;
A_tot1(n)=a_l(a);
Iy(n)=iy(a);
Iz(n)=iz(a);
Iyz_z=Iyz(a);
I1=((Iy(n)+Iz(n))./2)+sqrt(((Iy(n)-Iz(n))./2).^2+(Iyz_z).^2);
I2=((Iy(n)+Iz(n))./2)-sqrt(((Iy(n)-Iz(n))./2).^2+(Iyz_z).^2);
alpha1=atand(-(2.*Iyz_z)./(Iy(n)-Iz(n)))./2;
clc



end
end
clc
for ii=1:x
    disp(['For Shape(',num2str(ii),')'])
    disp(['Iy=',num2str(Iy(ii)/10^6),' e6 mm^4        ,                  Iz=',num2str(Iz(ii)/10^6),' e6 mm^4'])
    disp(['A_total=',num2str(A_tot1(ii)),' mm^2'])
    disp(['Iyz=',num2str(Iyz_z/10^6),' e6 mm^4 '])
    disp(['I1=',num2str(I1/10^6),' e6 mm^4 '])
    disp(['I2=',num2str(I2/10^6),' e6 mm^4 '])
    disp(['Alpha 1=',num2str(alpha1),' degrees '])
    
end

end

if UU==2||UU==4||UU==5
   
 for TT=c:(c+op-1)
   
n=input('Please enter the number of shapes that forms your cross section : ');

% Variables initialIZations
y=zeros(1,n);
z=zeros(1,n);
A=zeros(1,n);
A_t=zeros(1,n);
A_r=zeros(1,n);
A_b=zeros(1,n);
A_l=zeros(1,n);
h=zeros(1,n);
b=zeros(1,n);
IY_t=zeros(1,n);
IY_b=zeros(1,n);
IY_r=zeros(1,n);
IY_l=zeros(1,n);
IZ_t=zeros(1,n);
IZ_b=zeros(1,n);
IZ_r=zeros(1,n);
IZ_l=zeros(1,n);
IY=zeros(1,n);
IZ=zeros(1,n);
qy=zeros(1,n);
qy_t=zeros(1,n);
qy_b=zeros(1,n);
qy_l=zeros(1,n);
qy_r=zeros(1,n);
qz=zeros(1,n);
qz_t=zeros(1,n);
qz_b=zeros(1,n);
qz_l=zeros(1,n);
qz_r=zeros(1,n);
% Main shape

if n==1
    
 disp('For object no. 1')
 prompt='Please enter the width(b) = ';
 b(1)=input(prompt);
 prompt='Please enter the height(h) = ';
 h(1)=input(prompt);
 y(1)=0;
 z(1)=0;
 A(1)=b(1).*h(1);
 IY(1)=b(1).*(h(1).^3)./12;
 IZ(1)=h(1).*(b(1).^3)./12;
 qy(1)=0;
 qz(1)=0;
 

else
    
 disp('For object no. 1')
 prompt='Please enter the width(b) = ';
 b(1)=input(prompt);
 prompt='Please enter the height(h) = ';
 h(1)=input(prompt);
 y(1)=0;
 z(1)=0;
 A(1)=b(1).*h(1);
 IY(1)=b(1).*(h(1).^3)./12;
 IZ(1)=h(1).*(b(1).^3)./12;
 qy(1)=0;
 qz(1)=0;    

 len=n-1;

% Top Shapes
if len>=0
 disp('For Top shapes')
 prompt='Please enter the no. of Top shapes : ';
 i=input(prompt);
 len=len-i;
 i=i+1;
    
    for k=2:i
     
  
    if k==2
        
        disp('For object no. 2')
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        y(k)=0;
        z(k)=(h(1)+h(k))./2;
        A_t(k)=b(k).*h(k);
        IY_t(k)=(b(k).*(h(k).^3)./12)+(A_t(k).*z(k).^2);
        IZ_t(k)=h(k).*(b(k).^3)./12;
        qy_t(k)=A_t(k).*z(k);
        qz_t(k)=A_t(k).*y(k); 
        
     else
        
       disp(['For object no. ',num2str(k)])
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        
        %increment
        if k>2
            inch=0;
          for z=2:k-1 
                
                inch=inch+h(z);
               
          end
       end
        
        y(k)=0;
        z(k)=(h(1)./2) + inch + (h(k)./2);
        A_t(k)=b(k).*h(k);
        IY_t(k)=(b(k).*(h(k).^3)./12)+(A_t(k).*z(k).^2);
        IZ_t(k)=(h(k).*(b(k).^3)./12);
        qy_t(k)=A_t(k).*z(k);
        qz_t(k)=A_t(k).*y(k);
        
    end
    
    end
 end
  
% For Right shapes
if len >0
disp('For right shapes')
prompt='Please enter the no. of Right shapes : ';
i=input(prompt);
len=len-i;
i=i+1;

for k=2:i
    
      if k==2
        
        disp('For object no. 2')
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        y(k)=(b(1)+b(k))./2;
        z(k)=0;
        A_r(k)=b(k).*h(k);
        IY_r(k)=(b(k).*(h(k).^3))./12+(A_r(k).*z(k).^2);
        IZ_r(k)=h(k).*(b(k).^3)./12+(A_r(k)*y(k)^2);
        qy_r(k)=A_r(k).*z(k);
        qz_r(k)=A_r(k).*y(k);        
     else
        
       disp(['For object no. ',num2str(k)])
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        
        %increment
        if k>2
           incb=0;
            for z=2:k-1 
                
                
                incb=incb+b(z);
            end
       end
        
        y(k)=(b(1)./2) + incb + (b(k)./2);
        z(k)=0;
        A_r(k)=b(k).*h(k);
        IY_r(k)=(b(k).*(h(k).^3)./12);
        IZ_r(k)=(h(k).*(b(k).^3)./12)+(A_r(k).*y(k).^2);
        qy_r(k)=A_r(k).*z(k);
        qz_r(k)=A_r(k).*y(k);
      
    
    
      end
end
end


% For Left shapes
if len>0
disp('For left shapes')
prompt='Please enter the no. of Left shapes : ';
i=input(prompt);
len=len-i; 
i=i+1;
for k=2:i
    
      if k==2
        
        disp('For object no. 2')
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        y(k)=-1*((b(1)+b(k))./2);
        z(k)=0;
        A_l(k)=b(k).*h(k);
        IY_l(k)=(b(k).*(h(k).^3))./12+(A_l(k).*z(k).^2);
        IZ_l(k)=h(k).*(b(k).^3)./12+(A_l(k)*y(k)^2);
        qy_l(k)=A_l(k).*z(k);
        qz_l(k)=A_l(k).*y(k);        
     else
        
       disp(['For object no. ',num2str(k)])
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        
        %increment
        if k>2
           incb=0;
            for z=2:k-1 
                
                
                incb=incb+b(z);
            end
       end
        
        y(k)=-1*((b(1)./2) + inch + (b(k)./2));
        z(k)=0;
        A_l(k)=b(k).*h(k);
        IY_l(k)=(b(k).*(h(k).^3)./12);
        IZ_l(k)=(h(k).*(b(k).^3)./12)+(A_l(k).*y(k).^2);
        qy_l(k)=A_l(k).*z(k);
        qz_l(k)=A_l(k).*y(k);
  
      end
end
end

% For Bottom shapes
if len>0
    
disp('For bottom shapes')
prompt='Please enter the no. of Bottom shapes : ';
i=input(prompt);
len=len-i;
i=i+1;
  
for k=2:i
  
    if k==2
        
        disp('For object no. 2')
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        y(k)=0;
        z(k)=(h(1)+h(k))*-1/2;
        A_b(k)=b(k).*h(k);
        IY_b(k)=(b(k).*(h(k).^3)./12)+(A_b(k).*z(k).^2);
        IZ_b(k)=h(k).*(b(k).^3)./12;
        qy_b(k)=A_b(k).*z(k);
        qz_b(k)=A_b(k).*y(k); 
        
     else
        
       disp(['For object no. ',num2str(k)])
        prompt='Please enter the width(b) = ';
        b(k)=input(prompt);
        prompt='Please enter the height(h) = ';
        h(k)=input(prompt);
        
        %increment
        if k>2
            inch=0;
          for z=2:k-1 
                
                inch=inch+h(z);
               
          end
       end
        
        y(k)=0;
        z(k)=-1*((h(1)./2) + inch + (h(k)./2));
        A_b(k)=b(k).*h(k);
        IY_b(k)=(b(k).*(h(k).^3)./12)+(A_b(k).*z(k).^2);
        IZ_b(k)=(h(k).*(b(k).^3)./12);
        qy_b(k)=A_b(k).*z(k);
        qz_b(k)=A_b(k).*y(k);
        
    end
end
end
end


% Values intialIZation
IYz_bar=0;
IZ_bar=0;
qy_bar=0;
qz_bar=0;
A_bar=0;

% Total shape properties
for l=1:n
        
    A_bar=A_bar+A_t(l)+A_b(l)+A_r(l)+A_l(l);
    IYz_bar=IYz_bar+IY_t(l)+IY_l(l)+IY_r(l)+IY_b(l);
    IZ_bar=IZ_bar+IZ_t(l)+IZ_l(l)+IZ_r(l)+IZ_b(l);
    qy_bar=qy_bar+qy_r(l)+qy_l(l)+qy_b(l)+qy_t(l);
    qz_bar=qz_bar+qz_r(l)+qz_l(l)+qz_b(l)+qz_t(l);
end
A_tot1(TT)=A_bar+A(1);
IY_total=IYz_bar+IY(1);
IZ_total=IZ_bar+IZ(1);
qy_total=qy_bar;
qz_total=qz_bar;
y_bar=qz_total./A_tot1(TT);
z_bar=qy_total./A_tot1(TT);
    
%  Principal Axis  

    IY_total_princp=IY_total-A_tot1(TT).*z_bar.^2;
    IZ_total_princp=IZ_total-A_tot1(TT).*y_bar.^2;
    
% Outputs
clc
for ii=1:x
    disp(['For Shape(',num2str(ii),')'])
    disp(['Iy=',num2str(Iy(ii)/10^6),' e6  mm^4        ,                 Iz=',num2str(Iz(ii)/10^6),' e6 mm^4'])
    disp(['A_total=',num2str(A_tot1(ii)),' mm^2'])
    disp('  ')
    disp('  ')
    
end
disp(['For Custom Shape no.(',num2str(TT-x)])
  disp(['center of area position Y=  ',num2str(y_bar),' mm'])
  disp(['center of area position Z=  ',num2str(z_bar),' mm'])
  disp(['Iy =  ',num2str(IY_total_princp./10^6),' e6 mm^4'])
  disp(['Iz =  ',num2str(IZ_total_princp/10^6),' e6 mm^4'])
%   disp(['IY Total=  ',num2str(IY_bar/10^6),'e6 mm^4'])
%   disp(['IZ Total=  ',num2str(IZ_bar/10^6),'e6 mm^4'])
  disp(['A total=  ',num2str(A_tot1(TT)),'  mm^2'])
  disp(['Qy total=  ',num2str(qy_bar),' mm^3'])
  disp(['Qz total=  ',num2str(qz_bar),' mm^3'])
    
  end
    
end
    
if UU==3||UU==4||UU==5
    disp('   For the shown wing: ')
    disp(' *From double symmetry ; ')
    disp(' - centroidal Y & Z axes are identified ')
    disp('   Iyz=0')
    disp(' - only quarter of the cross-section can be considered')
    disp('   Note That : Solve using (IdealIzation of Thin-Walled Cross-section).')
    disp('  ')
disp('                                                                        ')
disp('                                    ^Z                                    ')
disp('                                    |                                    ')
disp('               __________ __________|__________ ___________                                                                      ')
disp('              |_O________O_________(|)_________O_________O_|                                                        ')
disp('              | |  ^                |                    | |                 ')
disp('              | |  |                |                    | |                 ')
disp('              | |  t2               |                    | |                 ')
disp('              | |                ---|--------------------|-|---------->Y                                 ')
disp('              | |                   |                    | |                 ')
disp('        t1--->| |<---                                    | |                 ')
disp('              |_|________ __________ __________ _________|_|                                                                       ')
disp('              |_O________O__________O__________O_________O_|                                                       ')
disp('                                    ^          ^         ^                 ')
disp('                                    |          |         |                ')
disp('                                    |          |         |                 ')
disp('                                    I          II       III                 ')
disp('                                                                        ')


h=input('Please Enter wing hight = ');
b=input('Please Enter wing width = ');
t1=input('Please Enter thickness 1 = ');
t2=input('Please Enter thickness 2 = ');

disp('Please Select The Flange number')
l=input('For Flange I choose Shape number:');
A_1=A_tot1(l);
ll=input('For Flange II choose Shape number:');
A_2=A_tot1(ll);
L=input('For Flange III choose Shape number:');
A_3=A_tot1(L);
Aw=(b-t1)*t2;
Ah=(h-t2)*t1;

a=zeros(1,6);
y=zeros(1,6);
z=zeros(1,6);
for k=1:5
    y(k)=(k-1)*(b-t1)/8;
    z(k)=(h-t2)/2;
end
y(6)=4*(b-t1)/8;
a(1)=(Aw/24)+A_1/2;
a(2)=(Aw/6);
a(3)=(Aw/12)+A_2;
a(4)=(Aw/6);
a(5)=(Aw/24)+A_3+(Ah/6);
a(6)=(1*Ah/3);
Iy=0;
Iz=0;
at=0;
for k=1:6
    Iy_S=a(k)*z(k)^2;
    Iz_S=a(k)*y(k)^2;
    Iy=Iy+Iy_S;
    Iz=Iz+Iz_S;
    at=at+a(k);
end
disp(['Iy = ',num2str(4*Iy/10^6),' e6  mm^4'])
disp(['Iz = ',num2str(4*Iz/10^6),' e6  mm^4'])
disp(['a=',num2str(at)])
end
disp('(1) End the programe ')
disp('(2) Resolve another Cross-section')
kk=input(' your choice is: ');
if kk==1
    disp('Thanks to using Sa3dia Calculator (^_^)')
elseif kk==2
    Sa3dia
else 
    clear
    clc
    disp('You are joking so we end the Programe')
end
    