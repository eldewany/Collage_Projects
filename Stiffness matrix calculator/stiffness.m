E=2*10^7;
I=0;
A=3;
L=483.7355;
cx=0.8682;
cy=0.4961;
format long g
k11=E*(((A/L)*cx^2)+(((12*I)/L^3)*cy^2));
k21=E*(((A/L)-((12*I)/L^3))*cx*cy);
k22=E*(((A/L)*cy^2)+(((12*I)/L^3)*cx^2));
k31=E*((-6*I*cy)/L^2);
k32=E*((6*I*cx)/L^2);
k33=E*(4*I/L);
mat=[k11 k21 k31 -k11 -k21 k31;
    k21 k22 k32 -k21 -k22 k32;
    k31 k32 k33 -k31 -k32 k33/2;
    -k11 -k21 -k31 k11 k21 -k31;
    -k21 -k22 -k32 k21 k22 -k32;
    k31 k32 k33/2 -k31 -k32 k33]
