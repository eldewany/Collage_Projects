clc
clear all
%% Givens
        C=1; % Assume unity chord
        c_c=0.04;
        t_c=0.13;
        % Thus
        b=C/4;
        e=t_c/1.3;
        beta=2*c_c;
        a=b*(1+e)/cos(beta);
        x0=-b*e;
        y0=a*beta;
%% Coordinates
        thetad=[0:0.01:2*pi]';
        rd=a;
        xd=rd*cos(thetad);
        yd=rd*sin(thetad);
        x=xd+x0;
        y=yd+y0;
        theta=atan2(y,x);
%% Plotting The airfoil 
        figure
        x1=2*b*cos(theta);
        y1=(2*b*e*(1-cos(theta)).*sin(theta))+(2*b*beta*sin(theta).*sin(theta));
        plot(x1,y1)
        ylim([-0.5,0.5])
        grid on
        xlabel("Airfoil X-Coordinates")
        ylabel("Airfoil Y-Coordinates")
        alpha=[-3*pi/180;6*pi/180;12*pi/180];
%% Plotting velocity distribution, pressure distribution and moment coefficients
        % Iterating for each alpha
        for i=1:length(alpha)
            vtd=-(2*sin(thetad-alpha(i))+2*sin(alpha(i)+beta));
            v=sqrt(vtd.^2./(2-2*cos(2*theta)));
            figure
            plot(x1,v)
            ylim([0,3])
            xlabel("Airfoil X-Coordinates")
            ylabel("Velocity distribution (V/V_inf)")
            title(["Velocity distribution over airfoil at AoA= " ,alpha(i)*180/pi])
            grid on
            cp=1-v.^2;
            figure
            plot(x1,cp)
            ylim([-1,1])
            xlabel("Airfoil X-Coordinates")
            ylabel("Pressure coefficient Cp")
            title(["Pressure distribution over airfoil at AoA= " ,alpha(i)*180/pi])
            grid on
            Cl(i)=2*pi*(1+e)*sin(alpha(i)+beta)
            alphaa=[-3*pi/180:0.01:12*pi/180]';
            Cll=2*pi*(1+e)*sin(alphaa+beta);
            figure
            plot(alphaa,Cll)
            grid on
            xlabel("Angle of attack (alpha)")
            ylabel("Coefficient of lift (Cl)")
            title("Cl Vs. Alpha curve")
            mid=length(thetad)/2;
        % Calculation of moments coeffiecients
          % dividing all the values by 0.5*rho*V_infinity^2 so we can use cp calues
          % directly
            for n=1:length(cp)-1
                fy(n)=(((cp(n)+cp(n+1))/2)*(x1(n+1)-x1(n)));  % Mean force in x1 direction
                fx(n)=-(((cp(n)+cp(n+1))/2)*(y1(n+1)-y1(n))); % Mean force in y1 direction
                Mom(n)= -fy(n)*(x1(n+1)+x1(n))/2 + fx(n)*(y1(n+1)+y1(n))/2 ;% moments due to mean forces at leading edge
            end

            cm(i)=sum(Mom)
            Xcp(i)= cm(i)/Cl(i)
        % Plotting stream lines over airfoil
        
            Vs=100; % let stream velocity be 100 m/s
            gamma = 4*pi*Vs*a*sin(alpha(i)+beta);
            rd1 = linspace(a,10*a,1000);
            [Rd,Td]=meshgrid(rd1,thetad);
        % simplyfing the stream equation    
            psi=gamma./(2*pi).*log(Rd)+Vs./Rd.*(sin(Td-alpha(i))).*(Rd.^2-a.^2);
            Xd = Rd .* cos(Td);
            Yd = Rd .* sin(Td);
            X = x0 + Xd;
            Y = y0 + Yd;
            X1=X.*(1+b^2./(X.^2+Y.^2));
            Y1=Y.*(1-b^2./(X.^2+Y.^2));
            contour(X1,Y1,psi,'b','levelstep',6,'LineWidth',2)
            hold on
            plot(x1,y1,'k','linewidth',2.5)
            xlim([-1 1])
            ylim([-1 1])
            title(['Streamlines at AoA = ',alpha(i)*180/pi],'fontweight','bold')
        end

            fit(alpha,cm','poly1')
            grid on
            xlabel("Angle of attack (alpha)")
            ylabel("Coefficient of moment (Cl)")
            title("Cm Vs. Alpha curve")
