clear all
clc


%% Variables Initialization
a=0;
b=1;
c=0;
d=1;
h_req=[1/3,0.1,0.01];
for r=1:length(h_req)
    h=h_req(r);
    k=h;
    M=1/h;
    N=1/k;
    tol=1e-8;
 %% Initializing Mesh
    u=zeros(M+1,N+1);
 %% Setting Boundry Condition
    u(M+1,:)=0;
    u(:,N+1)=0;
    for i=1:M+1
        u(1,i)=sin(pi*(i-1)*h);
    end

 %% Initializing interior points
    k1=mean(u(1,:));
    u(2:M,2:N)=k1;
 %% Gauss Seidel Solver
    err=1;
    count=0;
    while err>=tol 
        i=i+1;
        temp=u(2:M,2:N);
        for i=2:M-1    
            for j=2:N-1 
                u(i,j)=0.25*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));                
            end
        end
        u(2,:)=(1-6*h)*u(1,:);
 %% Averaging corner value
        u(1,1)=(u(2,1)+u(1,2))/2;
 %% Error calculations
        err(count+1)=abs((norm(u(2:M,2:N)-temp))/norm(u(2:M,2:N)));
        count=count+1; 
    end
%% Plotting contours
    figure
    contourf(u,50)
    colorbar
    colormap(jet)
    title(['For step k=h=',num2str(h)])
    figure
    plot(1:count,err)
    xlabel('Number of iterations')
    ylabel('Error')
end
