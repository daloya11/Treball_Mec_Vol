%% UNCONTROLLED LATERAL MOTION

%% PART 1
%Compute the eigenvalues and eigenvectors

%Lateral stability derivatives

%               Cy     Cl     Cn
coef_lateral=[-0.90 -0.160  0.160;   %beta
                0   -0.340 -0.026;   %p
                0    0.130 -0.280;]; %r
    
    
                %            Y                                   L                                  N
comp_lateral=[0.5*dens*u0*S*coef_lateral(1,1)    0.5*dens*u0*S*b*coef_lateral(1,2)    0.5*dens*u0*S*b*coef_lateral(1,3);         %v
              0.25*dens*u0*S*b*coef_lateral(2,1) 0.25*dens*u0*S*b^2*coef_lateral(2,2) 0.25*dens*u0*S*b^2*coef_lateral(2,3);       %p
              0.25*dens*u0*S*b*coef_lateral(3,1) 0.25*dens*u0*S*b^2*coef_lateral(3,2) 0.25*dens*u0*S*b^2*coef_lateral(3,3)];      %r
  
  
A=zeros(4);
%Fila 1
A(1,1)=comp_lateral(1,1)/m;
A(1,2)=comp_lateral(2,1)/m;
A(1,3)=(comp_lateral(3,1)/m)-u0;
A(1,4)=g;

%Fila 2
A(2,1)=comp_lateral(1,2)/Ixp+Ixzp*comp_lateral(1,3);
A(2,2)=comp_lateral(2,2)/Ixp+Ixzp*comp_lateral(2,3);
A(2,3)=comp_lateral(3,2)/Ixp+Ixzp*comp_lateral(3,3);

%Fila 3
A(3,1)=Ixzp*comp_lateral(1,2)+comp_lateral(1,3)/Izp;
A(3,2)=Ixzp*comp_lateral(2,2)+comp_lateral(2,3)/Izp;
A(3,3)=Ixzp*comp_lateral(3,2)+comp_lateral(3,3)/Izp;

%Fila 4
A(4,2)=1;

[eigenvectors_lateral,eigenvalues_lateral]=eig(A);
eigenvalues_lateral=diag(eigenvalues_lateral);

%% PART 2 
%Plot the eigenvalues for diffrent Cl_beta
figure;
hold on;
grid on;
xlabel("Re");
ylabel("Im");
Clbeta = zeros(1,5);
for i=1:5
    
Clbeta(i)=-0.375 + 0.075*i;
    
%               Cy      Cl        Cn
coef_lateral=[-0.90  Clbeta(i)   0.160;   %beta
                0     -0.340    -0.026;   %p
                0      0.130    -0.280;]; %r
    
                %            Y                                   L                                  N
comp_lateral=[0.5*dens*u0*S*coef_lateral(1,1)    0.5*dens*u0*S*b*coef_lateral(1,2)    0.5*dens*u0*S*b*coef_lateral(1,3);         %v
              0.25*dens*u0*S*b*coef_lateral(2,1) 0.25*dens*u0*S*b^2*coef_lateral(2,2) 0.25*dens*u0*S*b^2*coef_lateral(2,3);       %p
              0.25*dens*u0*S*b*coef_lateral(3,1) 0.25*dens*u0*S*b^2*coef_lateral(3,2) 0.25*dens*u0*S*b^2*coef_lateral(3,3)];      %r
      
    
  A_Clbeta=zeros(4);
%Fila 1
A_Clbeta(1,1)=comp_lateral(1,1)/m;
A_Clbeta(1,2)=comp_lateral(2,1)/m;
A_Clbeta(1,3)=(comp_lateral(3,1)/m)-u0;
A_Clbeta(1,4)=g;

%Fila 2
A_Clbeta(2,1)=comp_lateral(1,2)/Ixp+Ixzp*comp_lateral(1,3);
A_Clbeta(2,2)=comp_lateral(2,2)/Ixp+Ixzp*comp_lateral(2,3);
A_Clbeta(2,3)=comp_lateral(3,2)/Ixp+Ixzp*comp_lateral(3,3);

%Fila 3
A_Clbeta(3,1)=Ixzp*comp_lateral(1,2)+comp_lateral(1,3)/Izp;
A_Clbeta(3,2)=Ixzp*comp_lateral(2,2)+comp_lateral(2,3)/Izp;
A_Clbeta(3,3)=Ixzp*comp_lateral(3,2)+comp_lateral(3,3)/Izp;

%Fila 4
A_Clbeta(4,2)=1;

[eigenvectors_lat_Clbeta,eigenvalues_lat_Clbeta]=eig(A_Clbeta);
eigenvalues_lat_Clbeta=diag(eigenvalues_lat_Clbeta);    
    
Re=real(eigenvalues_lat_Clbeta);
Im=imag(eigenvalues_lat_Clbeta);
scatter(Re,Im,25,'filled') ;   
        
end

title('\fontsize{16}Eigenvalues evolution for different C_l_\beta')
legend('C_l_beta = -0.3','C_l_beta = -0.225','C_l_beta = -0.15','C_l_beta = -0.075','C_l_beta = -0');

%% PART 3
% Take the value of Cl_beta that gives a slightly unstable spiral mode, and
% integrate numerically and plot this mode corresponding to an initial
% perturbed mank angle phi=5º

tf=200; %[s]
Nstep=201;
tspan=linspace(0,tf,Nstep);
options=odeset("RelTol",1e-12,"AbsTol",1e-12);

xi=[0 0 0 0];

[t,x]=ode45(@section3_2,tspan,xi,options);
figure
hold on
grid on
plot(t(:),x(:,1)/u0)  %x(:,1)/u0 és aproximadament beta
xlabel('t [s]')
ylabel('\beta [rad]')

figure
hold on
grid on
plot(t(:),x(:,2))
xlabel('t [s]')
ylabel('p [rad/s]')

figure
hold on
grid on
plot(t(:),x(:,3))
xlabel('t [s]')
ylabel('r [rad/s]')

figure
hold on
grid on
plot(t(:),x(:,4))
xlabel('t [s]')
ylabel('\phi [rad]')


xhi=-x(:,1)/u0; %xhi=-beta
xe=u0.*t;
ye=u0.*xhi.*t+x(:,1);

figure
hold on
grid on
plot(xe,ye)
xlabel('x_E')
ylabel('y_E')
title('Trajectoria')



  