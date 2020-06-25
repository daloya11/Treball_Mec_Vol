%% UNCONTROLLED LATERAL MOTION
clc;
clear all;
    
%% Problem Data

m=288778;   % [kg]
S=511;      % [m^2]
MAC=8.32;   % [m]
b=59.74;    % [m
CD0=0.025; 
h=6069;     % [m]

Ix=2.47e7;  % [kgm2]
Iy=4.49e7;  % [kgm2]
Iz=6.74e7;  % [kgm2]
Ixz=1.32e6; % [kgm2]

Ixp=Ix-Ixz^2/Iz;
Izp=Iz-Ixz^2/Ix;
Ixzp=Ixz/(Ix*Iz-Ixz^2);

mach = 0.65;
a = sqrt(287*294*1.15);
theta0=0;
u0=a*mach;
g=9.81;
dens=1.225*exp(-(g/(287*(273+21)))*(h));
CL0=(2*m*g)/(dens*S*u0^2);
Cw0=(m*g)/(1/2*dens*u0^2*S);

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




  