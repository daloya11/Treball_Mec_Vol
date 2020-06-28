%Treball mecanica de vol
clear all 
clc
%% Problem Data

m=288778;   % [kg]
S=511;      % [m^2]
MAC=8.32;   % [m]
b=59.74;    % [m]
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

%% Uncontrolled Longitudinal Motion
%% QUESTION 1: 
%Find the eigenvalues and their corresponding eigenvectors

%      CT       CD      CL      Cm
data_1=[-0.055   0       0.13    0.013;  %û
        0        0.2     4.4     -1.0;   %alpha
        0        0       5.65    -20.5;  %^q
        0        0       6.6     -4.0];  %alpha_prima

Cxu = data_1(1,1)-data_1(1,2);  % Cxu = CTu - CDu
Czu = -data_1(1,3);             % Czu = -CLu
Cxa = CL0-data_1(2,2);          % Cxa = CL0 - CDa
Cza = -data_1(2,3)-CD0;         % Cza = -CLa - CD0
Cxq = -data_1(3,2);             % Cxq = -CDq
Czq = -data_1(3,3);             % Czq = -CLq
Cxap = -data_1(4,2);            % Cxap = -CDap
Czap = -data_1(4,3);            % Czap = -CLap


%           X                                 Z                                M
comp=[0.5*dens*u0*S*Cxu         -dens*u0*S*Cw0+0.5*dens*u0*S*Czu    0.5*dens*u0*S*MAC*data_1(1,4);      %u
      0.5*dens*u0*S*Cxa         0.5*dens*u0*S*Cza                   0.5*dens*u0*S*MAC*data_1(2,4);      %w
      0.25*dens*u0*S*MAC*Cxq    0.25*dens*u0*S*MAC*Czq              0.25*dens*u0*S*MAC^2*data_1(3,4);   %q
      0.25*dens*S*MAC*Cxap      0.25*dens*S*MAC*Czap                0.25*dens*S*MAC^2*data_1(4,4)];     %w·

A=zeros(4);

%Fila Increment u
A(1,1)=comp(1,1)/m;
A(1,2)=comp(2,1)/m;
A(1,4)=-g;

%Fila w·
A(2,1)=comp(1,2)/(m-comp(4,2));
A(2,2)=comp(2,2)/(m-comp(4,2));
A(2,3)=(comp(3,2)+m*u0)/(m-comp(4,2));

%Fila q·
A(3,1)=(1/Iy)*(comp(1,3)+comp(4,3)*comp(1,2)/(m-comp(4,2)));
A(3,2)=(1/Iy)*(comp(2,3)+comp(4,3)*comp(2,2)/(m-comp(4,2)));
A(3,3)=(1/Iy)*(comp(3,3)+comp(4,3)*(comp(3,2)+m*u0)/(m-comp(4,2)));

%Fila Increment Theta·
A(4,3)=1;

% Short Period Approximation
Ashort=zeros(3);

Ashort(1,1)=A(2,2);
Ashort(1,2)=A(2,3);
Ashort(2,1)=A(3,2);
Ashort(2,2)=A(3,3);
Ashort(3,2)=A(4,3);

% Fugoid Approximation
Afug=zeros(3);

Afug(1,1)=A(1,1)-(A(3,1)/A(3,3))*A(1,3);
Afug(1,2)=A(1,2)-(A(3,2)/A(3,3))*A(1,3);
Afug(1,3)=A(1,4);
Afug(2,1)=A(2,1)-(A(3,1)/A(3,3))*A(2,3);
Afug(2,2)=A(2,2)-(A(3,2)/A(3,3))*A(2,3);
Afug(3,1)=-A(3,1)/A(3,3);
Afug(3,2)=-A(3,2)/A(3,3);

[eigenvectors,eigenvalues]=eig(A);
eigenvalues=diag(eigenvalues);
[eigenvectors_short,eigenvalues_short]=eig(Ashort);
eigenvalues_short=diag(eigenvalues_short);
[eigenvectors_fug,eigenvalues_fug]=eig(Afug);
eigenvalues_fug=diag(eigenvalues_fug);

% Plot Q1
% Convert to polar in the report
disp('Question 1 results:');
display(eigenvalues);
display(eigenvalues_short);
display(eigenvalues_fug);


%% QUESTION 2 
% Plot in the complex plane the locus of the eigenvalues when the 
% stability derivative Cmq varies within an interval +- 50% about the
% nominal value, but the rest of stability derivatives remain unchanged

figure;
hold on;
grid on;
xlabel("Re");
ylabel("Im");
%axis([-0.9 0 -0.65 0.65]);
for i =1:11
    %      CT       CD      CL      Cm
    data_1=[-0.055   0       0.13    0.013;  %û
            0        0.2     4.4     -1.0;   %alpha
            0        0       5.65    -10.25-2.05*(i-1);  %^q
            0        0       6.6     -4.0];  %alpha_prima

    Cxu = data_1(1,1)-data_1(1,2);  % Cxu = CTu - CDu
    Czu = -data_1(1,3);             % Czu = -CLu
    Cxa = CL0-data_1(2,2);          % Cxa = CL0 - CDa
    Cza = -data_1(2,3)-CD0;         % Cza = -CLa - CD0
    Cxq = -data_1(3,2);             % Cxq = -CDq
    Czq = -data_1(3,3);             % Czq = -CLq
    Cxap = -data_1(4,2);            % Cxap = -CDap
    Czap = -data_1(4,3);            % Czap = -CLap
    
    %           X                       Z                           M
    comp_cmq=[0.5*dens*u0*S*Cxu      -dens*u0*S*Cw0+0.5*dens*u0*S*Czu   0.5*dens*u0*S*MAC*data_1(1,4);      %u
              0.5*dens*u0*S*Cxa      0.5*dens*u0*S*Cza                 0.5*dens*u0*S*MAC*data_1(2,4);      %w
              0.25*dens*u0*S*MAC*Cxq 0.25*dens*u0*S*MAC*Czq            0.25*dens*u0*S*MAC^2*data_1(3,4);   %q
              0.25*dens*S*MAC*Cxap   0.25*dens*S*MAC*Czap              0.25*dens*S*MAC^2*data_1(4,4)];     %w·

    A_cmq=zeros(4);
    A_cmq(1,1)=comp_cmq(1,1)/m;
    A_cmq(1,2)=comp_cmq(2,1)/m;
    A_cmq(1,4)=-g;
    A_cmq(2,1)=comp_cmq(1,2)/(m-comp_cmq(4,2));
    A_cmq(2,2)=comp_cmq(2,2)/(m-comp_cmq(4,2));
    A_cmq(2,3)=(comp_cmq(3,2)+m*u0)/(m-comp_cmq(4,2));
    A_cmq(3,1)=(1/Iy)*(comp_cmq(1,3)+comp_cmq(4,3)*comp_cmq(1,2)/(m-comp_cmq(4,2)));
    A_cmq(3,2)=(1/Iy)*(comp_cmq(2,3)+comp_cmq(4,3)*comp_cmq(2,2)/(m-comp_cmq(4,2)));
    A_cmq(3,3)=(1/Iy)*(comp_cmq(3,3)+comp_cmq(4,3)*(comp_cmq(3,2)+m*u0)/(m-comp_cmq(4,2)));
    A_cmq(4,3)=1;
    
    cmq(i)=-10.25-2.05*(i-1);
    [autovectors_var_cmq,autovalors_var_cmq]=eig(A_cmq);
    autovalors_var_cmq=diag(autovalors_var_cmq);
    x=real(autovalors_var_cmq);
    y=imag(autovalors_var_cmq);
    scatter(x,y,25,'filled')
end
title('Evolució del valors propis variant C_m_q')
legend('-23','-22,5','-22','-21,5','-21')



%% QUESTION 3
% For the lowest value of Cmq plot the short period mode corresponding to
% an initial perturbed angle of attack alpha = 3º

tf=200; %[s]
Nstep=201;
tspan=linspace(0,tf,Nstep);
options=odeset("RelTol",1e-12,"AbsTol",1e-12);

xi=[0 0 0 0];

[t,x]=ode45(@section3,tspan,xi,options);
figure
hold on
grid on
plot(t(:),x(:,1))    
xlabel('t [s]')
ylabel('\Deltau [m/s]')

figure
hold on
grid on
plot(t(:),x(:,2))
xlabel('t [s]')
ylabel('w [m/s]')

figure
hold on
grid on
plot(t(:),x(:,3))
xlabel('t [s]')
ylabel('q [rad/s]')

figure
hold on
grid on
plot(t(:),x(:,4))
xlabel('t [s]')
ylabel('\Delta\Theta [rad]')



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
% perturbed mank angle phi=5�

tf=200; %[s]
Nstep=201;
tspan=linspace(0,tf,Nstep);
options=odeset("RelTol",1e-12,"AbsTol",1e-12);

xi=[0 0 0 0];

[t,x]=ode45(@section3_2,tspan,xi,options);
figure
hold on
grid on
plot(t(:),x(:,1)/u0)  %x(:,1)/u0 �s aproximadament beta
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


