%Treball mecanica de vol

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
CD0=0.235;
CL0=(2*m*g)/(dens*S*u0^2);
Cw0=(m*g)/(1/2*dens*u0^2*S);

%% Longitudinal Motion

%      CT       CD      CL      Cm
data_1=[ -0.055   0       0.13    0.013;  %û
       0       	0.2     4.4     -1.0;   %alpha
       0        0       5.65    -20.5;  %^q
       0        0       6.6     -4.0];  %dalpha

% Question 1
alpha=0;

Cxu = data_1(1,1)*cosd(alpha)-data_1(1,2)*cosd(alpha)+data_1(1,3)*sind(alpha);
Czu = -data_1(1,1)*sind(alpha)-data_1(1,3)*cosd(alpha)-data_1(1,2)*sind(alpha);
Cxa = (-data_1(2,2)+CL0);
Cza = -(data_1(2,2)+CD0);
Cxq = 0;
Czq = -data_1(3,3);
Cxap = data_1(4,1)*cosd(alpha)-data_1(4,2)*cosd(alpha)+data_1(4,3)*sind(alpha);
Czap = -data_1(4,1)*sind(alpha)-data_1(4,3)*cosd(alpha)-data_1(4,2)*sind(alpha);

%Es considera Theta0 = 0
%           X                       Z                           M
comp=[0.5*dens*u0*S*Cxu      -dens*u0*S*Cw0+0.5*dens*u0*S*Czu   0.5*dens*u0*S*MAC*data_1(1,4);      %u
      0.5*dens*u0*S*Cxa      0.5*dens*u0*S*Cza                 0.5*dens*u0*S*MAC*data_1(2,4);      %w
      0.25*dens*u0*S*MAC*Cxq 0.25*dens*u0*S*MAC*Czq            0.25*dens*u0*S*MAC^2*data_1(3,4);   %q
      0.25*dens*S*MAC*Cxap   0.25*dens*S*MAC*Czap              0.25*dens*S*MAC^2*data_1(4,4)];     %w·

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


%Aproximació curt període
Acurt=zeros(3);

Acurt(1,1)=A(2,2);
Acurt(1,2)=A(2,3);
Acurt(2,1)=A(3,2);
Acurt(2,2)=A(3,3);
Acurt(3,2)=A(4,3);

%Aproximació fugoide
Afug=zeros(3);

Afug(1,1)=A(1,1)-(A(3,1)/A(3,3))*A(1,3);
Afug(1,2)=A(1,2)-(A(3,2)/A(3,3))*A(1,3);
Afug(1,3)=A(1,4);
Afug(2,1)=A(2,1)-(A(3,1)/A(3,3))*A(2,3);
Afug(2,2)=A(2,2)-(A(3,2)/A(3,3))*A(2,3);
Afug(3,1)=-A(3,1)/A(3,3);
Afug(3,2)=-A(3,2)/A(3,3);

[autovectors,autovalors]=eig(A);
autovalors=diag(autovalors);
[autovectors_curt,autovalors_curt]=eig(Acurt);
autovalors_curt=diag(autovalors_curt);
[autovectors_fug,autovalors_fug]=eig(Afug);
autovalors_fug=diag(autovalors_fug);

