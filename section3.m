function z=section3(t,x)
alpha = 3;
load('data_mv.mat');   % We load the input data from the main program

% Longitudinal Stability Derivatives
%          CT       CD      CL      Cm
data_1_3=[-0.055   0       0.13    0.013;  %û
        0        0.2     4.4     -1.0;   %alpha
        0        0       5.65    -10.25;  %^q (LOWEST VALUE OF Cmq)
        0        0       6.6     -4.0];  %alpha_prima

% Problema Longitudinal

Cxu = data_1_3(1,1)*cosd(alpha)-data_1_3(1,2)*cosd(alpha)+data_1_3(1,3)*sind(alpha);  % Cxu = CTu*cos(a) - CDu*cos(a) + CLu*sin(a)
Czu = -data_1_3(1,1)*sind(alpha)-data_1_3(1,3)*cosd(alpha)-data_1_3(1,3)*sind(alpha);  % Czu = CTu*sin(a) - CLu*cos(a) - CDu*sin(A)
Cxa = CL0-data_1_3(2,2);          % Cxa = CL0 - CDa
Cza = -data_1_3(2,3)-CD0;         % Cza = -CLa - CD0
Cxq = -data_1_3(3,2);             % Cxq = -CDq
Czq = -data_1_3(3,3);             % Czq = -CLq
Cxap = data_1_3(4,1)*cosd(alpha)-data_1_3(4,2)*cosd(alpha)+data_1_3(4,3)*sind(alpha); % Cxap = CTap*cos(a) - CDap*cos(a) + CLap*sin(a)
Czap = -data_1_3(4,1)*sind(alpha)-data_1_3(4,3)*cosd(alpha)-data_1_3(4,2)*sind(alpha);% Czap = -CTap*sin(a) - CLap*cos(a) - CDap*sin(a)


% Longitudinal Stability Derivatives
%                    X                                 Z                                M
ctrlderiv_3=[0.5*dens*u0*S*Cxu         -dens*u0*S*Cw0+0.5*dens*u0*S*Czu    0.5*dens*u0*S*MAC*data_1_3(1,4);      %u
             0.5*dens*u0*S*Cxa         0.5*dens*u0*S*Cza                   0.5*dens*u0*S*MAC*data_1_3(2,4);      %w
             0.25*dens*u0*S*MAC*Cxq    0.25*dens*u0*S*MAC*Czq              0.25*dens*u0*S*MAC^2*data_1_3(3,4);   %q
             0.25*dens*S*MAC*Cxap      0.25*dens*S*MAC*Czap                0.25*dens*S*MAC^2*data_1_3(4,4)];     %w·

A_3=zeros(4);

% Linearized Dynamic System
% Delta u Row
A_3(1,1)=ctrlderiv_3(1,1)/m;
A_3(1,2)=ctrlderiv_3(2,1)/m;
A_3(1,4)=-g;

% Delta w· Row
A_3(2,1)=ctrlderiv_3(1,2)/(m-ctrlderiv_3(4,2));
A_3(2,2)=ctrlderiv_3(2,2)/(m-ctrlderiv_3(4,2));
A_3(2,3)=(ctrlderiv_3(3,2)+m*u0)/(m-ctrlderiv_3(4,2));

% q· Row
A_3(3,1)=(1/Iy)*(ctrlderiv_3(1,3)+ctrlderiv_3(4,3)*ctrlderiv_3(1,2)/(m-ctrlderiv_3(4,2)));
A_3(3,2)=(1/Iy)*(ctrlderiv_3(2,3)+ctrlderiv_3(4,3)*ctrlderiv_3(2,2)/(m-ctrlderiv_3(4,2)));
A_3(3,3)=(1/Iy)*(ctrlderiv_3(3,3)+ctrlderiv_3(4,3)*(ctrlderiv_3(3,2)+m*u0)/(m-ctrlderiv_3(4,2)));

% Delta Theta· Row
A_3(4,3)=1;


%El vector B realment tindria més valors però donat a que la deflexió de
%l'estabilitzador horitzontal és zero i que simplifiquem b22 i b32 com a 0
%ja que són valors molt petits sense cap efecte important a l'avió, només
%calculem b12
B=zeros(4,2);
B(1,2)=g*CD0/CL0;

if t<=3
    z=A_3*x+B*[0; 0.5/3*t];
else
    z=A_3*x+B*[0; 0.5];
end

end

