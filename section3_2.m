function z=section3_2(t,x)

load('data_mv.mat');   % We load the input data from the main program

phi_0=5;

%Apartat 3
%Dades per al problema lateral
%               Cy     Cl     Cn
coef_lateral=[-0.90     0   0.160;   %beta 
                0   -0.340 -0.026;   %p
                0    0.130 -0.280; %r
                0     0     0; %delta a
                0     0     0];%delta r

%                  Y                              L                           N
comp_lateral=[0.5*dens*u0*S*coef_lateral(1,1)    0.5*dens*u0*S*b*coef_lateral(1,2)    0.5*dens*u0*S*b*coef_lateral(1,3);         %v
              0.25*dens*u0*S*b*coef_lateral(2,1) 0.25*dens*u0*S*b^2*coef_lateral(2,2) 0.25*dens*u0*S*b^2*coef_lateral(2,3);       %p
              0.25*dens*u0*S*b*coef_lateral(3,1) 0.25*dens*u0*S*b^2*coef_lateral(3,2) 0.25*dens*u0*S*b^2*coef_lateral(3,3)];      %r

Y_a=1/2*dens*u0^2*S*coef_lateral(4,1);
L_a=1/2*dens*u0^2*S*b*coef_lateral(4,2);
N_a=1/2*dens*u0^2*S*b*coef_lateral(4,3);

%Y_r=1/2*dens*u0^2*S*coef_lateral(5,1);
%L_r=1/2*dens*u0^2*S*b*coef_lateral(5,2);
N_r=1/2*dens*u0^2*S*b*coef_lateral(5,3);
          
AL=zeros(4,4);
%Fila 1
AL(1,1)=comp_lateral(1,1)/m;
AL(1,2)=comp_lateral(2,1)/m;
AL(1,3)=(comp_lateral(3,1)/m)-u0;
AL(1,4)=g*cosd(5); 

%Fila 2
AL(2,1)=comp_lateral(1,2)/Ixp+Ixzp*comp_lateral(1,3);
AL(2,2)=comp_lateral(2,2)/Ixp+Ixzp*comp_lateral(2,3);
AL(2,3)=comp_lateral(3,2)/Ixp+Ixzp*comp_lateral(3,3);

%Fila 3
AL(3,1)=Ixzp*comp_lateral(1,2)+comp_lateral(1,3)/Izp;
AL(3,2)=Ixzp*comp_lateral(2,2)+comp_lateral(2,3)/Izp;
AL(3,3)=Ixzp*comp_lateral(3,2)+comp_lateral(3,3)/Izp;

%Fila 4
AL(4,2)=1;
AL(4,3)=tand(5);

% if t<=40
%     sigma_a=8*pi/180*sind(9*t);
% else
%     sigma_a=0;
% end
% 
% B=zeros(4,1);
% B(1)=Y_a*sigma_a;
% B(2)=(L_a*sigma_a)/Ixp+Ixzp*(N_a);
% B(3)=Ixzp*L_a*sigma_a+(N_a*sigma_a)/Izp;
% 
% z=AL*x+B;

%El vector B realment tindria més valors però donat a que la deflexió de
%l'estabilitzador horitzontal és zero i que simplifiquem b22 i b32 com a 0
%ja que són valors molt petits sense cap efecte important a l'avió, només
%calculem b12
B=zeros(4,2);
B(1,2)=g*CD0/CL0;

if t<=3
    z=AL*x+B*[0; 0.5/3*t];
else
    z=AL*x+B*[0; 0.5];
end
end