clear all
close all

% -------------------------------------------------------------------------
% Grid parameter
% -------------------------------------------------------------------------
I_L1        = 15;                       %Current Cable in A
I_L2        = 0;                        %Current Cable in A
I_L3        = 0;                        %Current Cable in A

phi_L1      = 0*pi/180;
phi_L2      = 120*pi/180;
phi_L3      = 240*pi/180;

I_L         = zeros(3,1);               %Conductor matrix for I
I_L(1,1)    = I_L1*exp(1i*phi_L1);
I_L(2,1)    = I_L2*exp(1i*phi_L2);
I_L(3,1)    = I_L3*exp(1i*phi_L3);

f           = 50;                       %Frequency in Hz

% -------------------------------------------------------------------------
% Cable parameter
% -------------------------------------------------------------------------
l           = 1000;                     %Length in m
%rho_S       = 0.0265;                   %spec. resistance OHM*m/mm^2 (alu) for the shield
rho_S       = 0.01678;                  %spec. resistance OHM*m/mm^2 (copper)
rho_C       = 0.0265;                   %spec. resistance OHM*m/mm^2 (alu) for den conductor

% conductor
% -------------------------------------------------------------------------
A_C         = 50 *1e-6;                 %conductor area in m²
d_C         = 0.008;
r_C         = d_C/2;                    %Radius conductor outside in m

% shield
% -------------------------------------------------------------------------
A_S         = 16*1e-6;                  %cable shield area in m²
d_S         = 28*1e-3;
r_S         = d_S/2;                    %Radius shield outside in m

% insulation
% -------------------------------------------------------------------------
%d_ISO       =0.0349;                   %diameter insulation
%A_ISO       =pi/4*d_ISO^2;
%d_ISO       =d_S-2*0.2*1e-3;
A_ISO       = pi/4*d_S^2-A_S;
d_ISO       = sqrt(A_ISO*4/pi);



% other calculation
% -------------------------------------------------------------------------
r_aS        = r_S;                      %Radius shield outside in m
r_iS        = d_ISO/2;                  %Radius shield inside in m
r_miso      = (r_aS+r_iS)/2;            %medium radius insulation
% -------------------------------------------------------------------------
%Installation parameter
% -------------------------------------------------------------------------
V_Art       = true;                     %Laying type bundles
%V_Art      =false;                     %Laying type bundles
V_Abb       = 0.060;                    %Laying distance
h           = 1;                        %Installation depth in m

% -------------------------------------------------------------------------
%Grounding parameter
% -------------------------------------------------------------------------
rho_e       = 100;                       %specific electrical resistivity of the soil in Ohm x m
Z_E         = 10 ;                       %Resistance of the grounding (also complex possible)

% -------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------
mue_0       = 1.26e-06;                 %Induction constant in Vs/ Am


% -------------------------------------------------------------------------
% global variable
% -------------------------------------------------------------------------
omega       = 2*pi*f;
alpha       = sqrt(1i*omega*mue_0/rho_e);  
delta_E     = 1/alpha;             	     %Earth current depth for infinitely long conductor in m (DUBANTON)
             

R_C         = rho_C/(A_C*1e6);           %Resistance conductor  in ohm/m  
R_RL        = omega*mue_0/8;             %Resistance earth return in ohm/m           
R_S         = rho_S/(1e6*pi*(r_aS^2-r_iS^2));%Resistance shield in ohm/m  





% -------------------------------------------------------------------------
% cable-cable-loop
% -------------------------------------------------------------------------
Z_CC=zeros(3,3);   

    for j = 1:3
        for m = 1:3
            if j == m
                Z_CC(j,m)=l*R_C+l*R_RL+l*1i*omega*mue_0/2/pi*log((2*h+delta_E)/r_C+1/4);
            else
                if  V_Art==true             %Bundle
                   if j == 1
                       Z_CC(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h-sin(60/pi*180)*V_Abb+2*delta_E+V_Abb)^2)/V_Abb);    
                   else
                       Z_CC(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb)/V_Abb);
                   end
                else                        %Flat
                     Z_CC(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb*abs(j-m))/V_Abb*abs(j-m));
               end
    
            end
        end
    end
    
% -------------------------------------------------------------------------
% cable-shield-loop
% -------------------------------------------------------------------------
Z_CS=zeros(3,3);
   
    for j = 1:3
        for m = 1:3
            if j == m
                Z_CS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log((2*h+delta_E)/r_miso);
            else
               if  V_Art==true             %Bundle
                   if j == 1
                       Z_CS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h-sin(60/pi*180)*V_Abb+2*delta_E+V_Abb)^2)/V_Abb);
                   else
                       Z_CS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb)/V_Abb);
                   end
               else                        %Flat
                  Z_CS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb*abs(j-m))/V_Abb*abs(j-m));
               end
            end
        end
    end
    
% -------------------------------------------------------------------------
% shield-cable-loop
% -------------------------------------------------------------------------
Z_SC=transpose(Z_CS);
% -------------------------------------------------------------------------
% shield-shield-loop
% -------------------------------------------------------------------------
   

Z_SS=zeros(3,3);
    
    for j = 1:3
        for m = 1:3
            if j == m
                Z_SS(j,m)=l*R_S+l*R_RL+l*1i*omega*mue_0/2/pi*log((2*h+delta_E)/r_aS);
            else
               if  V_Art==true             %bundle
                   if j == 1
                       Z_SS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h-sin(60/pi*180)*V_Abb+2*delta_E+V_Abb)^2)/V_Abb);
                   else
                       Z_SS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb)/V_Abb);
                   end
               else                        %flat
               Z_SS(j,m)=l*R_RL+l*1i*omega*mue_0/2/pi*log(sqrt((2*h+2*delta_E)^2+V_Abb*abs(j-m))/V_Abb*abs(j-m));
               end
            end
        end
    end

% -------------------------------------------------------------------------
% impedance matrix
% -------------------------------------------------------------------------
Z_all   =[Z_CC Z_CS;Z_SC Z_SS];     

% -------------------------------------------------------------------------
% determiantion of the Currents
% -------------------------------------------------------------------------
I_Sm    = -1*(inv(Z_SS+Z_E))*(Z_SC+Z_E)*I_L;
I_Em    =-1*(I_L+I_Sm);
I_E     =sum(I_Em);

% -------------------------------------------------------------------------
% ground voltage
% -------------------------------------------------------------------------

U_E     =I_E*Z_E;
U_E_abs =abs(U_E);
disp(U_E_abs);
% -------------------------------------------------------------------------
% reduction factor
% -------------------------------------------------------------------------
r       =abs(I_E/(sum(I_L)));
r_abs   =abs(r);
disp(r_abs);
    



