% This matlab code solves the transmission coefficient (coeff) T and and reflection coefficient R for one dimensional (dim) potential shown below. 
% using a matching boundary conditions.  
%
%  Reference: R. L. Liboff, Introductory Quantum Mechanics, Addison-Wesley Publishing Company, New York (1080)  
%
% The atomic unit (au) is used in the calculation. 
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Nov 22, 2024 & University of North Dakota 
%
%            ________________                 
%           |                | V            
%   _______________________________________\ E (V > E)
%           |                |             /
% ______________________________________________ x
%          -a       0       -a  
% region I     region II       region III
%
%  region   I: psi_I(x)   = A*exp(i*k1*x) + B*exp(-i*k1*x)
%  region  II: psi_II(x)  = C*exp(k2*x)   + D*exp(-k2*x)
%  region III: psi_III(x) = F*exp(i*k1*x)
% 
%  boundary conditions: 
%  psi_I(-a) = psi_II(-a) & psi'_I(-a) = psi'_II(-a)
%  psi_II(a) = psi_III(a) * psi'_II(a) = psi'_III(a)
%
% Nov 22, 2024
function [] = one_dim_box_T_R_coeffs_1b
clc; format short
%
V =  0.500; % V > E case is presented. 
E0 = 0.100;
N_energy = 100;
delta_E = (V-E0)/(N_energy+1);
%
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('T_R_1d_test_2b.txt','w');
%
for ii = 1:N_energy
    %
    En = E0 + ii*delta_E;
    [T,R,T_plus_R,T_exact] = one_dim_box_scattering(En,V);
    %
    output = [ii, En, T, R, T_plus_R, T_exact];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %4.4f \t %4.12f \t %4.12f \t %4.12f \t %8.12f\n', output); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('T_R_1d_test_2b.txt', 'r');               % 
read_data  = textscan(read_data , '%f %f %f %f %f %f');
Energy = read_data{2};
T_coeff = read_data{3};
R_coeff = read_data{4};
TplusR_coeff = read_data{5};
T_coeff_exact = read_data{6};
%
figure(1)
hold on
plot(Energy, T_coeff, 'b', 'LineWidth', 1.5)
plot(Energy, R_coeff, 'g', 'LineWidth', 1.5)
plot(Energy, TplusR_coeff, 'r', 'LineWidth', 1.5)
%
plot(Energy(1:10:length(Energy)), T_coeff_exact(1:10:length(Energy)), 'ko','LineWidth', 1.5)
%
hold off
xlabel('\mbox{energy} (au)','Interpreter','latex') % ,'fontsize',16
ylabel('T \& R coeffs','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
%yscale log
box on
%%%
return
end

%%%
function [T,R,T_plus_R,T_exact] = one_dim_box_scattering(E,V)
%
%
ci = sqrt(-1.);
mass = 1.;
hbar2 = 1.;
%
aa = 1.;
%
k1 = 2*mass*E/hbar2;
k2 = 2*mass*(V - E)/hbar2;
%
d11 = exp(ci.*k1.*aa); d12 = -exp(-k2*aa); d13 = -exp(k2*aa); d14 = 0.;
d21 = -ci*k1*exp(ci.*k1*aa); d22 = -k2*exp(-k2*aa); d23 = k2*exp(k2*aa); d24 = 0.;
d31 = 0.; d32 = exp(k2*aa); d33 = exp(-k2*aa); d34 = -exp(ci*k1*aa);
d41 = 0.; d42 = k2*exp(k2*aa); d43 = -k2*exp(-k2*aa); d44 = -ci*k1*exp(ci*k1*aa);
%
D_mat = [d11, d12, d13, d14;
         d21, d22, d23, d24;
         d31, d32, d33, d34;
         d41, d42, d43, d44];
%
xx = [-exp(-ci*k1*aa);       % B/A
      -ci*k1*exp(-ci*k1*aa); % C/A
      0.;                    % D/A 
      0.];                   % F/A
%
coeff = D_mat\xx;
%
B_over_A = coeff(1);
F_over_A = coeff(4);
%
T = abs(F_over_A)^2 ;
R = abs(B_over_A)^2;
%
T_plus_R  = T + R;
%
T_exact = (1 + (1/4)*((k1^2+k2^2)/(k1*k2))^2 * sinh(2*k2*aa)^2)^(-1); % is consistent 
%%%
return
end

