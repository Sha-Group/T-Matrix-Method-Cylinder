% This Matlab code calculates the scattering field of multiple multilayer 2-D plasma cylinders using T-matrix method. This method can also be used
% for approximately estimating the scattering field of long and thin 3-D cylinders.
% The T-matrix of a single cylinder can be calculated by generalized reflection coefficient method, then a T-matrix equation is formed for the multiple scattering problem.
% Notice: Large number of Mie series is needed for larger-scale problem. T-matrix method is suitable for calculating the scattering of many compact (~10 wavelength) scatters.
% Reference: S. S. A. Yuan, Z. H. Lin, L. -B. Lv, S. -J. Hao and W. E. I. Sha, "Investigating the Scattering Characteristics of Artificial Field-Aligned Irregularities Based on T-Matrix Algorithm," in IEEE Journal on Multiscale and Multiphysics Computational Techniques, vol. 8, pp. 147-157, 2023.
% Shuai S. A. Yuan
% Email: shuaiyuan1997@zju.edu.cn


function T_matrix_for_Multiple_Cylinders
clc;clear
%% input
Wavelength=55; % Wavelength (m)
Polarization=1; % Polarizations (1-TM wave, 2-TE wave, 3-left circularly polarized wave, 4-right circularly polarized wave)
Radius=55;% Radius of a single cylinder, the cylinders are along the z axis
Layers=10;% Number of layers of a cylinder, the most outside layer is air layer
electron_density=6.5e11; 
central_density=5*electron_density; % central electron density
Angle=90; % Angle between the incident plane wave and z axis

% Number and positons of the cylinders
%X_position=[0,0,0,3,-3,-3,3,3,-3]*Radius; 
%Y_position=[0,3,-3,0,0,-3,3,-3,3]*Radius;
%X_position=[0]*Radius; % Single cylinder x  
%Y_position=[0]*Radius; % Single cylinder y
X_position=[-3,3]*Radius; 
Y_position=[-3,3]*Radius; 

%% Parameters
ind_x=Wavelength;
Pn=length(X_position);
%incident angle
incident_angle=Angle/180*pi; 
phi_angle=0/180*pi;
% basic parameter
epsc=1/(4*pi*9*10.^9);      %  permittivity in free space
murc=4*pi*10.^(-7);         %  permeability in free space
electron_charge=1.60217*10.^(-19);% electron charge (coulombs)
electron_mass=9.10938*10.^(-31); %electron mass (kg)
collision_frequency=1e3; % loss in drude model (rad/s)
num_layer=Layers; %number of layers
FAI_radius=Radius; %radius of FAI (m)
FAI_r=linspace(eps,eps+FAI_radius,num_layer);
FAI_w=FAI_radius*sqrt(-1/(log(electron_density/central_density))); % gaussian distribution coefficients
layer_density=central_density*exp(-FAI_r.^2/FAI_w^2);
%  relative permittivity
Air=ones(1,length(ind_x))*1;                 %  air
%Si_epr=3.4^2*ones(1,length(x))*1;       %  silicon

%% Wavenumber of incident wave
theta=incident_angle;
phi=phi_angle;
wavelength=ind_x;%m
k_0=2*pi/wavelength;
k_z=k_0*cos(theta);
k_rho=k_0*sin(theta);
k_x=k_rho*cos(phi);
k_y=k_rho*sin(phi);

omega=2*pi*3*10^8./wavelength;%
MAX=ceil(2*pi*FAI_radius/wavelength+4.05*(2*pi*FAI_radius/wavelength)^(1/3)+2); % Number of Mie series
L=num_layer;                 

% coefficient a
Z0=sqrt(murc/epsc);
% outmost layer coefficient for cylinder
% a at L layer, a 1x2 matrix amplitude for TE and TM
a=zeros(2,L,MAX);
for n=-MAX:1:MAX   % for all harmonics % - infinite to + infinite, different from spherical case
    %a(:,L,n+MAX+1)=[1j^(-n);1j^(-n)/Z0];  %
    a(:,L,n+MAX+1)=[0;(1j^(-n)/Z0)]; % TM (1j^(n)) and TE (1j^(n)/Z0)
end%
%% 广义反射系数
for fre=1:length(wavelength)
    fre
    %eps_array=[Si_epr(fre),Air(fre)]; %  relative permittivity (inner to outer)
    %eps_array=[eps_save(fre),Air(fre),eps_save(fre),Air(fre)]; %  relative permittivity (inner to outer)
    %eps_array=[eps_save(fre),Air(fre),eps_save(fre),Air(fre)]; %  relative permittivity (inner to outer)
    mu_array=ones(1,num_layer);                   %  relative permeability
    rho_array=FAI_r(2:end);                  %  radius (inner to outer)
    % drude model for ionosphere plasma
    plasma_fre=sqrt(electron_charge.^2*layer_density/epsc/electron_mass);
    eps_plasma=1-(electron_charge.^2*layer_density/epsc/electron_mass/(omega(fre).^2+1i*omega(fre)*collision_frequency));%not relative
    %eps_plasma=ones(1,num_layer);
    %eps_plasma(1:L-1)=3;
    eps_plasma(L)=1; %use air for most out layer

    k_array=sqrt(mu_array.*eps_plasma)*2*pi/(wavelength(fre));
    kz_array=k_array*cos(incident_angle);
    krho_array=k_array*sin(incident_angle);
    eps_array=eps_plasma*epsc;
    %eps_array=eps_array*epsc;
    mu_array=mu_array*murc;
    rho_array=rho_array;%m
    
    %  inverse the order
    inv_k_arr=k_array(end:-1:1);
    inv_kz_arr=kz_array(end:-1:1);
    inv_krho_arr=krho_array(end:-1:1);
    inv_eps_arr=eps_array(end:-1:1);
    inv_mu_arr=mu_array(end:-1:1);
    inv_rho_arr=rho_array(end:-1:1);
    %  Outgoing wave R and T coef matrix (inner to outer)
    for m=1:L-1;
        for n=-MAX:1:MAX
            R_Out(:,:,m,n+MAX+1)=...  %inv(...) = D matrix
                inv(J_matrix(n,krho_array(m),kz_array(m),rho_array(m),eps_array(m),omega(fre),mu_array(m))...
                *cb2(n,krho_array(m+1),rho_array(m))-H_matrix(n,krho_array(m+1),kz_array(m+1),rho_array(m),...
                eps_array(m+1),omega(fre),mu_array(m+1))*cb1(n,krho_array(m),rho_array(m)))*( cb2(n,...
                krho_array(m),rho_array(m))*H_matrix(n,krho_array(m+1),kz_array(m+1),rho_array(m),...
                eps_array(m+1),omega(fre),mu_array(m+1))- cb2(n,krho_array(m+1),rho_array(m))*H_matrix(n,...
                krho_array(m),kz_array(m),rho_array(m),eps_array(m),omega(fre),mu_array(m)) );
            T_Out(:,:,m,n+MAX+1)=...
                inv(J_matrix(n,krho_array(m),kz_array(m),rho_array(m),eps_array(m),omega(fre),mu_array(m))...
                *cb2(n,krho_array(m+1),rho_array(m))-H_matrix(n,krho_array(m+1),kz_array(m+1),rho_array(m),...
                eps_array(m+1),omega(fre),mu_array(m+1))*cb1(n,krho_array(m),rho_array(m)))*( cb2(n,...
                krho_array(m),rho_array(m))*J_matrix(n,krho_array(m),kz_array(m),rho_array(m),...
                eps_array(m),omega(fre),mu_array(m))- cb1(n,krho_array(m),rho_array(m))*H_matrix(n,...
                krho_array(m),kz_array(m),rho_array(m),eps_array(m),omega(fre),mu_array(m)) );
            T_Out_test(:,:,m,n+MAX+1)=...
                inv(J_matrix(n,krho_array(m),kz_array(m),rho_array(m),eps_array(m),omega(fre),mu_array(m))...
                *cb2(n,krho_array(m+1),rho_array(m))-H_matrix(n,krho_array(m+1),kz_array(m+1),rho_array(m),...
                eps_array(m+1),omega(fre),mu_array(m+1))*cb1(n,krho_array(m),rho_array(m)))*[eps_array(m),0;...
                0,-mu_array(m)]*2*omega(fre)/pi/krho_array(m)^2/rho_array(m);
        end
    end
    
    %  Standing wave R and T coef matrix (outer to inner)
    %  R_43  R_32  R_21
    for m=1:L-1;
        for n=-MAX:1:MAX
            R_In_inv(:,:,m,n+MAX+1)=...  %inverse the eps mu k
                inv(J_matrix(n,inv_krho_arr(m+1),inv_kz_arr(m+1),inv_rho_arr(m),inv_eps_arr(m+1),omega(fre),inv_mu_arr(m+1))...
                *cb2(n,inv_krho_arr(m),inv_rho_arr(m))-H_matrix(n,inv_krho_arr(m),inv_kz_arr(m),inv_rho_arr(m),...
                inv_eps_arr(m),omega(fre),inv_mu_arr(m))*cb1(n,inv_krho_arr(m+1),inv_rho_arr(m)))*( cb1(n,...
                inv_krho_arr(m+1),inv_rho_arr(m))*J_matrix(n,inv_krho_arr(m),inv_kz_arr(m),inv_rho_arr(m),...
                inv_eps_arr(m),omega(fre),inv_mu_arr(m))- cb1(n,inv_krho_arr(m),inv_rho_arr(m))*J_matrix(n,...
                inv_krho_arr(m+1),inv_kz_arr(m+1),inv_rho_arr(m),inv_eps_arr(m+1),omega(fre),inv_mu_arr(m+1)) );
            T_In_inv(:,:,m,n+MAX+1)=...
                inv(J_matrix(n,inv_krho_arr(m+1),inv_kz_arr(m+1),inv_rho_arr(m),inv_eps_arr(m+1),omega(fre),inv_mu_arr(m+1))...
                *cb2(n,inv_krho_arr(m),inv_rho_arr(m))-H_matrix(n,inv_krho_arr(m),inv_kz_arr(m),inv_rho_arr(m),...
                inv_eps_arr(m),omega(fre),inv_mu_arr(m))*cb1(n,inv_krho_arr(m+1),inv_rho_arr(m)))*[inv_eps_arr(m),0;...
                0,-inv_mu_arr(m)]*2*omega(fre)/pi/inv_krho_arr(m)^2/inv_rho_arr(m);
            %notice the region and inverse, layer order exchange
        end
    end
    %  Standing wave R and T coef matrix（inner to outer）
    %  R_21  R_32  R_43
    R_In=R_In_inv(:,:,end:-1:1,:);
    T_In=T_In_inv(:,:,end:-1:1,:);
    R_In_test=R_In(:,:,end,6);
    
    %  Generalized R and T coef matrix Ref
    Gen_R(:,:,1,:)=R_In(:,:,1,:);%first layer
    for m=2:L-1
        for n=-MAX:1:MAX
            Gen_R(:,:,m,n+MAX+1)=R_In(:,:,m,n+MAX+1)+(T_Out(:,:,m,n+MAX+1)*Gen_R(:,:,m-1,n+MAX+1)...
                *inv(eye(2)-R_Out(:,:,m,n+MAX+1)*Gen_R(:,:,m-1,n+MAX+1))*T_In(:,:,m,n+MAX+1));
        end
    end
    %test_Tout(:,:,:)=T_Out(:,:,:,1);
    test_Gen_R(:,:)=Gen_R(:,:,end,10);
    % get coefficients ampititude
    for m=L-1:-1:2
        for n=-MAX:1:MAX
            a(:,m,n+MAX+1)=inv(eye(2)-R_Out(:,:,m,n+MAX+1)*Gen_R(:,:,m-1,n+MAX+1))*T_In(:,:,m,n+MAX+1)*a(:,m+1,n+MAX+1);
        end
    end
    for n=-MAX:1:MAX %most inner interface
        a(:,1,n+MAX+1)=T_In(:,:,1,n+MAX+1)*a(:,2,n+MAX+1);
    end
    test_a(:,:,:)=a(:,1,1);
    %% b vector
    b=zeros((2*MAX+1)*2,1); % Including TM and TE
    b_vector=zeros((2*MAX+1)*2*Pn,1);% Final b vector
    for n=-MAX:1:MAX   % for all harmonics
        pol_index=[1,0,-1j,1j;0,1,1,1];
        b(n+MAX+1)=1j^(-n)*pol_index(1,Polarization);%1j^(-n)
        b(n+3*MAX+2)=1j^(-n)/Z0*pol_index(2,Polarization);%1j^(-n)/Z0
        %     b_test(n+MAX+1)=a(1,L,n+MAX+1);
        %     a_test(n+MAX+1)=Gen_R(1,1,L-1,n+MAX+1)*a(1,L,n+MAX+1);
    end
    for p=1:Pn 
        b_vector((p-1)*(4*MAX+2)+1:p*(4*MAX+2))=b*exp(-1j*(k_x*X_position(p)+k_y*Y_position(p)));
    end
    %% T matrix
    T_i=zeros((2*MAX+1)*2,(2*MAX+1)*2);
    T_matrix=zeros((2*MAX+1)*2*Pn,(2*MAX+1)*2*Pn);
    for n=-MAX:1:MAX   % for all harmonics
        T_MM(n+MAX+1)=Gen_R(1,1,L-1,n+MAX+1);
        T_MN(n+MAX+1)=Gen_R(1,2,L-1,n+MAX+1);
        T_NM(n+MAX+1)=Gen_R(2,1,L-1,n+MAX+1);
        T_NN(n+MAX+1)=Gen_R(2,2,L-1,n+MAX+1);
    end
    T_i=[diag(T_MM),diag(T_MN);diag(T_NM),diag(T_NN)];% T matrix of each scatterer
    for pp=1:Pn % Phase difference between each scatterer
        for qq=1:Pn
            if pp==qq
                T_matrix((pp-1)*(4*MAX+2)+1:pp*(4*MAX+2),(qq-1)*(4*MAX+2)+1:qq*(4*MAX+2))=T_i;
            end
        end
    end
    %% Translation matrix
    Alpha_matrix=zeros((2*MAX+1)*2*Pn,(2*MAX+1)*2*Pn);
    for pp=1:Pn 
        for qq=1:Pn
            if pp~=qq
                Alpha_matrix((pp-1)*(4*MAX+2)+1:pp*(4*MAX+2),(qq-1)*(4*MAX+2)+1:qq*(4*MAX+2))=...
                    translation_matrix(X_position(pp),Y_position(pp),X_position(qq),Y_position(qq),MAX,k_rho);
                T_test=translation_matrix(X_position(pp),Y_position(pp),X_position(qq),Y_position(qq),MAX,k_rho);
            end
        end
    end 
    %% Matrix inversion
    matrix_test=T_matrix*Alpha_matrix;
%    rank_test=rank(matrix_test);
    a_vector=(eye((2*MAX+1)*2*Pn)-T_matrix*Alpha_matrix)\(T_matrix*b_vector);
    a_vector_t=a_vector(1:(4*MAX+2));

    a_vector_sum=zeros(4*MAX+2,1);
    for tt=1:Pn
        a_vector_sum=a_vector_sum+translation_matrix2(0,0,X_position(tt),Y_position(tt),MAX,k_rho)*...
            a_vector((tt-1)*(4*MAX+2)+1:tt*(4*MAX+2));
    end%

    %% RCS calculation
    % get near field (sample)
    rrho=linspace(eps,2*rho_array(end)+eps,300);
    %zz=linspace(-2*rho_array(end)+eps,2*rho_array(end)+eps,round(4*rho_array(end)/0.2*10^9));
    zz=[0];
    phi=linspace(eps,2*pi+eps,360);
    %phi=[0,pi];
    
    % e-field components
    Field_z=zeros(2,length(rrho),length(zz),length(phi));
    Farfield_z=zeros(2,length(phi));
    Field_phi=zeros(2,length(rrho),length(zz),length(phi));
    Field_rho=zeros(2,length(rrho),length(zz),length(phi));
    % calculate Ez and Hz of each layer
    
    % most outer part  % air
    s=L;
    R_test=Gen_R(:,:,s-1,1);
    % farfield E_z for RCS calculation
    s=L;
    for q=-MAX:1:MAX  % order
        for t=1:length(phi)% phi %
            %             Farfield_z(:,t)=Farfield_z(:,t)+2/(sqrt(krho_array(s)))*exp(1j*q*phi(t)-1j*q/2*pi)...
            %                 *Gen_R(:,:,s-1,q+MAX+1)*a(:,s,q+MAX+1);
            Farfield_z(:,t)=Farfield_z(:,t)+exp(1j*q*phi(t)+1j*q/2*pi)...
                *Gen_R(:,:,s-1,q+MAX+1)*a(:,s,q+MAX+1);
            exp_test(q+MAX+1)=exp(1j*q*phi(t)-1j*q/2*pi);
            % the tangential field is very small when rho -> infinite
        end
    end
end
%% RCS， T matrix method
RCS_TM=zeros(1,(2*MAX+1));
RCS_TE=zeros(1,(2*MAX+1));
for t=1:length(phi)  % phi %
    for q=-MAX:1:MAX   %
        if abs(q)<=MAX
        RCS_TM(q+MAX+1)=exp(1j*q*phi(t)+1j*q/2*pi);
        RCS_TE(q+MAX+1)=exp(1j*q*phi(t)+1j*q/2*pi)*Z0;
        end
    end
    RCS(t)=2/k_rho/sin(theta)*(abs(RCS_TM*a_vector_sum(1:2*MAX+1)).^2+...
        abs(RCS_TE*a_vector_sum(2*MAX+2:4*MAX+2)).^2);
    
end
length_cy=20*ind_x; % Length of cylinder for approximately calculating the RCS of finite cylinders 
RCS=10.*log10(RCS*2*sin(theta)*(length_cy)^2/ind_x);%
plot(phi,RCS)
xlabel('Phi')
ylabel('RCS (dB)')
hold on

%% Cylindrical harmonic functions
% notice i and j
% Cylindrical Bessel 1 (J_n)
function output=cb1(degree,k_rho,rho)
output=besselj(degree,k_rho*rho);

% Derivative of Cylindrical Bessel 1 (J_n')
% function output=cb_der1(degree,k_rho,rho)
% output=(degree*cb1(degree-1,k_rho,rho)-(degree+1)*cb1(degree+1,k_rho,rho))/(2*degree+1);
function output=cb_der1(degree,k_rho,rho)
output=(cb1(degree-1,k_rho,rho)-cb1(degree+1,k_rho,rho))/2;
%different from the book?j_{n-1}-n/(k*rho)*j_{n}%

% Cylindrical Bessel 2 (H^1_n)
function output=cb2(degree,k_rho,rho)
output=besselh(degree,1,k_rho*rho);
% Derivative of Cylindrical Bessel 2 (H^1_n')
% function output=cb_der2(degree,k_rho,rho)
% output=(degree*cb2(degree-1,k_rho,rho)-(degree+1)*cb2(degree+1,k_rho,rho))/(2*degree+1);
function output=cb_der2(degree,k_rho,rho)
output=(cb2(degree-1,k_rho,rho)-cb2(degree+1,k_rho,rho))/2;

% J matrix
function output=J_matrix(degree,k_rho,k_z,rho,eps,omega,mu)

output=1/(k_rho^2*rho)*[1j*omega*eps*k_rho*rho*cb_der1(degree,k_rho,rho),-1*degree*k_z*cb1(degree,k_rho,rho);
     -1*degree*k_z*cb1(degree,k_rho,rho), -1j*omega*k_rho*rho*mu*cb_der1(degree,k_rho,rho) ];%

% H matrix
function output=H_matrix(degree,k_rho,k_z,rho,eps,omega,mu)

 output=1/(k_rho^2*rho)*[1j*omega*eps*k_rho*rho*cb_der2(degree,k_rho,rho),-1*degree*k_z*cb2(degree,k_rho,rho);
     -1*degree*k_z*cb2(degree,k_rho,rho), -1j*omega*mu*k_rho*rho*cb_der2(degree,k_rho,rho) ];
% modified J and H matrix for the calculation of rho component of E field

% modified J matrix
function output=J_matrix2(degree,k_rho,k_z,rho,eps,omega,mu)

output=1/(k_rho^2*rho)*[omega*eps*degree*cb1(degree,k_rho,rho),1j*k_z*k_rho*rho*cb_der1(degree,k_rho,rho);
    1j*k_z*k_rho*rho*cb_der1(degree,k_rho,rho), -1*omega*mu*degree*cb1(degree,k_rho,rho) ];

% modified H matrix
function output=H_matrix2(degree,k_rho,k_z,rho,eps,omega,mu)

output=1/(k_rho^2*rho)*[omega*eps*degree*cb2(degree,k_rho,rho),1j*k_z*k_rho*rho*cb_der2(degree,k_rho,rho);
    1j*k_z*k_rho*rho*cb_der2(degree,k_rho,rho), -1*omega*mu*degree*cb2(degree,k_rho,rho) ];

