% phase matching stuff for RKTP, assuming first order QPM.
function[QPM,n_p,n_s,n_i,T] = phase_match()

clear
clc
close all
clear global

T = 300;  % Temperature, kelvin.
L = 12000;
sig_range_pos = 4;   %um
sig_range_neg = 1.2; %um

lamp_mu = 1.0642;
lams_mu = linspace(sig_range_neg*lamp_mu,sig_range_pos*lamp_mu,1000);
lami_mu = zeros(1,length(lams_mu));

n_i = zeros(1,length(lams_mu));
n_p = n_i;
n_s = n_i;

% Sellmeier stuff (Andrius thesis KTP/RKTP).

if lamp_mu >= 1
    
    A = 2.12725; B = 1.18431; C = 0.0514852; D = 0.6603; E = 100.00507; F = 0.00968956;
    
    a1_0 =  9.9587*10^-6; a1_1 = 9.9228*10^-6; a1_2 = -8.9603*10^-6; a1_3 = 4.1010*10^-6;  % n1 param.
    a2_0 = -1.1882*10^-8; a2_1 = 10.459*10^-8; a2_2 = -9.8136*10^-8; a2_3 = 3.1481*10^-8;  % n2 param.   
    
    % del_n = n1*(T - 25C) + n2*(T-25)^2;   % Eq 4.6.3
    
    n_p  = sqrt(A + (B/(1-C*(lamp_mu^-2))) + (D/(1-E*(lamp_mu^-2))) - F*lamp_mu^2);
    
    n1 = a1_0 + (a1_1/lamp_mu) + (a1_2/(lamp_mu^2)) + (a1_3/(lamp_mu^3));
    n2 = a2_0 + (a2_1/lamp_mu) + (a2_2/(lamp_mu^2)) + (a2_3/(lamp_mu^3)); 
    
    deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.
    
    n_p = n_p + deln;
    
else
    
    A = 2.25411; B = 1.06543; C = 0.05486; D = 0; E = 0; F = 0.02140;
    
    a =  1.2415*10^-5;  b = -4.4414*10^-5; c = 5.9129*10^-5; d = -1.2101*10^-5;
    
    % deln = (a/lam3 + b/lam2 + c/lam1 + d)*dT.  % Eq 4.6.2
    
    n_p  = sqrt(A + (B/(1-C*(lamp_mu^-2))) + (D/(1-E*(lamp_mu^-2))) - F*lamp_mu^2);
    deln = ((a/(lamp_mu^3)) + (b/(lamp_mu^2)) + (c/lamp_mu) + d)*T;
    n_p  = n_p + deln;
    
end

for i = 1:length(lams_mu)
    
    lami_mu(i) = ((1./lamp_mu) - (1./lams_mu(i))).^-1;

    if lams_mu >=1;
        
        A = 2.12725; B = 1.18431; C = 0.0514852; D = 0.6603; E = 100.00507; F = 0.00968956;

        a1_0 =  9.9587*10^-6; a1_1 = 9.9228*10^-6; a1_2 = -8.9603*10^-6; a1_3 = 4.1010*10^-6;  % n1 param.
        a2_0 = -1.1882*10^-8; a2_1 = 10.459*10^-8; a2_2 = -9.8136*10^-8; a2_3 = 3.1481*10^-8;  % n2 param.

        n_s(i) = sqrt(A + (B/(1-C*(lams_mu(i)^-2))) + (D/(1-E*(lams_mu(i)^-2))) - F*lams_mu(i)^2);
        n1 = a1_0 + (a1_1/lams_mu(i)) + (a1_2/(lams_mu(i)^2)) + (a1_3/(lams_mu(i)^3));
        n2 = a2_0 + (a2_1/lams_mu(i)) + (a2_2/(lams_mu(i)^2)) + (a2_3/(lams_mu(i)^3));    
        deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.
        n_s(i) = n_s(i) + deln;

        n_i(i) = sqrt(A + (B/(1-C*(lami_mu(i)^-2))) + (D/(1-E*(lami_mu(i)^-2))) - F*lami_mu(i)^2);
        n1 = a1_0 + (a1_1/lami_mu(i)) + (a1_2/(lami_mu(i)^2)) + (a1_3/(lami_mu(i)^3));
        n2 = a2_0 + (a2_1/lami_mu(i)) + (a2_2/(lami_mu(i)^2)) + (a2_3/(lami_mu(i)^3)); 

        deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.   
        n_i(i) = n_i(i) + deln; 

    
    else
    
        A = 2.25411; B = 1.06543; C = 0.05486; D = 0; E = 0; F = 0.02140;

        a =  1.2415*10^-5;  b = -4.4414*10^-5; c = 5.9129*10^-5; d = -1.2101*10^-5;

        % deln = (a/lam3 + b/lam2 + c/lam1 + d)*dT.  % Eq 4.6.2

        n_s(i)  = sqrt(A + (B/(1-C*(lams_mu(i)^-2))) + (D/(1-E*(lams_mu(i)^-2))) - F*lams_mu(i)^2);
        deln = ((a/(lams_mu(i)^3)) + (b/(lams_mu(i)^2)) + (c/lams_mu(i)) + d)*T;
        n_s(i)  = n_s(i) + deln;
        
            if lami_mu >=1
                
                n_i(i) = sqrt(A + (B/(1-C*(lami_mu(i)^-2))) + (D/(1-E*(lami_mu(i)^-2))) - F*lami_mu(i)^2);
                n1 = a1_0 + (a1_1/lami_mu(i)) + (a1_2/(lami_mu(i)^2)) + (a1_3/(lami_mu(i)^3));
                n2 = a2_0 + (a2_1/lami_mu(i)) + (a2_2/(lami_mu(i)^2)) + (a2_3/(lami_mu(i)^3)); 

                deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.   
                n_i(i) = n_i(i) + deln; 
                
            else
                
                n_i(i)  = sqrt(A + (B/(1-C*(lami_mu(i)^-2))) + (D/(1-E*(lami_mu(i)^-2))) - F*lami_mu(i)^2);
                deln = ((a/(lami_mu(i)^3)) + (b/(lami_mu(i)^2)) + (c/lami_mu(i)) + d)*T;
                n_i(i)  = n_i(i) + deln;

            end

    
    end
    
    QPM(i) = ((n_p./lamp_mu)-(n_s(i)./lams_mu(i))-(n_i(i)./lami_mu(i))).^-1;
    deltaK(i) = (n_p./lamp_mu)-(n_s(i)./lams_mu(i))-(n_i(i)./lami_mu(i))-(1./(QPM(i)));
%     Int(i) = sinc(deltaK(i).*L/2);
    
end


% for i = 1:length(lams_mu)
%     for j = 1:length(QPM)
%         
%         deltaKK(i,j) = (n_p./lamp_mu)-(n_s(i)./lams_mu(i))-(n_i(i)./lami_mu(i))-(1./(QPM(j)));
% 
%         Int(i,j) = sinc(deltaKK(i,j).*L/2).^2;
%     
%     end
% end

figure1 = figure('Position', [100, 100, 1500, 800]);
figure1;

subplot(2,3,1)
% [x,y] = meshgrid(QPM,lams_mu);
% mesh(x,y,Int); colorbar;
% xlabel('QPM period');
% ylabel('Wavenlength ({\mu}m)')
% view([0 90])

subplot(2,3,2)
plot(QPM,lami_mu,'r','Linewidth',2)
hold on
grid on
plot(QPM,lams_mu,'b','Linewidth',2)
xlabel('QPM period ({\mu}m)')
ylabel('Wavelength ({\mu}m)')
legend('Idler','Signal')
title('Phase matching curve for pump at 1.0642 {\mu}m')


%% need to do Sinc function plot, QPM is kept fixed, deltaK is calculated
% for fixed QPM and varing pump wavelength.

% QPM = 38.8603; % um, for perfect degeneracy signal/idler.

QPM = 38.8603;

lamp_mu = linspace(1.0642-0.05,1.0642+0.05,3000);
lams_mu = linspace(2.125-0.05,2.130+0.05,3000);
lami_mu = zeros(1,length(lamp_mu));

% just assume all wavelengths > 1 um for now.

for i = 1:length(lamp_mu)
    
    lami_mu(i) = ((1./lamp_mu(i)) - (1./lams_mu(i))).^-1;

    A = 2.12725; B = 1.18431; C = 0.0514852; D = 0.6603; E = 100.00507; F = 0.00968956;

    a1_0 =  9.9587*10^-6; a1_1 = 9.9228*10^-6; a1_2 = -8.9603*10^-6; a1_3 = 4.1010*10^-6;  % n1 param.
    a2_0 = -1.1882*10^-8; a2_1 = 10.459*10^-8; a2_2 = -9.8136*10^-8; a2_3 = 3.1481*10^-8;  % n2 param.   

    % del_n = n1*(T - 25C) + n2*(T-25)^2;   % Eq 4.6.3

    n_p(i)  = sqrt(A + (B/(1-C*(lamp_mu(i)^-2))) + (D/(1-E*(lamp_mu(i)^-2))) - F*lamp_mu(i)^2);

    n1 = a1_0 + (a1_1/lamp_mu(i)) + (a1_2/(lamp_mu(i)^2)) + (a1_3/(lamp_mu(i)^3));
    n2 = a2_0 + (a2_1/lamp_mu(i)) + (a2_2/(lamp_mu(i)^2)) + (a2_3/(lamp_mu(i)^3)); 

    deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.

    n_p(i) = n_p(i) + deln;

    n_s(i) = sqrt(A + (B/(1-C*(lams_mu(i)^-2))) + (D/(1-E*(lams_mu(i)^-2))) - F*lams_mu(i)^2);
    n1 = a1_0 + (a1_1/lams_mu(i)) + (a1_2/(lams_mu(i)^2)) + (a1_3/(lams_mu(i)^3));
    n2 = a2_0 + (a2_1/lams_mu(i)) + (a2_2/(lams_mu(i)^2)) + (a2_3/(lams_mu(i)^3));    
    deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.
    n_s(i) = n_s(i) + deln;
    
    n_i(i) = sqrt(A + (B/(1-C*(lami_mu(i)^-2))) + (D/(1-E*(lami_mu(i)^-2))) - F*lami_mu(i)^2);
    n1 = a1_0 + (a1_1/lami_mu(i)) + (a1_2/(lami_mu(i)^2)) + (a1_3/(lami_mu(i)^3));
    n2 = a2_0 + (a2_1/lami_mu(i)) + (a2_2/(lami_mu(i)^2)) + (a2_3/(lami_mu(i)^3)); 

    deln = n1*(T-298.15) + n2*(T-298.15)^2;  % converted 25 degree celcius to kelvin here - 298.15.   
    n_i(i) = n_i(i) + deln; 
    
    deltaK(i) = (n_p(i)./lamp_mu(i))-(n_s(i)./lams_mu(i))-(n_i(i)./lami_mu(i))-(1./(QPM));
    
end

L = 12000; % Degen crystal length (um), 12 mm.

subplot(2,3,3); plot(lamp_mu,sinc((deltaK.*L)/2).^2,'Linewidth',1.5);

L = 9000; % Degen crystal length (um), 9 mm.

hold on; plot(lamp_mu,sinc((deltaK.*L)/2).^2,'r','Linewidth',1.5);

L = 6000;

hold on; plot(lamp_mu,sinc((deltaK.*L)/2).^2,'g','Linewidth',1.5);

L = 3000;

hold on; plot(lamp_mu,sinc((deltaK.*L)/2).^2,'k','Linewidth',1.5);
xlabel('Pump wavelength (um)')
ylabel('phase match efficiency')
title('Phase matching effiency for QPM period = 38.8603 um')

str = sprintf('Phase matching effiency for QPM period = %d',QPM);
title(str);
legend('L = 12 mm','L = 9 mm','L = 6 mm','L = 3 mm')


L = 7000; % Degen crystal length (um), 12 mm.

subplot(2,3,4); plot(deltaK,sinc((deltaK.*L)/2).^2,'Linewidth',1.5);


subplot(2,3,5)
plotyy(lamp_mu,deltaK.*10^3,lamp_mu,lams_mu);
grid on
xlabel('Pump wavelength (um)')
ylabel('{\Delta}k (1/mm)')

subplot(2,3,6)
plot(lamp_mu,deltaK.*10^3,'k','Linewidth',1.5);
grid on
xlabel('Pump wavelength (um)')
ylabel('{\Delta}k (1/mm)')



% If we assume pump is single frequency, but signal and idler multimode,
% then will need to calculate the deltaK value for fixed pump value,
% varying signal/idler wavelengths.

1;




















end


