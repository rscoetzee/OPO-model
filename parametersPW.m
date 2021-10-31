function[h,deff,L,Lcav,c,lam_p,lam_s,lam_i,deltaK,w0x,w0y,w0sx,w0sy,w0ix,w0iy, ...
 AreaP,AreaS,AreaI,AreaCM,tp,PRF,vp,vs,vi,hplank,hbar,kp,ks,ki,n_p,n_s,n_i,e0, ...
 z,tRT,R1p,R2p,R1s,R2s,R1i,R2i,AR1p,AR2p,AR1s,AR2s,AR1i,AR2i,t,Pphase,  ...
 Sphase,Iphase,gg1,RT,w0cmx,w0cmy] = parametersPW()

lam_p = 1064.2*10^-9;                      % pump wavelength (m)
lam_s = 1856*10^-9;                      % signal wavelength (m)
invlam_i = (1./lam_p) - (1./lam_s);     
lam_i = 1/(invlam_i);                    % idler wavelength (m)
deltaK = 0;                              % wave vector mismatch (1/m)
c = 2.99792458*10^8;                          % speed of light. m/s
L = 0.007;                               % length of crystal (m).c
Lcav = 0.010;                            % length of cavity (m)   
deff = 10*10^-12; %(16.9*2/pi)*10^-12;               % m/V    %SNLO is a good source for these values.
e0 = 8.854*10^-12;                       % F/m
mu = 4*pi*10^-7;
w0x = 300*10^-6;                         % pump beam radius (m)
w0y = 330*10^-6;
w0cmx = 0.0300;
w0cmy = 0.0330;
w0sx = 300*10^-6;
w0sy = 330*10^-6;
w0ix = w0sx;
w0iy = w0sy;
tp = 10*10^-9;                           % pulse duration (s)
PRF = 100;                               % Pulse Repition rate H                           
vp = 2*pi*c/lam_p;                       % angular frequency
vs = 2*pi*c/lam_s;
vi = 2*pi*c/lam_i;
hplank = 6.62606957*10^-34;              % J.s
hbar = hplank/(2*pi);
AreaP = (pi*w0x*w0y)/2;
AreaS = (pi*w0sx*w0sy)/2;
AreaI = (pi*w0ix*w0iy)/2;
AreaCM = (pi*w0cmx*w0cmy)/2;
lamp_mu = 1.0642;
lams_mu = 1.856;
lami_mu = ((1./lamp_mu) - (1./lams_mu))^-1;

n_p = sqrt(4.59423 + (0.06206/((lamp_mu^2) - 0.04763)) ...     %Sellmeier KTP
    + (110.80672/((lamp_mu^2) - 86.12171)));           %Refractive indicies

n_s = sqrt(4.59423 + (0.06206/((lams_mu^2) - 0.04763)) ...
    + (110.80672/((lams_mu^2) - 86.12171)));

n_i = sqrt(4.59423 + (0.06206/((lami_mu^2) - 0.04763)) ...
    + (110.80672/((lami_mu^2) - 86.12171)));

kp = 2*pi*n_p/(lam_p);    % wave numbers
ks = 2*pi*n_s/(lam_s);
ki = 2*pi*n_i/(lam_i);        
h = L*0.01;    % Step size   
% h = 38.86*10^-6;
z = 0:h:L;
tRT = 2*(Lcav*n_i)/c;       % Calculate round trip time for given cavity L.

R1s = 0.999;     % Mirror reflectivites.
R2s = 0.5;
R1i = 0.999;   
R2i = 0.5;
R1p = 0.001; 
R2p = 0.001; 

AR1p = 0.086;                       % Reflectivity on 1st and 2nd side of crystal (if no AR-coating).
AR1s = 0.086;                       % Set to zero if crystal is coated.
AR1i = 0.086;                       % pump ~ 8.6 %, signal idler ~ 8.2 % reflectivity

AR2p = 0.086;
AR2s = 0.086;
AR2i = 0.086;

Pphase = 2*pi*rand;
Sphase = 2*pi*rand;
Iphase = 2*pi*rand;

RT = 300;
gg1 = 10;                %(number of slices per round trip)
dq = 1/gg1;               % 1/(number of slices per round trip)
t = (-RT/2:dq:(RT/2-dq))*tRT;
