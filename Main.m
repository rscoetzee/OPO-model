clc
close all force
close all
clear  
tic

[h,deff,L,Lcav,c,lam_p,lam_s,lam_i,deltaK,w0x,w0y,w0sx,w0sy,w0ix,w0iy, ...
AreaP,AreaS,AreaI,AreaCM,tp,PRF,vp,vs,vi,hplank,hbar,kp,ks,ki,n_p,n_s,n_i,e0, ...
z,tRT,R1p,R2p,R1s,R2s,R1i,R2i,AR1p,AR2p,AR1s,AR2s,AR1i,AR2i,t,Pphase,  ...
Sphase,Iphase,gg1,RT,w0cmx,w0cmy] = parametersPW();


% 8-02-2016, This model works really well. Quantum noise initiation done
% time vector problem solved. 

%% OPO part %%%%

%% Determining a slope efficiency %%    

% % This data set is crystal 5, degen, 12 mm long, have excel sheet.
% Pumpvalues = [1.39864 2.88255 3.87183 5.85038 7.82893 9.31285 11.2914 13.76459 15.2485 16.73242 18.71097];
% MeasuredOut = [0 0 0.05552 0.30538 0.80386 1.49205 2.89009 4.63315 5.78995 6.84641 8.27025];
% % % Calculated SNLO Fluences (J/cm^2) for crystal 5 (factor 2.* is because for degen,
% % % signal and idler have same energy so add them up for total output).
% SNLO = 2.*[0 0.0024 0.0075 0.0172 0.0263 0.0328 0.0411 0.051 0.0568 0.0624 0.0697]; 
% SNLO2 = 2.*[0 0 0.0000792 0.00384 0.0114 0.0170 0.0242 0.0330 0.0382 0.0433 0.0499];      

% Pumpvalues = [1 2 4.43 4.62 6.08 7.5 8.75 10.24 11.65 12.9 14.2 15.58 16.9 18.18 19.55 21.1 22.3 23.87];
% MeasuredOut = [0 0 0.kalle02 0.1 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7]; % degen VBG OPO measured 9/02/2016
% 
Pumpvalues = linspace(0,1,10);
MeasuredOut = linspace(0,1,10);

% Pumpvalues = 15;
% MeasuredOut = 1;

SimulatedOut = zeros(1,length(Pumpvalues));
OutOut = zeros(1,length(Pumpvalues));
SNLO = zeros(1,length(Pumpvalues));
SNLO2 = SNLO;

feature('jit', 'on'); feature('accel', 'on');

wait = waitbar(0,'Computing');

xmat = repmat(t'./2,1,length(Pumpvalues));
ymat = repmat(Pumpvalues,numel(t),1);
zmat = xmat;

xmat1 = xmat;
ymat1 = ymat;
zmat1 = zmat;

for ii = 1:length(Pumpvalues)
    
Eseed = 9.3330*10^-20; %J  % page 425 in Yariv, gives expression, for a single mode. (photon energy*bandwidth)
                             % Multimode - see page 432 Yariv.
                             % So assume bandwidth ~ cm^-1 ~ 3*10^8 Hz ~ 0.3 Ghz.
                             % E = hv = 9.3330e -20 J
                             % Therefore input power ~ 2.799*10^-9 W, 2.799 nW.

Ppseed = sqrt(4*log(2)/pi)*(Eseed/tp);
Ep(ii) = Pumpvalues(ii)/1000;                 % pulse energy J
Pp(ii) = sqrt(4*log(2)/pi)*(Ep(ii)/tp);


initial_pumpenvelope = sqrt(Pp(ii)).*exp(-1.*(t.^2)./(tp/sqrt(2*log(2))).^2);
Seed_env = sqrt(Ppseed).*exp(-1.*(t.^2)./(tp/sqrt(2*log(2))).^2);

trapz(t,initial_pumpenvelope.^2);         % integrate under envelope = pulse energy.
trapz(t,initial_pumpenvelope.^2).*PRF;    % integrate under envelope * PRF = pulse power.

Signal_power_out = zeros(length(t),1);
Idler_power_out = zeros(length(t),1);
Total_output = zeros(length(t),1);
Signal_power_out_left = zeros(length(t),1);
Idler_power_out_left = zeros(length(t),1);
DepPump_left = zeros(length(t),1);
OutputPulse = zeros(length(t),1);
DepPump = zeros(length(t),1);
Theta1 = zeros(length(t),1);
Theta2 = zeros(length(t),1);
A0p = zeros(1,gg1);       
A0s = zeros(1,gg1);      
A0i =  zeros(1,gg1);       
A0pr = zeros(1,gg1)'; 
A0ptrans = zeros(1,gg1)';
A0sr = A0pr; A0ir = A0pr;
A0strans = A0ptrans; A0itrans = A0ptrans;

d = 1;
Noshots = 1;
Outvalue = zeros(1,Noshots);

E0 = zeros(1,length(t));
E0s = E0; E0i = E0;
    
    for SN = 1:Noshots    % Shot number

    E0  = sqrt(2.*(initial_pumpenvelope.^2)./(AreaP.*e0.*n_p.*c)).*exp(-1i.*Pphase);    % Field amplitude
    E0s = sqrt(2.*(Seed_env.^2)./(AreaS.*e0.*n_s.*c)).*exp(-1i.*Sphase);
    E0i = sqrt(2.*(Seed_env.^2)./(AreaI.*e0.*n_i.*c)).*exp(-1i.*Iphase);

        for kk = 0:RT-1 

        A0p = sqrt((1-AR1p)).*sqrt((1-R1p)).*(E0(1+kk*gg1:kk*gg1+gg1))  + sqrt((1-AR1p)).*transpose(A0pr);   
        A0s = sqrt((1-AR1s)).*sqrt((1-R1s)).*(E0s(1+kk*gg1:kk*gg1+gg1)) + sqrt((1-AR1s)).*transpose(A0sr); 
        A0i = sqrt((1-AR1i)).*sqrt((1-R1i)).*(E0i(1+kk*gg1:kk*gg1+gg1)) + sqrt((1-AR1i)).*transpose(sort(A0ir,'descend')); 

        % Forward prop
        [As,Ai,Ap] = rk4CK(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);
%         [As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);

        Theta2(1+kk*gg1:kk*gg1+gg1)    = phase(Ap(:,end))-phase(As(:,end))-phase(Ai(:,end));
        pumpphase(1+kk*gg1:kk*gg1+gg1) = phase(Ap(:,end));
        sigphase(1+kk*gg1:kk*gg1+gg1)  = phase(As(:,end));
        idphase(1+kk*gg1:kk*gg1+gg1)   = phase(Ai(:,end));

        Signal_power_out(1+kk*gg1:kk*gg1+gg1) = ((1-AR2s)).*(1-R2s).*((abs(As(:,end)).^2)).*(AreaS.*e0.*n_s.*c./2);
        Idler_power_out (1+kk*gg1:kk*gg1+gg1) = ((1-AR2i)).*(1-R2i).*((abs(Ai(:,end)).^2)).*(AreaI.*e0.*n_i.*c./2);
        DepPump         (1+kk*gg1:kk*gg1+gg1) = ((1-AR2p)).*(1-R2p).*((abs(Ap(:,end)).^2)).*(AreaP.*e0.*n_p.*c./2);

        Total_output(1+kk*gg1:kk*gg1+gg1) = Signal_power_out(1+kk*gg1:kk*gg1+gg1) + Idler_power_out(1+kk*gg1:kk*gg1+gg1);

        A0i = ((1-AR2i)).*(sqrt(R2i).*abs(Ai(:,end))).*exp(-1i.*(phase(Ai(:,end)))); 
        A0p = ((1-AR2p)).*(sqrt(R2p).*abs(Ap(:,end))).*exp(-1i.*(phase(Ap(:,end)))); 
        A0s = ((1-AR2s)).*(sqrt(R2s).*abs(As(:,end))).*exp(-1i.*(phase(As(:,end))));  
        
        d = -d; 
        z = sort(z,'descend');

        % Backward Prop
        [As,Ai,Ap] = rk4CK(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);
%         [As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);
        
        A0pr = sqrt((1-AR1p)).*(sqrt(R1p).*abs(Ap(:,end))).*exp(-1i.*(phase(Ap(:,end))));                             
        A0sr = sqrt((1-AR1s)).*(sqrt(R1s).*abs(As(:,end))).*exp(-1i.*(phase(As(:,end))));  
        A0ir = sqrt((1-AR1i)).*(sqrt(R1i).*abs(Ai(:,end))).*exp(-1i.*(phase(Ai(:,end))));  
        
        Signal_power_out_left(1+kk*gg1:kk*gg1+gg1) = ((1-AR1s)).*(1-R1s).*((abs(As(:,end)).^2)).*(AreaS.*e0.*n_s.*c./2);

        d = -d;
        z = sort(z,'ascend');

        Theta1(1+kk*gg1:kk*gg1+gg1) = phase(Ap(:,end))- phase(As(:,end))-phase(Ai(:,end));

        end

      OutputPulse = abs(Total_output);  % Power right side (signal+idler).
      
      OutputPulse_left = abs(Total_output);
          
    end    % for shot number > 1.
    
    waitbar(ii/length(Pumpvalues));
    
    Es = trapz(t,OutputPulse);
    Outvalue(SN) = Es*1000; 
    OutOut(ii) = mean(Outvalue);
    
    Es_left = trapz(t,OutputPulse_left);
    Outvalue_left(SN) = Es_left*1000;
    OutOut_left(ii) = mean(Outvalue_left);
    
    Undep(ii) = 1000*trapz(t,DepPump);
    
    zmat(:,ii) = OutputPulse;
    zmat1(:,ii) = DepPump;
    
end

close(wait)
 
subplot(1,2,1); plot(Pumpvalues,OutOut_left,'b+');

subplot(1,2,2); plot(t,Signal_power_out_left);

1;

figure;
plot3(xmat,ymat,zmat,'r','Linewidth',1.2);
xlabel('time(s)')
ylabel('Pump energy (mJ)')
zlabel('Output peak power (arb)')
grid on

figure;
plot3(xmat1,ymat1,zmat1,'b','Linewidth',1.2);
xlabel('time(s)')
ylabel('Pump energy (mJ)')
zlabel('Peak power (arb)')
grid on
        
Es = trapz(t,OutputPulse);  
e = 10;        % B = find(A <= (max(A)/2)+e & (max(A)/2)-e <= A)
MaxO = find(OutputPulse >= max(OutputPulse)/2);
% FWHM = t(MaxO(end))-t(MaxO(1));

Output_envelope = OutputPulse;

figure;

plot(t.*10^9,abs(Signal_power_out)/max(abs(initial_pumpenvelope.^2)),'g','Linewidth',1.5)
hold on
plot(t.*10^9,DepPump/max(abs(initial_pumpenvelope.^2)),'b-.','Linewidth',1.5) 
hold on
plot(t.*10^9,abs(initial_pumpenvelope).^2/max(abs(initial_pumpenvelope.^2)),'r-','Linewidth',1.5)
hold on
plot(t.*10^9,abs(Idler_power_out)/max(abs(initial_pumpenvelope.^2)),'k','Linewidth',1.5)
hold on
plot(t.*10^9,abs(Total_output)/max(abs(initial_pumpenvelope.^2)),'m:','Linewidth',1.5)
xlabel('Time (ns)');
ylabel('Normalized Instantaneous Power (W)');
legend('Signal','DepPump','Pump','Idler')
grid on
%axis([t(1).*10^9 t(end).*10^9 0 1.2.*10^9])

figure1 = figure('Position', [100, 100, 1500, 800]);
figure1;
grid on

subplot(2,4,1)
plot(t,initial_pumpenvelope.^2,'r','Linewidth',2);
xlabel('time (s)');
ylabel('Instantaneous Power (W)');
grid on

subplot(2,4,2)
hold on    
plot(Pumpvalues,OutOut./(1000*AreaCM),'g*','Linewidth',2);
hold on
plot(Pumpvalues,MeasuredOut./(1000*AreaCM),'ro','Linewidth',2);
hold on
plot(Pumpvalues,SNLO,'k+','Linewidth',2);
xlabel('Pump energy (mJ)');
ylabel('Fluence (J/cm^2)');
legend('Num','Exp','SNLO','Location','NorthWest');
grid on;

subplot(2,4,3)
plot(Theta1,'b');
hold on
plot(Theta2,'r');
xlabel('Index number');
ylabel('Phase (rad)');
legend('p-s-i','p-s-i (right mirror)');
axis([0 length(Theta1) -3*pi 3*pi])

subplot(2,4,4)
plot(abs(Total_output),'k:','Linewidth',2);
hold on
plot(abs(DepPump),'b-.','Linewidth',2); 
hold on
plot(abs(initial_pumpenvelope.^2),'r','Linewidth',2);
hold on
xlabel('Index number');   % Convert this to time (ns) etc.
ylabel('Instantaneous Power (W)');
legend('Signal','DepPump','Pump','Location','NorthWest')

% SNLOout = 1000.*(0.5).*(SNLO.*pi*((0.2360/2)^2));  %mJ
SNLOout = 1000.*(0.5).*(SNLO.*pi*((w0cmx*w0cmy))); 
% SNLOout = 1000.*(SNLO.*pi*((w0cmx)^2));
% SNLOout = zeros(1,length(Pumpvalues));
SNLOout2 = 1000.*(0.5).*(SNLO2.*pi*((w0cmx*w0cmy)));

subplot(2,4,5);
plot(Pumpvalues,MeasuredOut,'ro','Linewidth',2);
hold on    
plot(Pumpvalues,OutOut,'g*','Linewidth',2);
hold on
plot(Pumpvalues,SNLOout,'k+','Linewidth',2);
hold on
plot(Pumpvalues,Undep,'bo','Linewidth',2);
xlabel('Pump energy (mJ)');
ylabel('Total Output energy (mJ)');
legend('Exp','Num','SNLO','Undep','Location','NorthWest');
grid on;

subplot(2,4,6)
plot(pumpphase,'r');
hold on
plot(sigphase,'g')
hold on
plot(idphase,'k')
xlabel('index number');
ylabel('Phase(rad)');
legend('Pump','Signal','Idler');
axis([0 length(Theta1) -3*pi 3*pi])

subplot(2,4,7)
plot(Pumpvalues,(100*MeasuredOut./Pumpvalues),'r+','Linewidth',2);
xlabel('Pump energy (mJ)');
ylabel('Conversion Efficiency (%)');
hold on
plot(Pumpvalues,(100*OutOut./Pumpvalues),'g+','Linewidth',2);
hold on
plot(Pumpvalues,(100*SNLOout./Pumpvalues),'k+','Linewidth',2);
legend('Exp','Num','SNLO','Location','NorthWest');

ffpumpin = padarray(initial_pumpenvelope.^2,[0 1000000]);
ffdep = padarray(DepPump',[0 1000000]);
ffsig = padarray(Signal_power_out',[0 1000000]);
ffi = padarray(Idler_power_out',[0 1000000]);

t1 = padarray(t,[0 1000000]);
dt = abs(t1(round(length(t1)/2))-(t1(round(1+length(t1)/2))));
nfft = length(t1);
%frequency spectrum
fs = 1/dt;
df = fs/(nfft);
f = -fs/2:df:(fs/2)-df;

wc = max(f);            % Carrier frequency of pulse 
subplot(2,4,8)
% FF1 = fftshift(abs(fft(ffpumpin)));
% plot(f,FF1,'r');
hold on;
FF2 = fftshift(abs(fft(ffsig)));
plot(f,FF2,'g','Linewidth',1.2);
FF3 = fftshift(abs(fft(ffi)));
hold on;
plot(f,FF3,'k','Linewidth',1.2);
FF4 = fftshift(abs(fft(ffdep)));
hold on;
plot(f,FF4,'b','Linewidth',1.2);
axis([-6*10^8 6*10^8 0 1.1*max(FF4)]);
xlabel('Frequency (Hz)','Linewidth',2);

figure;
plot(Pumpvalues,(100*MeasuredOut./Pumpvalues),'r+','Linewidth',2);
xlabel('Pump energy (mJ)');
ylabel('Conversion Efficiency (%)');
hold on
plot(Pumpvalues,(100*OutOut./Pumpvalues),'g*','Linewidth',2);
hold on
plot(Pumpvalues,(100*SNLOout2./Pumpvalues),'k+','Linewidth',2);
legend('Exp','Num','SNLO2','Location','NorthWest');

figure;
plot(Pumpvalues,MeasuredOut,'ro','Linewidth',2);
hold on    
plot(Pumpvalues,OutOut,'g*','Linewidth',2);
hold on
plot(Pumpvalues,SNLOout2,'k+','Linewidth',2);
hold on
plot(Pumpvalues,Undep,'bo','Linewidth',2);
xlabel('Pump energy (mJ)');
ylabel('Total Output energy (mJ)');
legend('Exp','Num','SNLO2','Undep','Location','NorthWest');
grid on;

% plot(Outvalue,'r+'); 
% xlabel('Shot number');
% ylabel('Energy output');

% %%
% 
% dt = abs(t(2)-t(1));
% nfft = length(t);
% %frequency spectrum
% fs = 1/dt;
% df = fs/(nfft);
% f = -fs/2:df:(fs/2)-df;
% 
% wc = max(f);            % Carrier frequency of pulse 
% figure;
% FF1 = fftshift(abs(fft(initial_pumpenvelope.^2)));
% plot(f,FF1,'r');
% xlabel('Frequency (Hz)');
% hold on;
% FF2 = fftshift(abs(fft(Signal_power_out)));
% plot(f,FF2,'g');
% FF3 = fftshift(abs(fft(Idler_power_out)));
% hold on;
% plot(f,FF3,'k');
% FF4 = fftshift(abs(fft(DepPump)));
% hold on;
% plot(f,FF4,'b');
% axis([-6*10^8 6*10^8 0 1.1*max(FF1)]);


% figure; plot(t,initial_pumpenvelope.^2)

%%

%% Slope Efficiency end %%

%           if ~isempty(find(isnan(Total_output)))
%               1;
%           end

toc
