
function[As,Ai,Ap] = rk4CK(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi)

N = floor(abs(L/h));

As = zeros(length(A0p),N);      % 2D matrix, rows -> Pulse values, coloumns -> N position.
Ai = zeros(length(A0p),N);
Ap = zeros(length(A0p),N);

Ap(:,1) = transpose(A0p);  
As(:,1) = transpose(A0s);
Ai(:,1) = transpose(A0i);

b21 =  1/5;
b31 =  3/40;  b32 =  9/40;
b41 =  3/10;  b42 = -9/10; b43 =  6/5;
b51 = -11/54; b52 =  5/2;  b53 = -70/27; b54 = 35/27;
b61 =  1631/55296; b62 = 175/512; b63 = 575/13824; b64 = 44275/110592; b65 = 253/4096;

a2 = 1/5; a3 = 3/10; a4 = 3/5; a5 = 1; a6 = 7/8;

c1 = 37/378; c2 = 0; c3 = 250/621; c4 = 125/594; c5 = 0; c6 = 512/1771;

for jj = 1:N+1 
        
        [k1p,k1s,k1i] = RHS(h,As(:,jj),Ap(:,jj),Ai(:,jj),z(jj),d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);        
                   
        [k2p,k2s,k2i] = RHS(h,As(:,jj)+b21*k1s,Ap(:,jj)+b21*k1p,Ai(:,jj)+b21*k1i,z(jj)+a2*h,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi); 
        
        [k3p,k3s,k3i] = RHS(h,As(:,jj)+b31*k1s+b32*k2s,Ap(:,jj)+b31*k1p+b32*k2p,Ai(:,jj)+b31*k1i+b32*k2i,z(jj)+a3*h, ... 
                            d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi); 
        
        [k4p,k4s,k4i] = RHS(h,As(:,jj)+b41*k1s+b42*k2s+b43*k3s,Ap(:,jj)+b41*k1p+b42*k2p+b43*k3p,Ai(:,jj)+b41*k1i+b42*k2i+b43*k3i, ... 
                            z(jj)+a4*h,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi); 
         
        [k5p,k5s,k5i] = RHS(h,As(:,jj)+b51*k1s+b52*k2s+b53*k3s+b54*k4s,Ap(:,jj)+b51*k1p+b52*k2p+b53*k3p+b54*k4p, ... 
                            Ai(:,jj)+b51*k1i+b52*k2i+b53*k3i+b54*k4i,z(jj)+a5*h,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi); 
        
        [k6p,k6s,k6i] = RHS(h,As(:,jj)+b61*k1s+b62*k2s+b63*k3s+b64*k4s+b65*k5s,Ap(:,jj)+b61*k1p+b62*k2p+b63*k3p+b64*k4p+b65*k5p, ...
                            Ai(:,jj)+b61*k1i+b62*k2i+b63*k3i+b64*k4i+b65*k5i,z(jj)+a6*h,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi); 
        
        Ap(:,jj+1) = Ap(:,jj) + c1*k1p + c2*k2p + c3*k3p + c4*k4p + c5*k5p + c6*k6p;
        As(:,jj+1) = As(:,jj) + c1*k1s + c2*k2s + c3*k3s + c4*k4s + c5*k5s + c6*k6s;
        Ai(:,jj+1) = Ai(:,jj) + c1*k1i + c2*k2i + c3*k3i + c4*k4i + c5*k5i + c6*k6i;

end


%         [k1p,k1s,k1i] = RHSDiff(h,FAp,FAs,FAi,FMixP,FMixS,FMixI,z,sx,sy,d);
%         
%         [k2p,k2s,k2i] = RHSDiff(h,FAp+k1p/2,FAs+k1s/2,FAi+k1i/2,FMixP+k1p/2,FMixS+k1s/2,FMixI+k1i/2,z+h/2,sx,sy,d);
%         
%         [k3p,k3s,k3i] = RHSDiff(h,FAp+k2p/2,FAs+k2s/2,FAi+k2i/2,FMixP+k2p/2,FMixS+k2s/2,FMixI+k2i/2,z+h/2,sx,sy,d);
%         
%         [k4p,k4s,k4i] = RHSDiff(h,FAp+k3p,FAs+k3s,FAi+k3i,FMixP+k3p,FMixS+k3s,FMixI+k3i,z+h,sx,sy,d);

% Algorithm for Diff model:

% (1) Start with Ap,As,Ai -> Calculate Mixing term(s)-> MixP, MixS, MixI.
% (2) FT(MixP,MixS,MixI) -> FMixP, FMixS, FmixI; FT(Ap,As,Ai) -> FAp,FAs,FAi
% (3) Use FMixP, FMixS, FmixI & FAp,FAs,FAi in RK4 to calculate FAp, etc at
% next step.
% (4) ifft(FAp(z+h)) -> Ap(z+h), etc.



% Old Rk4 for just mixing terms, i.e. no diffraction terms etc.

% N = floor(abs(L/h));
% 
% As = zeros(length(A0p),N);      % 2D matrix, rows -> Pulse values, coloumns -> N position.
% Ai = zeros(length(A0p),N);
% Ap = zeros(length(A0p),N);
% 
% Ap(:,1) = transpose(A0p);  
% As(:,1) = transpose(A0s);
% Ai(:,1) = transpose(A0i);
% 
% for j = 1:N+1
%         
%     [k1p,k1s,k1i] = RHS(h,As(:,j),Ap(:,j),Ai(:,j),z(j),d);  
%     
%     [k2p,k2s,k2i] = RHS(h,As(:,j)+k1s/2,Ap(:,j)+k1p/2,Ai(:,j)+k1i/2,z(j)+h/2,d);
%     
%     [k3p,k3s,k3i] = RHS(h,As(:,j)+k2s/2,Ap(:,j)+k2p/2,Ai(:,j)+k2i/2,z(j)+h/2,d);
%     
%     [k4p,k4s,k4i] = RHS(h,As(:,j)+k3s,Ap(:,j)+k3p,Ai(:,j)+k3i,z(j)+h,d);
%  
%     As(:,j+1) = As(:,j) + (k1s + 2*k2s +2*k3s + k4s)/6;
%     Ai(:,j+1) = Ai(:,j) + (k1i + 2*k2i +2*k3i + k4i)/6;
%     Ap(:,j+1) = Ap(:,j) + (k1p + 2*k2p +2*k3p + k4p)/6;   
%     
% end
% 

%% Note about fftshift etc:

        %  y=exp(-t.^2);
        %  figure;plot(abs(fftshift(fft(y))))
        %  figure;plot(phase(fftshift(fft(y))))
        %  figure;plot(phase(fftshift(fft(ifftshift(y)))))
        
