function[k1p,k1s,k1i] = RHS(h,As,Ap,Ai,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi) 


alpha_s = 0;  %mu*c/n_p; % 0.001; 
alpha_i = 0;  %mu*c/n_i; % 0.001; 
alpha_p = 0;  %mu*c/n_s; % 0.001; 

%% amplitude coupled part

k1p = d*h.*((1i.*deff.*vp./(c.*n_p)).*Ai.*As.*exp(-d*1i.*deltaK.*z) -(alpha_p/2).*Ap);

k1s = d*h.*((1i.*deff.*vs./(c.*n_s)).*Ap.*conj(Ai).*exp(d*1i.*deltaK.*z) -(alpha_s/2).*As); 

k1i = d*h.*((1i.*deff.*vi./(c.*n_i)).*Ap.*conj(As).*exp(d*1i.*deltaK.*z) -(alpha_i/2).*Ai);

% d = 1;
% 
% k1p = d*h.*((1i.*deff.*vp./(c.*n_p)).*Ai.*As.*exp(-d.*1i.*deltaK.*z) -(alpha_p/2).*Ap);
% 
% k1s = -d*h.*((1i.*deff.*vs./(c.*n_s)).*Ap.*conj(Ai).*exp(+d.*1i.*deltaK.*z) -(alpha_s/2).*As); 
% 
% k1i = d*h.*((1i.*deff.*vi./(c.*n_i)).*Ap.*conj(As).*exp(+d.*1i.*deltaK.*z) -(alpha_i/2).*Ai);



%%
