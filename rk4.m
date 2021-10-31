function[As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi)

N = floor(abs(L/h));

As = zeros(length(A0p),N);      % 2D matrix, rows -> Pulse values, coloumns -> N position.
Ai = zeros(length(A0p),N);
Ap = zeros(length(A0p),N);

Ap(:,1) = transpose(A0p);  
As(:,1) = transpose(A0s);
Ai(:,1) = transpose(A0i);

for j = 1:N+1
    
    [k1p,k1s,k1i] = RHS(h,As(:,j),Ap(:,j),Ai(:,j),z(j),d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);  
    
    [k2p,k2s,k2i] = RHS(h,As(:,j)+k1s/2,Ap(:,j)+k1p/2,Ai(:,j)+k1i/2,z(j)+h/2,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);  
    
    [k3p,k3s,k3i] = RHS(h,As(:,j)+k2s/2,Ap(:,j)+k2p/2,Ai(:,j)+k2i/2,z(j)+h/2,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);  
    
    [k4p,k4s,k4i] = RHS(h,As(:,j)+k3s,Ap(:,j)+k3p,Ai(:,j)+k3i,z(j)+h,d,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi);  
 
    As(:,j+1) = As(:,j) + (k1s + 2*k2s +2*k3s + k4s)/6;
    Ai(:,j+1) = Ai(:,j) + (k1i + 2*k2i +2*k3i + k4i)/6;
    Ap(:,j+1) = Ap(:,j) + (k1p + 2*k2p +2*k3p + k4p)/6;   
    
end



