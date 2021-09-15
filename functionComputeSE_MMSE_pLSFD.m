%-----  Compute SE with MMSE combining    -----%
%-----  Shuaifei Chen December 12 2019        -----%
% ==================================================================
%     This function outputs the SE.
% ==================================================================
%     This function was developed as a part of the paper:
% 
%     Shuaifei Chen, Jiayi Zhang, Emil Bj?rnson, Jing Zhang, 
%     and Bo Ai, "Structured Massive Access for Scalable 
%     Cell-Free Massive MIMO Systems," IEEE Journal on Selected 
%     Areas in Communications, vol. 39, no. 4, pp. 1086-1100, April 2021.
% 
%     License: This code is licensed under the GPLv2 license. 
%     If you in any way use this code for research that results 
%     in publications, please cite our paper as described above.
% ==================================================================


function [SE_MMSE] = functionComputeSE_MMSE_pLSFD(Hhat,H,R,B,MatA,tau_c,tau_p,nbrOfRealizations,N,K,L,p)

%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
    p = p*ones(K,1);
end

dim = length(B(1,1,1,1,:));

UEsServedByAPs = zeros(L,K);
for l = 1:L
    i = 1;
    for k = 1:K
        if MatA(l,k) == 1
            UEsServedByAPs(l,i) = k;
            i = i + 1;
        end
    end
end

UEsServedByCommonAPs = zeros(K,K);
MatB = MatA'*MatA;
for k = 1:K
    j = 1;
    for i = 1:K
        if MatB(k,i) ~= 0
            UEsServedByCommonAPs(k,j) = i;
            j = j + 1;
        end
    end
end

%Store identity matrices of different sizes
eyeN = eye(N);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MMSE = zeros(K,dim);

%Compute sum of all estimation error correlation matrices at every BS
C_tot = zeros(N,N,L,dim);

for l = 1:L
    Dl = UEsServedByAPs(l,:);
    Dl(Dl == 0) = [];
    for k = Dl
        for d = 1:dim
            C_tot(:,:,l,d) = C_tot(:,:,l,d) + p(k)*(R(:,:,l,k)-B(:,:,l,k,d));
        end
    end
end

%Diagonal matrix with transmit powers and its square root
Dp = diag(p);
Dp12 = diag(sqrt(p));



%Prepare to save simulation results
signal_MMSE = zeros(L,K,dim);
scaling_MMSE = zeros(L,K,dim);
G_MMSE_1 = zeros(L,L,K,dim);
G_MMSE_2 = zeros(L,L,K,dim);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    for d = 1:dim
        %Levels 1-3
        gp_MMSE = zeros(L,K,K);
        
        
        %Go through all APs
        for l = 1:L
            
            %Extract channel realizations from all UEs to AP l
            Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
            
            %Extract channel estimate realizations from all UEs to AP l
            Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:,d),[N K]);
            
            Dl = UEsServedByAPs(l,:);
            Dl(Dl == 0) = [];
            Hhatallj_l = Hhatallj(:,Dl);
            Dpl = diag(p(Dl));
            
            %Compute LP-MMSE combining
            V_MMSE = ((Hhatallj_l*Dpl*Hhatallj_l')+C_tot(:,:,l,d)+eyeN)\(Hhatallj*Dp);
            
            
            %Go through all UEs
            for k = 1:K
                
                
                
                %%MMSE combining
                v = V_MMSE(:,k);
                
                Dkl = diag(repmat(MatA(l,k),[1 N]));
                %Level 2 and Level 3
                signal_MMSE(l,k,d) = signal_MMSE(l,k,d) + (v'*Dkl*Hallj(:,k))/nbrOfRealizations;
                gp_MMSE(l,:,k) = gp_MMSE(l,:,k) + (v'*Dkl*Hallj)*Dp12;
                scaling_MMSE(l,k,d) = scaling_MMSE(l,k,d) + norm(Dkl*v).^2/nbrOfRealizations;
                
                
            end
            
        end
        
        for k = 1:K
            
            G_MMSE_1(:,:,k,d) = G_MMSE_1(:,:,k,d) + gp_MMSE(:,:,k)*gp_MMSE(:,:,k)'/nbrOfRealizations;
            Pk = UEsServedByCommonAPs(k,:);
            Pk(Pk == 0) = [];
            G_MMSE_2(:,:,k,d) = G_MMSE_2(:,:,k,d) + gp_MMSE(:,Pk,k)*gp_MMSE(:,Pk,k)'/nbrOfRealizations;
        end
        
    end
    
end

for d = 1:dim
    for k = 1:K
        
        %With LP-MMSE combining
        b = signal_MMSE(:,k,d);
        A = G_MMSE_1(:,:,k,d) + diag(scaling_MMSE(:,k,d)) - p(k)*(b*b');
        A4w = G_MMSE_2(:,:,k,d) + diag(scaling_MMSE(:,k,d));
        LSFDweight = pinv(A4w)*b;
        numerator = p(k)*(LSFDweight'*b)^2;
        denominator = LSFDweight'*A*LSFDweight;
        SE_MMSE(k,d) = prelogFactor*real(log2(1+numerator/denominator));
        
    end
end
