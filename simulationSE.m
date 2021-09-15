%This Matlab script generates spectral efficiency results in the paper:
%
% Shuaifei Chen, Jiayi Zhang, Emil Bj?rnson, Jing Zhang, and Bo Ai,
% "Structured Massive Access for Scalable Cell-Free Massive MIMO Systems,"
% IEEE Journal on Selected Areas in Communications, vol. 39, no. 4,
% pp. 1086-1100, April 2021. https://ieeexplore.ieee.org/abstract/document/9174860)
%
%This is version 1.0 (Last edited: 2019-12-12)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.

close all;
clear;


%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
L = 100;

%Number of UEs in the network
nbrOfUEs = [100];

%Number of antennas per AP
N = 4;

%Compute number of pilots per coherence block
tau_p = 5;

%Length of the coherence block
tau_c = 200;

%Uplink transmit power per UE (mW)
pmax = 100;

SE_tot_all = zeros(max(nbrOfUEs),1,nbrOfSetups,length(nbrOfUEs));

for nbrOfUE = 1:length(nbrOfUEs)
    
    K = nbrOfUEs(nbrOfUE);
    
    
    %Prepare to save simulation results
    SE_tot = zeros(K,1,nbrOfSetups);
    
    n = 1;
    
    %% Go through all setups
    while n <= nbrOfSetups
        
        %Display simulation progress
        disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
        tic
        %Generate one setup with UEs at random locations
        [R,beta,UEpositions,APpositions] = functionGenerateR(L,K,N);
        
        [MatA,~,lstOfWeakUEs] = functionAPselection(beta,tau_p);
        
        [pilotIndex1,MatA_ug,delta,warnPilot,Groups] = functionUEgroup(beta,MatA,lstOfWeakUEs,tau_p);
        
        [p_ug] = functionPowerControl(beta,MatA_ug,pmax,1);
        
        if warnPilot == 1
            
            disp(['User-Grouping failure in setup ' num2str(n)]);
            
            continue
            
        else
            [Hhat,H,B] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex1,pmax);
            [SE] = functionComputeSE_MMSE_pLSFD(Hhat,H,R,B,MatA_ug,tau_c,tau_p,nbrOfRealizations,N,K,L,p_ug);
            
            %Save SE values
            SE_tot(:,:,n) = SE;
            
            n = n + 1;
            
            %Remove large matrices at the end of analyzing this setup
            clear H_AP R B Hhat;
            
        end
        toc
        disp([datetime]);
    end
    
    %Save SE values
    SE_tot_all(1:K,:,:,nbrOfUE) = SE_tot;
    
end

% save('data_SE_MMSE_K.mat','SE_tot_all');

