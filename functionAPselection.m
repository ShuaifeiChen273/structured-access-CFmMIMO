%-----          Initial Access Algorithm (v.1.0)           -----%
%-----             Shuaifei Chen December 12 2019                -----%
% ==================================================================
%     L*K matrix beta:          Large-Scale fading coefficients
%     L                         number of access points (APs)
%     K                         number of user equipments (UEs)
% ==================================================================
%     This function outputs the AP selection result.
% ==================================================================

function [MatA,BlackListOfUEs,lstOfWeakUEs] = functionAPselection(beta,tau_p)

[L,K] = size(beta);

%List of UEs which haven't finished AP selection.
lstOfUEs = 1:K;

%List of UEs which are served by one AP.
lstOfWeakUEs = [];

%Entry A(l,k) = 1 of matrix A means that UE k is served by AP l.
MatA = zeros(L,K);

SelectListOfUEs = zeros(L,K);

ServeListOfAPs = zeros(L,K);



%Entries in kth row of this matrix are indicies of the APs that will not be
%selected by the UE k.
BlackListOfUEs = zeros(L,K);

while lstOfUEs
    
    k = lstOfUEs(1);
    

    
    while 1
        
        BlackListOfUE = find(BlackListOfUEs(:,k)~=0);
        SelectListOfUE = find(MatA(:,k)~=0);
        
        %List of APs can offer service to UE (row vector)
        lstOfAPs = setdiff((1:L)', [BlackListOfUE; SelectListOfUE], 'stable');
        
        
        if isempty(lstOfAPs)
            break
        else
            
            if any(lstOfWeakUEs == k)
                MatA(lstOfAPs,k) = 1;
                break
            else
                
                [~,APIndex] = max(beta(lstOfAPs,k));
                l = lstOfAPs(APIndex);
                MatA(l,k) = 1;
                
                if sum(MatA(l,:)==1) > tau_p
                    
                    ServeListOfAP = find(MatA(l,:)~=0);
                    [~,UEIndex] = min(beta(l,ServeListOfAP));
                    kWeakest = ServeListOfAP(UEIndex);
                    
                    if any(lstOfWeakUEs == kWeakest)
                        kWeakest = k;
                    end
                    
                    BlackListOfUEs(l,kWeakest) = 1;
                    if sum(BlackListOfUEs(:,kWeakest) == 1) == L-1
                        lstOfWeakUEs = [lstOfWeakUEs kWeakest];
                    else
                        lstOfUEs = union(lstOfUEs, kWeakest, 'sorted');
                    end
                    
                    MatA(l,kWeakest) = 0;
                end
            end
        end
    end
    lstOfUEs(lstOfUEs == k) = [];
end



end





