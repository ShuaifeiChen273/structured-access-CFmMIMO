%-----  User-Grouping Pilot Assignment (v.4.2)    -----%
%-----  Shuaifei Chen December 12 2019        -----%

% ==================================================================
%     Users served by the least likely clusters of APs share
%     one pilot sequence and barely cause inter-user interference.
%
%     Based on this, all users can be seperated into groups, and
%     each group is to be assigned one orthogonal pilot sequence.
% ==================================================================
%     M*K matrix beta                 all large-scale coefficients
%     M                               number of access points (APs)
%     K                               number of legitimate users (LUs)
%     delta                           threshold
%     M*K matrix MatrixA              AP-UE serving relation
%     K*K matrix MatrixB              UE-UE interfering relation
%     (K-1)*(K-1) matrix MatrixC      UE-UE interfering relation
%     K*K matrix Groups               UE groups where each row reprents a group
%     nbrOfGroups                     number of the LU groups
% ==================================================================
%     This function outputs the UE group and the number of it.
% ==================================================================

function [pilotIndex,MatA,delta,warnPilot,Groups] = functionUEgroup(beta,MatA,lstOfWeakUEs,tau_p)
%% Use bisection to find the appropriate AP selection threshold
%Number of groups,
nbrOfGroups = 0;
[L,K] = size(MatA);

%Threshold section in AP selection
deltaSect = [0.1 0.6];

%AP selection threshold
delta = 0.3;

flag = 0;

warnPilot = 0;

%% Prepare matrix B, and C
betaxMatA = beta.*MatA;
listBetaxMatA = reshape(betaxMatA,[],1);
[listBetaxMatA,~] = sort(listBetaxMatA,'descend');
listBetaxMatA(listBetaxMatA==0)=[];
% Indices(1:(L*K - length(listBetaxMatA))) = [];

while nbrOfGroups~=tau_p
    
    %Number of groups,
    nbrOfGroups = 0;
    
    %% user-group Pilot Assignment
    %beta-based pilot assignment which outputs the pilotIndex, number of
    %groups, and the groups
    
    %Each row is a UE group and the entries are the number of this group
    Groups = zeros(K, K);
    
    %Index that shows the pilot assigned to each UE
    pilotIndex = zeros(K,1);
    
    %List of the UE candidates which have not been grouped into a specific group
    %0 entry at the end of the list is a index will be used in the following progress to indicate all the UEs are grouped
    ListCand = 1:K;
    ListCand = [ListCand 0];
    
    MatA_ug = zeros(L,K);
    
    thresholdOfBeta = listBetaxMatA(ceil(delta*length(listBetaxMatA)));
    for k = lstOfWeakUEs
        MatA_ug(:,k) = MatA(:,k);
    end
    
    kp = 1;
    while listBetaxMatA(kp) >= thresholdOfBeta
        [i,j] = find(betaxMatA == listBetaxMatA(kp));
        MatA_ug(i,j) = 1;
        kp = kp+1;
    end
    %Symmetric matrix B shows the interference relationship among APs
    MatB = MatA_ug'*MatA_ug;
    %Compute matrix C
    MatC = zeros(K - 1,K - 1);
    
    for i = 1:K-1
        k = 1;
        for j = i+1:K
            if MatB(i,j) ==0
                MatC(i,k) = j;
                k = k + 1;
            end
        end
    end
    
    % ==================================================================
    %     Matrix B is a symmetric matrix, thus we only focus on upper
    %     triangular except the main diagonal.
    %
    %     B(i,j) = 0 means ith UE and jth UE are served by totally
    %     different clusters of APs; thus we assum that they can
    %     share same pilot sequence and barely interfere with each other.
    %
    %     Seperate UEs into different groups based on matrix C.
    % ==================================================================
    
    % Step1, select the groups that only have one mumber
    
    %Index that shows the sequence number of a group
    groupIndex = 1;
    
    %List that includes the groups that only have one mumber
    SoloGroup = [];
    
    %Flag used to select the single-member groups in Step1
    soloFlag = 0;
    
    %Repeat until all single-member groups in Step1 are selected
    while 1
        
        %Only select in the list of the UE candidates
        for i = setdiff(ListCand, [K 0], 'stable')
            
            %Corresponding row of this UE is a null vector means that the only member in its group is itself
            if MatC(i,:)==zeros(1,K-1)
                
                Groups(groupIndex, 1) = i;
                groupIndex = groupIndex + 1;
                
                %Add this group into to the list of single-member groups
                SoloGroup = [SoloGroup i];
                
            end
            
        end
        
        %Decide whether the single-member groups in Step1 are all selected
        if soloFlag == sum(SoloGroup)
            break
        end
        
        %Remove this UE from the rows whcih are in the list of the UE candidates
        for i = setdiff(ListCand, [K 0], 'stable')
            rowMatCLeft = setdiff(MatC(i,:), SoloGroup, 'stable');
            MatC(i,:) = [rowMatCLeft zeros(1,K-1-length(rowMatCLeft))];
        end
        
        
        soloFlag = sum(SoloGroup);
        
        %Remove the UEs constructing groups with themself from the list of the UE candidates
        ListCand = setdiff(ListCand, SoloGroup, 'stable');
    end
    
    % ==================================================================
    %     Since matrix C is £¨K-1£©-dimensional matrix, Kth LU candidate will be
    %     dicussed in the end if it still in the list of the UE candidates
    % ==================================================================
    
    % Step2, select the groups left
    
    %First member in the list of the UE candidates
    ListCand_1st = ListCand(1);
    
    %Repeat until all UEs are grouped
    while ListCand_1st > 0
        
        %No UE selects the Kth UE in its group
        if ListCand_1st == K
            
            %Kth UE constructs its group comprising only itself
            Groups(groupIndex, 1) = K;
            
            %Remove the Kth UE from the list of the UE candidates, and the progress ends
            ListCand(1)=[];
            %         disp(['Kth LU candidate is left in the end']);
            break
            
        else
            
            %First UE in the list of the UE candidates constructs its group with itself as the first member
            Groups(groupIndex, 1) = ListCand_1st;
            
            %Remove this UE from the list of the UE candidates
            ListCand(1) = [];
            
            %Extract the corresponding row of this UE in matrix C
            rowOfListCand = MatC(ListCand_1st, :);
            
            %         MatC(ListCand_1st, :) = zeros(1, K-1);
            
            %Remove this UE from each row corresponding to the UEs in the list of the UE candidates
            for i = setdiff(ListCand, [K 0], 'stable')
                rowMatCLeft = setdiff(MatC(i,:), ListCand_1st, 'stable');
                MatC(i,:) = [rowMatCLeft zeros(1,K-1-length(rowMatCLeft))];
            end
            
            %If the corresponding row of this UE is a null vector, which means
            %this UE can only construct its group with itself, then there is no
            %need to continue the following process
            if sum(rowOfListCand) == 0
                
                SoloGroup = [SoloGroup ListCand_1st];
                
            else
                
                %Enties in the row corresponding to this UE indicat the UEs
                %could be the members in the group constructed by this UE
                rowOfListCand_1st = rowOfListCand(1);
                
                %Index that shows the sequence number of a member in a group
                memberIndex = 2;
                
                %Repeat untill all UEs in this row are considered
                while rowOfListCand_1st > 0
                    
                    %Kth LU is the last member in this group
                    if rowOfListCand_1st == K
                        Groups(groupIndex, memberIndex) = K;
                        
                        %Remove the Kth UE from the list of the UE candidates
                        ListCand(ListCand == K) = [];
                        
                        
                        %Remove this UE from each row corresponding to the UEs
                        %in the list of the UE candidates
                        for i = setdiff(ListCand, [K 0], 'stable')
                            rowMatCLeft = setdiff(MatC(i,:), K, 'stable');
                            MatC(i,:) = [rowMatCLeft zeros(1,K-1-length(rowMatCLeft))];
                        end
                        
                        %Progress in this group ends
                        break
                        
                    else
                        
                        %First UE left in the row was added as a member of the group
                        Groups(groupIndex, memberIndex) = rowOfListCand_1st;
                        
                        %Remove this UE from the list of the UE candidates
                        ListCand(ListCand == rowOfListCand_1st) = [];
                        
                        %Remove the UEs in this row which could not be members
                        %of the UE which is just removed from the list of the
                        %UE candidates
                        candLeft = intersect(rowOfListCand, MatC(rowOfListCand_1st, :),  'stable');
                        
                        
                        %Remove this UE from each row corresponding to the UEs
                        %in the list of the UE candidates
                        for i = setdiff(ListCand, [K 0], 'stable')
                            rowMatCLeft = setdiff(MatC(i,:), rowOfListCand_1st, 'stable');
                            MatC(i,:) = [rowMatCLeft zeros(1,K-1-length(rowMatCLeft))];
                        end
                        
                        %If no members left in this row, then there is no need
                        %to continue the process of this group
                        if candLeft == 0
                            break
                        else
                            
                            rowOfListCand = [candLeft zeros(1,K-1-length(candLeft))];
                            
                        end
                        
                    end
                    
                    memberIndex = memberIndex+1;
                    
                    %Update the first candidate of this group
                    rowOfListCand_1st = rowOfListCand(1);
                    
                end
                
            end
            
        end
        
        groupIndex = groupIndex + 1;
        
        %Update the first candidate in the list of the UE candidates
        ListCand_1st = ListCand(1);
    end
    
    
    %Assign pilots to each UE and compute the number of groups
    for i = 1: K
        
        %Indicates that all the groups are considered
        if Groups(i,1)==0
            break
        end
        
        j = 1;
        
        %UEs in one group are assigned with the same pilot
        while Groups(i,j) ~= 0
            pilotIndex(Groups(i,j)) = i;
            j = j + 1;
        end
        
        %Compute the number of groups
        nbrOfGroups = nbrOfGroups + 1;
        
    end
    
    %% Use bisection to find the appropriate AP selection threshold (continued)
    
    Delta = delta;
    if nbrOfGroups == tau_p
        break
    elseif nbrOfGroups<tau_p
        deltaSect(1) = delta;
        delta = mean(deltaSect);
    else
        deltaSect(2) = delta;
        delta = mean(deltaSect);
    end
    
    if abs(delta - Delta) < 0.0001
        flag = flag + 1;
        if flag > 10
            nbrOfGroups = tau_p;
            warnPilot = 1;
        end
        
    end
    
end

end


