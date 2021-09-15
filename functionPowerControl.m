%-----  Fractional Power Control    -----%
%-----  Shuaifei Chen December 12 2019        -----%
% ==================================================================
%     This function outputs the power coefficients.
% ==================================================================
%     This function was developed as a part of the paper:
% 
%     S. Chen, J. Zhang, E. Bj?rnson, J. Zhang, and B. Ai, 
%     "Structured massive access for scalable cell-free 
%     massive MIMO systems," IEEE J. Sel. Areas Commun., 
%     Early Access, 2020.
% 
%     License: This code is licensed under the GPLv2 license. 
%     If you in any way use this code for research that results 
%     in publications, please cite our paper as described above.
% ==================================================================

function [p] = functionPowerControl(beta,MatA,pmax,theta)

%Large-scale fading coefficients that count
beta_count = beta.*MatA;

%Sum of the large-scale fading coefficients that count for all UEs
SumOfBetaPerUE = (sum(beta_count, 1)).^theta;

scaling = min(SumOfBetaPerUE);

p = pmax*scaling./SumOfBetaPerUE;
end





