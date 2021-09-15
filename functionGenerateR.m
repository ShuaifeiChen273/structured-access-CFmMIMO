%-----  Generate Setup    -----%
%-----  Shuaifei Chen November 30 2019        -----%
% ==================================================================
%     This function outputs the spatial correlation matrix
%     and large-scale coefficients
% ==================================================================
%INPUT:
%L                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per AP
%
%OUTPUT:
%R_AP               = Matrix with dimension N x N x L x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%                     power
%gainOverNoisedB_AP = Matrix with dimension L x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by
%                     the noise power) between AP l and UE k in setup n
%UEpositions        = Vector of length K with UE positions, where the real
%                     part is the horizontal position and the imaginary
%                     part is the vertical position
%APpositions        = Vector of length L with the AP locations, measured in
%                     the same way as UEpositions
% ==================================================================
%-----  references    -----%
% ==================================================================
%     Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%     Competitive With MMSE Processing and Centralized Implementation,"
%     Submitted for publication, https://arxiv.org/abs/1903.10611
% ==================================================================

function [R,beta,UEpositions,APpositions] = functionGenerateR(L,K,N)
%% Define simulation setup

%Size of the coverage area (as a square with wrap-around)
squareLength = 500; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 5;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading
sigma_sf = 4;

%Decorrelation distance of the shadow fading
decorr = 9;

%Height difference between an AP and a UE
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Angular standard deviation around the nominal angle (measured in degrees)
ASDdeg = 15;

%Number of APs per dimension on the grid
nbrAPsPerDim = sqrt(L);

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

% Deploy APs on the grid
% locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
% locationsGridVertical = locationsGridHorizontal';
% APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

% Generate a random AP locations in the area
APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;

%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';

%Compute alternative AP locations by using wrap around
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Prepare to save results
gainOverNoisedB_AP = zeros(L,K);
R = zeros(N,N,L,K);
beta = zeros(L,K);


%% Generate Spatial Correlation Matrix and Large-Scale Coefficients
%Prepare to compute UE locations
UEpositions = zeros(K,1);

%Prepare to store shadowing correlation matrix
shadowCorrMatrix = sigma_sf^2*ones(K,K);

nbrOfUEs = 0;

%Add UEs until the each BS has tau_p UEs to serve
while nbrOfUEs < K
    
    %Generate a random UE location in the area
    UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
    
    %If this is not the first UE
    if nbrOfUEs>0
        
        %Compute distances from the new prospective UE to all other UEs
        shortestDistances = zeros(nbrOfUEs,1);
        
        for i = 1:nbrOfUEs
            shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
        end
        
        %Compute conditional mean and standard deviation necessary to
        %obtain the new shadow fading realizations, when the previous
        %UEs' shadow fading realization have already been generated.
        %This computation is based on Theorem 10.2 in "Fundamentals of
        %Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
        
    else %If this is the first UE
        
        %Add the UE and begin to store shadow fading correlation values
        newcolumn = [];
        
    end
    nbrOfUEs = nbrOfUEs + 1;
    
    %Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:nbrOfUEs-1,nbrOfUEs) = newcolumn;
    shadowCorrMatrix(nbrOfUEs,1:nbrOfUEs-1) = newcolumn';
    
    UEpositions(nbrOfUEs) = UEposition;
    
end

%Create correlated shadow fading realizations
shadowAPrealizations = sqrtm(shadowCorrMatrix)*randn(K,L);

%Go through all UEs
for k = 1:K
    
    %Compute distances assuming for the UE to the APs
    [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
    distances = sqrt(distanceVertical^2+distanceAPstoUE.^2);
    
    %Compute the channel gain divided by the noise power (in dB)
    gainOverNoisedB_AP(:,k) = constantTerm - alpha*log10(distances) + shadowAPrealizations(k,:)' - noiseVariancedBm;
    
    %Go through all APs
    for l = 1:L
        
        %Compute nominal angle between UE k and AP l
        angletoUE = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l)));
        
        %Generate normalized spatial correlation matrix using the local
        %scattering model
        R(:,:,l,k) = db2pow(gainOverNoisedB_AP(l,k))*functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
        
        %Compute the large-sacle coefficients
        beta(l,k) = trace(R(:,:,l,k))/N;
        
    end
%     k
end

end









