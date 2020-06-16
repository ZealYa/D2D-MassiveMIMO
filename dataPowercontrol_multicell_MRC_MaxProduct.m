clc;
clear all;%close all;

rng('shuffle');%first curvesshuffle

%% Model parameters
M = 100; %number of Antenna at the BS
K = 5; %number of cellular users per BS
NumberD2dTX = 10; %number of D2D transmitters
NumberD2dRX = NumberD2dTX; %number of D2D receivers
L = 9; % number of BS in the area
NumberTotalUsers =K*L+NumberD2dTX+NumberD2dRX;
tau = K+NumberD2dTX/2;
tauCohe = 200;
NumberSmallScaleIteration = 1000;
NumberLargeScaleIteration = 1000;
D2DPilotIndex = 1:5;
maxDistance = 0.5;%Km
D2dDistnce = 0.01;%D2D distance in meter
PowerCellular = 23;%dBm transmit power of cellular pilot and data transmission
PowerD2D = 23;%dBm transmit powet of D2D transmitters
CarrierFreq = 2000;%Mhz
d1 = 0.05;
d0 = 0.01;
Bandwidth = 20e6; % Bandwidth in Hz
NoiseFigure = 9;%
NoiseTemp = 290;%Kelvin
EffectiveNoise = -94;%10*log10(Bandwidth*1.38e-23*NoiseTemp)+ 30+ NoiseFigure;
Rho_CU = PowerCellular - EffectiveNoise;
Rho_DD = PowerD2D - EffectiveNoise;
Rho_CU_lin = 10.^(Rho_CU/10);
Rho_DD_lin = 10.^(Rho_DD/10);
hight = 1.65;%m
hightAP = 15;
LBS =  46.3+33.9*log10(CarrierFreq) - 13.82*log10(hightAP) - (1.1*log10(CarrierFreq)-.7)*hight + (1.56*log10(CarrierFreq)-0.8);
LDD =  46.3+33.9*log10(CarrierFreq) - 13.82*log10(hight) - (1.1*log10(CarrierFreq)-.7)*hight + (1.56*log10(CarrierFreq)-0.8);


%Set the length in meters of the total square area
squareLength = 1;

%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);


%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%% main
rateCellular = zeros(NumberLargeScaleIteration,K*L);
sumRate = zeros(NumberLargeScaleIteration,1);

for iterLarge= 1:NumberLargeScaleIteration
    disp(['The iteration number is ',num2str(iterLarge)])
    posD2dRX = zeros(1,NumberD2dRX);
    posXD2dTx = squareLength*rand(1,NumberD2dTX);
    posYD2dTx = squareLength*rand(1,NumberD2dTX);
    posD2dTX =posXD2dTx+ 1i*posYD2dTx;
    a= 2*pi*rand(1,NumberD2dRX);
    poskr = (D2dDistnce*cos(a)+posXD2dTx) + 1i*(D2dDistnce*sin(a) +posYD2dTx);
    for kk =1:NumberD2dRX
        
        while abs(poskr-BSpositions(5)) > abs(squareLength+1i*squareLength)
            a= 2*pi*rand(1);
            poskr(kk) = (D2dDistnce*cos(a)+posXD2dTx(kk)) + 1i*(D2dDistnce*sin(a) +posYD2dTx(kk));
        end
        posD2dRX(1,kk) =poskr(kk);
    end
    
    distancesD2dTxD2DRX = abs(repmat(posD2dTX,NumberD2dTX,1)-repmat(posD2dRX.',1,NumberD2dRX));
    BetaD2dTxD2dRx = zeros(size(distancesD2dTxD2DRX));
    BetaD2dTxD2dRx(distancesD2dTxD2DRX > d1)= -LDD-35*log10(distancesD2dTxD2DRX(distancesD2dTxD2DRX > d1));
    BetaD2dTxD2dRx(distancesD2dTxD2DRX > d0)= -LDD-15*log10(d1)-20*log10(distancesD2dTxD2DRX(distancesD2dTxD2DRX > d0));
    BetaD2dTxD2dRx(distancesD2dTxD2DRX <= d0)= -LDD-15*log10(d1)-20*log10(d0);
    BetaD2dTxD2dRxLinear = (10.^(BetaD2dTxD2dRx/10));
    
    
    %Prepare to put out UEs in the cells
    UEpositions = zeros(K,L);
    perBS = zeros(L,1);
    
    %Prepare to store average channel gain numbers (in dB)
    channelGaindB = zeros(K,L,L);
    BetaCellularBsLinear = zeros(K,L,L);
    for l = 1:L
        
        %Put out K UEs in the cell, uniformly at random. The procedure is
        %iterative since UEs that do not satisfy the minimum distance are
        %replaced with new UEs
        while perBS(l)<K
            
            %Put out new UEs
            UEremaining = K-perBS(l);
            posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
            posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
            posXY = posX + 1i*posY;
            
            %Store new UEs
            UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
            perBS(l) = perBS(l)+length(posXY);
            
        end
        
        distancesCellularD2DRX(:,:,l)=abs(repmat(UEpositions(:,l),[1 size(posD2dRX,2)]) - repmat(posD2dRX,[K 1]));
        %Go through all BSs
        for j = 1:L
            
            %Compute the distance from the UEs in cell l to BS j with a wrap
            %around topology, where the shortest distance between a UE and the
            %nine different locations of a BS is considered
            [distancesBSj] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
            
            %Compute average channel gain using the large-scale fading model in
            for user = 1:K
                if distancesBSj(user) > d1
                    channelGaindB(user,l,j)= -LBS-35*log10(distancesBSj(user));
                elseif distancesBSj(user)> d0
                    channelGaindB(user,l,j)= -LBS-15*log10(d1)-20*log10(distancesBSj(user));
                else
                    channelGaindB(user,l,j)=-LBS-15*log10(d1)-20*log10(d0);
                end
            end
            BetaCellularBsLinear(:,l,j) =  10.^(channelGaindB(:,l,j)./10);
        end
        
    end
    
    BetaCellularD2dRxLinearCell= zeros(K,NumberD2dTX,L);
    for l = 1:L
        BetaCellularD2dRx = zeros(size(distancesCellularD2DRX));
        BetaCellularD2dRx(distancesCellularD2DRX > d1)= -LDD-35*log10(distancesCellularD2DRX(distancesCellularD2DRX > d1));
        BetaCellularD2dRx(distancesCellularD2DRX > d0)= -LDD-15*log10(d1)-20*log10(distancesCellularD2DRX(distancesCellularD2DRX > d0));
        BetaCellularD2dRx(distancesCellularD2DRX <= d0)= -LDD-15*log10(d1)-20*log10(d0);
        BetaCellularD2dRxLinearCell(:,:,l) = (10.^(BetaCellularD2dRx(:,:,l)/10));
    end
    
    distancesD2dTxBS = zeros(NumberD2dTX,L);
    BetaD2dTxBS = zeros(NumberD2dTX,L);
    for uu1 = 1:NumberD2dTX
        for kk1= 1:L
            distancesD2dTxBS(uu1,kk1) = abs(posD2dTX(uu1)- BSpositions(kk1));
            if distancesD2dTxBS(uu1,kk1) > d1
                BetaD2dTxBS(uu1,kk1)= -LBS-35*log10(distancesD2dTxBS(uu1,kk1));
            elseif distancesD2dTxBS(uu1,kk1) > d0
                BetaD2dTxBS(uu1,kk1)= -LBS-15*log10(d1)-20*log10(distancesD2dTxBS(uu1,kk1));
            else
                BetaD2dTxBS(uu1,kk1)=-LBS-15*log10(d1)-20*log10(d0);
            end
        end
    end
    
    BetaD2dTxBSLinear = (10.^(BetaD2dTxBS/10));
    
    %random pilot allocation to users
    for ll = 1:NumberD2dTX
        if ll <= length(D2DPilotIndex)
            D2DPilot(ll) = D2DPilotIndex(ll);
        else
            D2DPilot(ll) = D2DPilotIndex(randperm(numel(D2DPilotIndex),1));
        end
    end
    
    gammaCuBs = zeros(K,L,L);
    for l = 1:L %received
        for j =1:L %location
            for k = 1:K
                gammaCuBs(k,j,l) = (tau*Rho_CU_lin*((BetaCellularBsLinear(k,j,l)).^2))/(1+Rho_CU_lin*tau*sum(BetaCellularBsLinear(k,:,l)));
            end
        end
    end
    
    
    gammaD2dTxBs = zeros(NumberD2dTX,L);
    for l = 1:L
        for ll = 1:NumberD2dTX
            gammaD2dTxBs(ll,l) = (tau*Rho_DD_lin*(BetaD2dTxBSLinear(ll,l).^2))/(1+Rho_DD_lin*tau*sum(BetaD2dTxBSLinear(D2DPilot == D2DPilot(ll),l)'));
        end
    end
    
    gammaCuD2dRx = zeros(K,NumberD2dRX,L);
    for d = 1:NumberD2dRX
        for l=1:L
            for k = 1:K
                gammaCuD2dRx(k,d,l) = (tau*Rho_CU_lin*(BetaCellularD2dRxLinearCell(k,d,l).^2))...
                    /(1+Rho_CU_lin*tau*sum(BetaCellularD2dRxLinearCell(k,d,:)));
            end
        end
    end
    gammaD2dTxRx = zeros(NumberD2dRX,NumberD2dTX);
    for ll = 1:NumberD2dRX
        for dd = 1:NumberD2dTX
            gammaD2dTxRx(ll,dd) = (tau*Rho_DD_lin*(BetaD2dTxD2dRxLinear(ll,dd).^2))/(1+Rho_DD_lin*tau*sum(BetaD2dTxD2dRxLinear(ll,D2DPilot == D2DPilot(dd))));
        end
    end
    
    
    cvx_solver('mosek');%SDPT3
    cvx_precision medium;
    cvx_begin gp
    cvx_quiet(true);
    variable lambda_mid(K,L);
    variable zeta(1,NumberD2dTX);
    variable p_d(1,NumberD2dTX) nonnegative;
    variable p_c(K,L) nonnegative;
    expression intereferenceCellularD2dRx(1,NumberD2dRX);
    expression cohInterferenceD2d(1,NumberD2dTX)
    expression intereferenceCellularBS(1,L)
    expression interferenceD2DTXBS(1,L)
    
    maximize prod(zeta)*prod(prod(lambda_mid))
    subject to
    for dd = 1:NumberD2dRX
        intereferenceCellularD2dRxpercell= cvx(zeros(L,K));
        for ll = 1:L
            for kk = 1:K
                intereferenceCellularD2dRxpercell(ll,kk) = p_c(kk,ll).*(BetaCellularD2dRxLinearCell(kk,dd,ll));
            end
        end
        intereferenceCellularD2dRx(dd) = sum(intereferenceCellularD2dRxpercell(:));
        cohInterferenceD2d(dd) = p_d(1,dd)*(BetaD2dTxD2dRxLinear(dd,dd) - gammaD2dTxRx(dd,dd));
    end
    
    
    
    for d =1:NumberD2dRX
        zeta(1,d)*(1+Rho_DD_lin*cohInterferenceD2d(d)+Rho_CU_lin*intereferenceCellularD2dRx(d) ...
            + Rho_DD_lin.*p_d(d).*(BetaD2dTxD2dRxLinear(d,d)-gammaD2dTxRx(d,d))) <= (Rho_DD_lin)*p_d(d)*(gammaD2dTxRx(d,d)) ;
    end
    
    for l = 1:L
        CellularBS = cvx(zeros(K,L));
        for ll = 1:L
            for k = 1:K
                CellularBS(k,ll) = p_c(k,ll).*(BetaCellularBsLinear(k,ll,l));
            end
        end
        intereferenceCellularBS(1,l) = sum(CellularBS(:));
        loopD2DTXBS = cvx(zeros(NumberD2dTX,L));
        for tt= 1:NumberD2dTX
            loopD2DTXBS(tt,:) = p_d(tt).*(BetaD2dTxBSLinear(tt,:));
        end
        interferenceD2DTXBS(1,l) = sum(loopD2DTXBS(:));
        
    end
    
    
    for l = 1:L
        for k =1:K
            lambda_mid(k,l)*(((1+Rho_CU_lin*intereferenceCellularBS(1,l)...
                + Rho_DD_lin*interferenceD2DTXBS(1,l)))+ Rho_CU_lin*M*sum(p_c(k,[1:l-1,l+1:end]).*gammaCuBs(k,[1:l-1,l+1:end],l)))...
                <= M*p_c(k,l)*gammaCuBs(k,l,l)*Rho_CU_lin;
        end
    end
    p_c <= 1;
    p_d <= 1;
    cvx_end
    cvx_status
    
    
    Ypprimecell= zeros(K,NumberD2dRX,NumberSmallScaleIteration,L);
    for l= 1:L
        Ypprimecell(:,:,:,l)=sqrt(Rho_CU_lin*tau)*((randn(K,NumberD2dRX,NumberSmallScaleIteration)+ 1i*randn(K,NumberD2dRX,NumberSmallScaleIteration))/sqrt(2)).*...
            repmat(sqrt(BetaCellularD2dRxLinearCell(:,:,l)),[1 1  NumberSmallScaleIteration])+ ...
            (randn(K,NumberD2dRX,NumberSmallScaleIteration)+...
            1i*randn(K,NumberD2dRX,NumberSmallScaleIteration))/sqrt(2);
    end
    estimatorFactorscell =zeros(K,NumberD2dRX,L);
    for ll=1:L
        for dd = 1:NumberD2dRX
            for kk = 1:K
                estimatorFactorscell(kk,dd,ll) = (sqrt(tau*Rho_CU_lin)*(BetaCellularD2dRxLinearCell(kk,dd,ll)))/(1+Rho_CU_lin*tau*sum(BetaCellularD2dRxLinearCell(kk,dd,:)));
            end
        end
    end
    for l = 1:L
        GCellularD2dRxhat(:,:,:,l) = Ypprimecell(:,:,:,l).*repmat((estimatorFactorscell(:,:,l)),[1 1 NumberSmallScaleIteration]);
    end
    GtotD2dTxD2dRx = (randn(NumberD2dRX,NumberD2dTX,NumberSmallScaleIteration)+1i*randn(NumberD2dRX,NumberD2dTX,NumberSmallScaleIteration))/sqrt(2);
    GD2dTxD2dRx = GtotD2dTxD2dRx.*repmat(sqrt(BetaD2dTxD2dRxLinear),[1 1 NumberSmallScaleIteration]);
    
    W = (randn(NumberD2dRX,NumberD2dTX,NumberSmallScaleIteration)+1i*randn(NumberD2dRX,NumberD2dTX,NumberSmallScaleIteration))/sqrt(2);
    
    Ypprime = sqrt(Rho_DD_lin*tau)*GD2dTxD2dRx + W;
    
    estimatorFactors = zeros(NumberD2dRX,NumberD2dTX);
    for ll = 1:NumberD2dRX
        for dd = 1:NumberD2dTX
            estimatorFactors(ll,dd) = (sqrt(Rho_DD_lin*tau)*BetaD2dTxD2dRxLinear(ll,dd)) /(1+sum(Rho_DD_lin*tau*BetaD2dTxD2dRxLinear(ll,D2DPilot == D2DPilot(dd))));
        end
    end
    GD2dTxD2dRxhat = Ypprime.*repmat((estimatorFactors),[1 1 NumberSmallScaleIteration]);
    
    firstTerm = zeros(NumberD2dRX,NumberSmallScaleIteration);
    secondTerm = zeros(NumberD2dRX,NumberSmallScaleIteration);
    thirdTerm = zeros(NumberD2dRX,NumberSmallScaleIteration);
    forthTerm = zeros(NumberD2dRX,NumberSmallScaleIteration);
    randomRateSmallScale = zeros(NumberD2dRX,NumberSmallScaleIteration);
    
    for iter = 1:NumberSmallScaleIteration
        for d = 1: NumberD2dRX
            firstTerm(d,iter) = p_d(d)*Rho_DD_lin*abs(GD2dTxD2dRxhat(d,d,iter))^2;
            secondTerm(d,iter) = p_d(d)*Rho_DD_lin*(BetaD2dTxD2dRxLinear(d,d)-gammaD2dTxRx(d,d));
            for ll=1:L
                thridTerm_interfrence(1,ll) = sum(p_c(:,ll).*(abs(GCellularD2dRxhat(:,d,iter,ll)).^2 + (BetaCellularD2dRxLinearCell(:,d,ll) - gammaCuD2dRx(:,d,ll))));
            end
            thirdTerm(d,iter) = Rho_CU_lin*sum(thridTerm_interfrence);
            forthTerm(d,iter) = Rho_DD_lin*sum(abs(GD2dTxD2dRxhat(d,[1:d-1,d+1:end],iter)).^2 + (BetaD2dTxD2dRxLinear(d,[1:d-1,d+1:end]) - gammaD2dTxRx(d,[1:d-1,d+1:end])));
            randomRateSmallScale(d,iter) = log2(1+(firstTerm(d,iter)./(1+secondTerm(d,iter)+thirdTerm(d,iter)+forthTerm(d,iter))));
        end
        
    end
    
    
    for l = 1:L
        CellularBS = zeros(K,L);
        for ll = 1:L
            for k = 1:K
                CellularBS(k,ll) = p_c(k,ll).*(BetaCellularBsLinear(k,ll,l));
            end
        end
        intereferenceCellularBS(1,l) = sum(CellularBS(:));
        loopD2DTXBS =zeros(NumberD2dTX,L);
        for tt= 1:NumberD2dTX
            loopD2DTXBS(tt,:) = p_d(tt).*(BetaD2dTxBSLinear(tt,:));
        end
        interferenceD2DTXBS(1,l) = sum(loopD2DTXBS(:));
    end
    
    
    SINR_c = zeros(K,L);
    for l = 1:L
        for k =1:K
            SINR_c(k,l)= (M*p_c(k,l)*gammaCuBs(k,l,l)*Rho_CU_lin)./ (((1+Rho_CU_lin*intereferenceCellularBS(1,l) + Rho_DD_lin*interferenceD2DTXBS(1,l)))+ Rho_CU_lin*M*sum(p_c(k,[1:l-1,l+1:end]).*gammaCuBs(k,[1:l-1,l+1:end],l)));
        end
    end
    rateCellular(iterLarge,:) = (log2(1+SINR_c(:)))';
    sumRate(iterLarge,:) = sum(sum(log2(1+SINR_c)))+sum(mean(randomRateSmallScale,2));
    
end


