function [SNR,XYbc,XYbcnr,CI99] = Function_DSW(XY,lbws,ck)
%% Statistical method for baseline correction and noise removal based on double moving windows
% Statistical methods were used for identification of peaks, estimation of noise, and calculation of SBR
% For baseline, two potential baselines were generated, potential baselines were selected adaptively based on existance of peak
% AND: No peak: small window; With peak: large window

% 2022/12/01 By Zijiang Yang
% v2 note: moving averageing index - the flunctuation of baseline reflects
% Baseline is not smooth due to uncertainty
% There is a trade-off between smoothness and uncertainty
% 2022/12/03 Adding factors for skewness of the corrected noise distribution (the corrected noise is not normally distributed but skewed)
% 2022/12/03 Adding calculation code for estimating fup and flw
% 2022/12/13 Skewness correction factor was added for better estimating of sigma_ns

% 2023/01/09 v4 Reorder XY first - ascending order

%% Iput data and parameters
% ck      = 12; % Index of spectrum for checking, range is 400-4000, if may change if spectra are of other ranges
% XY      = [Wavenumber Intensity-1, intensity-2...] Spectrum data, array
% lbws    = 25; % % small window size and stepsize for baseline correction, this is also used for calculating envelope

%% Preperation
% Make sure X is of ascending order
if XY(1,1) > XY(end,1)
    XY = flip(XY);
end

X  = XY(:,1);
Y  = XY(:,2:end);

% Input parameters - may be adjusted for optimazation
lbws1  = lbws;     % small window size and stepsize for baseline correction, this is also used for calculating envelope
lbws2  = 5*lbws1;  % large window size and stepsize for baseline correction (factor of 5 was empirical, at prilimilary trials, we found 25 was good. 5*25=125 ~ 128 was from Renner et al., 2016, then factor of 5 was used to lower the number of parameters)

ex     = 1.3143;% expantion factor. ex=1 for 95%CI; ex=1.3143 for 99%CI, ex=1.4321 for 99.5%CI, ex=1.6791 for 99.9%. Reason using 95% CI to estimate others for more stable convergent of 95%CI
sn     = 4;     % sensitivity of peak detection (if we have more than 4 points outside noise range in a window (@lbws1) at the same time, then it is considered as a peak), and 4 is empirical; larger sn may misidentify peak as noise, small sn may misidentify noise as peak, so it is a trade off.
pr     = lbws1; % smoothing radium when combine the two potential baselines, let it = lbws1 (small window) for fewer parameters
fs     = 1;     % smoothing factor - to presenting sharp changes when combining two baselines: small fs: sensitive to narrow BC hill, but may made BC under peak bend, =1 means lbws1 was used

% Input parameters - based on mathmatical deriviation   
% Skewness correction factors: parameters that was estimated by simulation (at the end of script)
% Rn is the estimated 95% CI length of noise by assuming normal distribution but the real distribution of baseline corrected pure noise was not positively distributed
fup   = 3.777  * lbws1^-0.5296 + 0.8164;   % 1.5 % CIup - CImd = 0.5*Rn*fup, and fup = f1(lbws1) by fitting simulated data
flw   = 0.9594 * lbws1^-0.6008 + 0.9407;   % 0.7 % CImd - CIlw = 0.5*Rn*fmd2lw, and fmd2lw = f2(lbws1) by fitting simulated data
fsk   = 0.2246 * lbws1^-0.1126 + 0.2369;   % ~0.4 at lbws = 25; empirical

%% Detction of peak and noise
% Nagetive signal correction
Y     = Y - min(Y);

% Estimation of noise range
% Rn: most likely noise range; XYbn: baseline corrected & noise removed spectrum based on small window size
[Rn,~,XYbn1] = Function_PeakIdentify_f4(XY,lbws1,ck);

% Double baseline 1: small window size for noise
Ybc1   = msbackadj(X,Y,'WindowSize',lbws1,'StepSize',lbws1,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',0,'ShowPlot',0); 
BC1    = Y - Ybc1 + flw*0.5*Rn; % baseline corrected by most likely noise range, 0.7 is the distance between Q0.005 and Q50

% Double baseline 2: large window size for peaks
Ybc2   = msbackadj(X,Y,'WindowSize',lbws2,'StepSize',lbws2,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',0,'ShowPlot',0); 
BC2    = Y - Ybc2 + flw*0.5*Rn; % baseline corrected by most likely noise range

% Selection index based on candidate peaks
index0  = XYbn1(:,2:end) > 0;
index1  = movmean(1*index0,fs*pr / abs(X(2) - X(1))); % for index ~= 1, 1.5 factor is for smoothing
index2  = index1 > sn/pr;

% Selected baselines
index   = movmean(1*index2,lbws2 / abs(X(2) - X(1))); % smooth the indext to have smoothed baseline
BCnp    = BC1 .* (1 - index); % Selected baseline at non-peak region

BCpk1   = BC1 .* index; % Small BC at peak region
BCpk2   = BC2 .* index; % Large BC at peak region
BCpk    = NaN(height(BCpk2),width(BCpk1));

for i = 1:1:width(BCpk1)
    BCpk(:,i) = min(BCpk1(:,i),BCpk2(:,i));
end

BCfinal = BCnp+BCpk;

% Baseline corrected spectrum
Ybc     = Y - BCfinal;
XYbc    = [X Ybc];

%% Noise removal and SNR estimation
Ybcnr   = ones(height(Ybc),width(Ybc));

for i = 1:1:width(Ybc)
    ybc = Ybc(:,i);
    ybc(ybc < 0+ex*fup*0.5*Rn(i))=0; % 95% CI * ex; ex is CI expantion factor for other CI's
    Ybcnr(:,i) = ybc;
end

XYbcnr = [X Ybcnr]; % noise removed & baseline corrected spectrum - this is not used for now

% SNR estimation
sigma = fsk*(Rn.* (fup + flw) ./ (1.96*2)); % Here, Rn* (fup + flw) is 95% CI, but since it is not symmetric, cant use the relation between CI and sigma to estimate sigma directly
Hpeak = max(Ybcnr);
snr   = Hpeak ./ sigma; % fup - flw is 95% CI, and it equals to 4*σ; max(Ybcnr) = max peak height
SNR   = [Hpeak; sigma; snr];

% CI data
Y995  = BCfinal + ex*fup*0.5*Rn;
Y005  = BCfinal - ex*flw*0.5*Rn;
CI99  = [Y005 Y995];

%% Plot for checking
if ck > 0
    figure
    
    % Spectrum with candidate baselines & detected peak
    subplot(3,1,1)
    area(X,2*max(Y(:,ck))*index2(:,ck),'FaceColor','#CCECFF','EdgeColor','none') % candidate peak
    hold on
    plot(X,Y(:,ck),'lineWidth',1.2,'Color','#6699FF') % original spectra
    plot(X,BCfinal(:,ck),'lineWidth',1.2,'Color','#FF7C80') % estimated BC
    
    plot(X,BCfinal(:,ck) + ex*fup*0.5*Rn(ck),'lineWidth',0.5,'Color','#FF7C80') % Upper boundary of BC
    plot(X,BCfinal(:,ck) - ex*flw*0.5*Rn(ck),'lineWidth',0.5,'Color','#FF7C80') % Lower boundary of BC
    
    title('Spectrum with estimated baseline (99% CI) & potential peak')
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity/Absorption')
    xlim([400 4000])
    ylim([min(Y(:,ck)) 1.2*max(Y(:,ck))])
    set(gca, 'XDir','reverse')
    
    % Baseline corrected spectrum
    subplot(3,1,2)
    plot(X,Ybc(:,ck),'lineWidth',1.2,'Color','#00CC95')
    hold on
    plot([min(X) max(X)],[0 0],'lineWidth',1.2,'Color','#FF7C80')
    plot([min(X) max(X)],[0+ex*fup*0.5*Rn(ck) 0+ex*fup*0.5*Rn(ck)],'lineWidth',0.5,'Color','#FF7C80')
    plot([min(X) max(X)],[0-ex*flw*0.5*Rn(ck) 0-ex*flw*0.5*Rn(ck)],'lineWidth',0.5,'Color','#FF7C80') 
    title('Baseline (99% CI) corrected spectrum ')
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity/Absorption')
    xlim([400 4000])
    set(gca, 'XDir','reverse')
    
    subplot(3,1,3)
    plot(X,Ybcnr(:,ck),'lineWidth',1.2,'Color','#00CC99')
    hold on
    title('Baseline corrected & noise removed spectrum')
    xlim([400 4000])
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Probability density')
    set(gca, 'XDir','reverse')
    
    set(gcf, 'Position', [1            1       1/2*2194.3       1234.3])
end

end


%% Below is used for estimating pdf (and fup and flw) of Rn
% PeakIdentify_f4.m is needed

% %% Prepare Workspace
% clc, clear, close all
% 
% tic
% %% Parameters
% lbws  = (10:1:100)';   % window size and stepsize for baseline correction
% n     = 10;
% 
% X  = (400:1:4000)';
% 
% mCIlw = NaN(length(lbws),1);
% mCIup = NaN(length(lbws),1);
% mCImd = NaN(length(lbws),1);
% sYbc  = NaN(length(lbws),1);
% 
% fup   = NaN(length(lbws),1);
% flw   = NaN(length(lbws),1);
% 
% for i = 1:1:length(lbws)
% 
%     tic
%     lbws1  = lbws(i);
% 
%     % Simulated dataset
%     XY = [X normrnd(0,1,height(X),n)];
%     X  = XY(:,1);
%     Y  = XY(:,2:end);
%                 
%     %% Simulation
%         
%     [Rn1,XYbc,~] = PeakIdentify_f4(XY,lbws1,0);
%     
%     Rn = Rn1';
%     CIlw = quantile(XYbc(:,2:end),0.005)';
%     CIup = quantile(XYbc(:,2:end),0.995)';
%     CImd = mean(XYbc(:,2:end))';
% 
%     mCIup(i,1)   = mean(CIup);
%     mCIlw(i,1)   = mean(CIlw);
%     mCImd(i,1)   = mean(CImd);
%     
%     sYbc(i,1)    = mean(std(XYbc(:,2:end)))';
%     
%     fup(i,1)  = (mean(CIup) - 0) / mean(0.5*Rn); % Using Rn and lwbs to estimate upper 95% bound
%     flw(i,1)  = (0 - mean(CIlw)) / mean(0.5*Rn); % Using Rn and lwbs to estimate lower 95% bound
% 
%     toc
%     display(lbws1)
% end