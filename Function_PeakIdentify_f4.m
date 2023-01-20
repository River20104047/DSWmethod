function [Rn,XYbc,XYbn] = Function_PeakIdentify_f4(XY,lbws,ck)
    %% This is used for (rough) peak identification, especially for Raman and IR spectrum
    % This basic idea is to used estimated noize range/std, and then transfer
    % these range into smooth line
    
    % 2022/11/30 By Zijiang Yang
    % 2022/12/01 Output is XY format
    % 2022/12/01 Adding rough estimated noise range as output
    
    %% Input parameters
    % XY = Tsmp{:,1:20};   % XY = [Wavenumver Signals * n]; matrix format
    % ck  = 9;             % column of Y to check
    % lbws= 25;            % window size and stepsize for envelope, related to peak width
    
    es  = lbws;  % window size and stepsize for baseline correction = es
    fn  = 1.5;   % Skewness correction factor: since baseline corrected noise is not normally distributed, it is nagatively skewed, and the factor is ~1.5 to obtain 95% CI
    % Thus, fn = 1.5 around lbws = 10:1:100

    %% Processing
    % check X orders
    if XY(1,1) > XY(1:end)
       XY = flip(XY);
    end
        
    X   = XY(:,1);
    Y   = XY(:,2:end);
    
    % Envelope - upper bound
    Yeu   = msbackadj(X,Y,'WindowSize',es,'StepSize',es,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',1,'ShowPlot',0); 
    Yup   = Y - Yeu;
    
    % Envelope - lower bound
    Yel   = msbackadj(X,Y,'WindowSize',es,'StepSize',es,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',0,'ShowPlot',0); 
    Ylw   = Y - Yel;
    
    % Envelop - middle
    % Ymd   = 1/2*(Yup + Ylw);
    
    % Difference between lower bound and upper bound
    dYe   = Yup - Ylw;
    
    % Noise calculation
    % sigma = NaN(1,width(dYe));
    Rn    = NaN(1,width(dYe));
    F     = NaN(5001,width(dYe));
    Xi    = NaN(5001,width(dYe));
    
    for i = 1:1:width(dYe)
        dye   = dYe(:,i); % choose the i-th spectrum
        
        index9= dye < quantile(dye,0.99); 
    
        % in case 99% quantile = min
        if sum(index9) == 0
            index9(1) = 1;
        end
    
        dye = dye(index9); % remove top 1% for faster calculation
        
        [f,xi] = ksdensity(dye,(min(dye):((max(dye) - min(dye))/5000):max(dye))'); % Kernel pdf of dye based on 0-99%
        F(:,i)  = f;
        Xi(:,i) = xi;
        
        index  = f == max(f);
        
        % in case multiple elements in f = max(f)
        if sum(index) > 1
            Rn(i) = 0;
        else
        Rn(i)  = xi(index);    % range of noise
        % sigma(i)  = Rn(i)/3;      % std of noise
        end
    end

    % Baseline correction
    Ybc   = msbackadj(X,Y,'WindowSize',lbws,'StepSize',lbws,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',0,'ShowPlot',0); 
    BC    = Y - Ybc;          % estimated baseline
    BC2   = BC + 0.5 * Rn;    % unbiased estimate of baseline
    Ybc   = Y  - BC2;         % unbiased estimate of baseline corrected spectrum
    
    % Noise removal
    Ybn0  = Ybc - fn.* Rn;   % fn.* Rn = radiu of noise removal range, fn = 1.5 results ~99% noise
    index0 = Ybn0 > 0;
    Ybn   = Ybn0 .* index0;
    
    % Peak height calculation and signal-noise ratio
    % Hp    = max(Ybc) - min(Ybc);    % peak heigh = max peak hight after baseline correction
    % snr   = Hp ./ sigma;
    
    % Output
    XYbc  = [X,Ybc];
    XYbn  = [X,Ybn];

%% Check plots
    if ck > 0
    figure
%         
%         % Plot of original spectrum and estimated baseline
%         subplot(3,1,1)
%         plot(X,Y(:,ck),'lineWidth',1.2,'Color','#6699FF')
%         hold on
%         plot(X,BC2(:,ck),'lineWidth',1.2,'Color','#FF7C80')
%         title('Spectrum with estimated baseline')
%         xlim([min(X) max(X)])
%         set(gca, 'XDir','reverse')
%         
%         % Plot of dYe and estimated noise range (Rn)
%         subplot(3,1,2)
%         plot(X,Ybc(:,ck),'lineWidth',1.2,'Color','#00CC99')
%         hold on
%         plot([min(X) max(X)],[0 0],'lineWidth',1.2,'Color','#FF7C80')
%         plot([min(X) max(X)],[fn*Rn(ck) fn*Rn(ck)],'lineWidth',0.5,'Color','#FF7C80')
%         plot([min(X) max(X)],[-fn*Rn(ck) -fn*Rn(ck)],'lineWidth',0.5,'Color','#FF7C80')
%         title('Baseline corrected spectrum')
%         xlim([min(X) max(X)])
%         set(gca, 'XDir','reverse')
    
        % pdf of noise range (3 std of noise)
        % subplot(3,1,3)
        plot(Xi(:,ck),F(:,ck),'lineWidth',1.2,'Color','#CC99FF')

        title('PDF of noise range (d_n_s)')
        xlabel('Noise range')
        ylabel('Probability density')
        
    
        
        % Baseline corrected spectrum
    %     subplot(3,1,3)
    %     plot(X,Ybn(:,ck),'lineWidth',1.2,'Color','#6699FF')
    %     set(gca, 'XDir','reverse')
    %     title('Baseline corrected & noise removed spectrum')
    %     xlim([min(X) max(X)])
        
        % set(gcf, 'Position', [1            1       1/2*2194.3       1234.3])
    end
end

