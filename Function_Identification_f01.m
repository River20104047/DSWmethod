function [T_HQI,T_Rsmp,T_copolymer] = Function_Identification_f01(Tpre2smp,Tpre2std,hq,di)
% For version note see Identification control file

%% Preparation
% Alignment
xsmp    = Tpre2smp{:,1};
xstd    = Tpre2smp{:,1};
if xsmp ~= xstd
    msgbox('X values of standard and sample are not aligned',"Error","error");
end
clear xsmp xstd

% prepare spectra data
% x     = Tpre2smp{:,1}; No need for x -  correlation coeffecient is independent from x
ysmp    = Tpre2smp{:,2:end};
ystd    = Tpre2std{:,2:end};

ID.smp  = Tpre2smp.Properties.VariableNames(2:end);         % Spectrum ID
ID.std  = Tpre2std.Properties.VariableNames(2:end); 

%% correlation coeffecient
Rsmp    = zeros([width(ysmp),width(ystd)]);

rs01  = strings; rn01 = strings; r01 = zeros([width(ysmp),1]);  % r of string
rs02  = strings; rn02 = strings; r02 = zeros([width(ysmp),1]);
rs03  = strings; rn03 = strings; r03 = zeros([width(ysmp),1]);

for i = 1:1:width(ysmp)    
    % Zero-spectra detection - if spectrum is of pure noise, then processed signal may be zero straight line, and -1 was assigned to it (the least likely spectrum)
    if sum(ysmp(:,1)) == 0 
        ysmp(1,i) = -1;     
    end

    r  = corrcoef([ystd,ysmp(:,i)]);
    R = r(end,1:end-1);                  % used for further analysis
    Rsmp(i,:) = R;

    % Get R ( = HQI) of the first 3 greast type
    R_sorted = sort(R,'descend');

    r01(i,1) = R_sorted(1); % largest R
    r02(i,1) = R_sorted(2); % second largest R
    r03(i,1) = R_sorted(3); % third largest R

    rn01(i,1) = num2str(100*R_sorted(1),4); % r of number into string format
    rn02(i,1) = num2str(100*R_sorted(2),4);
    rn03(i,1) = num2str(100*R_sorted(3),4);

    index1   = find(R == R_sorted(1)); rs01(i,1) = ID.std{index1};
    index2   = find(R == R_sorted(2)); rs02(i,1) = ID.std{index2};
    index3   = find(R == R_sorted(3)); rs03(i,1) = ID.std{index3};
end

VarNames = {'Sample ID','Most likely','HQI1','Second likely','HQI2','Third likely','HQI3'};


%% co-polymer (up to two)
rcopolymer = rs01; % HQI of copolymer = HQI for the most likely type
for i = 1:1:width(ysmp)
    if r01(i) < hq
        rcopolymer(i) = 'unknown';                
    elseif r01(i) - r02(i) < di
        rcopolymer(i) = append(rs01(i),'/',rs02(i));
    end
end

%% Organize results
T_HQI   = table(ID.smp',rs01,rn01,rs02,rn02,rs03,rn03,'VariableNames',VarNames);%,'VariableNames',VarNames);
T_Rsmp      = array2table(Rsmp,'RowNames',ID.smp,'VariableNames',ID.std);
T_copolymer = table(ID.smp',rcopolymer,rn01,rn02,'VariableNames',{'Sample ID','Type','HQI-01','HQI-02'});

end