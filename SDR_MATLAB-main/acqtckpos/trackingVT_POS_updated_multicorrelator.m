function [TckResultVT_mltCorr, navSolutionsVT_mltCorr] = trackingVT_POS_updated_multicorrelator(file,signal,track,cmn,solu, Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT,navSolutionsCT)
%Purpose:
%   Vector tracking and positioning
%Inputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%	TckResultCT             - conventional tracking results
%	navSolutionsCT          - navigation solutions in conventional tracking mode
%Outputs:
%	TckResultVT         - vector tracking results
%	navSolutionsVT   	- navigation solutions in vector tracking mode
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
%
% Copyright (C) X X
% Written by X X

%%
% Spacing = -0.6:0.05:0.6;  % [-0.5,0,0.5];%
Spacing = 0.6:-0.05:-0.6;  % [-0.5,0,0.5];%
% Spacing = 0.7:-0.05:-0.7;  % [-0.5,0,0.5];%
datalength = track.msToProcessVT; 

sv  = Acquired.sv;
svlength    = length(Acquired.sv);
sv_clk      = zeros(1,32);
sv_clk_pos  = zeros(1,32);
eph_idx     = ones(1,svlength);

% Kalman Filter Parameter
num_state = 8;
error_state = zeros(num_state,1);
Dynamic_Model = diag([0,0,0,0,0,0,0,0]);
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model * track.pdi * signal.ms;

state_cov = 1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

process_noise(1:3,1:3) = diag(ones(1,3)*1e0);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;


% for experimental variance estimation
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 200/track.pdi;

% PVT initialization based on STL
estPos      = navSolutionsCT.usrPos(file.skiptimeVT/(solu.navSolPeriod),:);  
estVel      = navSolutionsCT.usrVel(file.skiptimeVT/(solu.navSolPeriod),:); 
clkBias     = navSolutionsCT.clkBias(file.skiptimeVT/(solu.navSolPeriod));   
clkDrift    = navSolutionsCT.clkDrift(file.skiptimeVT/(solu.navSolPeriod));
total_state = [estPos,estVel,clkBias,clkDrift]';

oldCarrError    = zeros(1,svlength);
carrError       = zeros(1,svlength);
codeError       = zeros(1,svlength);
oldCarrNco      = zeros(1,svlength);

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

% parameter for updating the iono and tropo correction.
corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(track.pdi*signal.ms);
counter_corr    = corrUpt-1 * ones(svlength,1);

% Find the start of VTL, in unit of samples  
sampleStart = zeros(1, svlength); % the first subframe 
for svindex = 1:svlength
    prn = sv(svindex);
    sampleStart(svindex) = ...
        TckResult_Eph(prn).absoluteSample(sbf.nav1(prn)+eph(prn).sfb(1)*20); 
end
sampleStartMeaCT = max(sampleStart) + 1; % start of the measurement of STL, in unit of samples
measSampleStepCT = fix(signal.Fs * solu.navSolPeriod/1000)*1;%file.dataType; % measurement step of STL, in unit of samples
measStartCT = sampleStartMeaCT;% + measSampleStepCT; % first measurement epoch of STL, in unit of samples 

sampleStartTckVT = measStartCT + file.skiptimeVT/(solu.navSolPeriod); % First sample for VT, same to all channels % 22 Jun 2021
msStartTckVT = sbf.nav1(prn)+eph(prn).sfb(1)*20 + file.skiptimeVT/solu.navSolPeriod;


% Tracking initialization based on STL
for svindex = 1:svlength
    prn = sv(svindex);    
    codetemp                = generateCAcode(Acquired.sv(svindex));
    Code(svindex,:)         = [codetemp(end) codetemp(end-1) repmat(codetemp,1,track.pdi) codetemp(1) codetemp(2)];
    codeFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).codeFreq(msStartTckVT); %fix(sampleStartTckVT/2/signal.Fs)
    codePhaseStep(svindex)  = codeFreq(svindex)/signal.Fs;
    remChip(svindex)        = TckResultCT(Acquired.sv(svindex)).remChip(msStartTckVT);
    carrFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT);
    oldCarrFreq(svindex)    = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT);
    remCarrPhase(svindex)   = TckResultCT(Acquired.sv(svindex)).remCarrPhase(msStartTckVT);
    file_ptr(svindex)       = TckResultCT(Acquired.sv(svindex)).absoluteSample(msStartTckVT);
    carrFreqBasis(svindex)  = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT);
    codedelay_tck(svindex)  = TckResultCT(Acquired.sv(svindex)).codedelay(msStartTckVT);
    oldCarrError(svindex)   = TckResultCT(Acquired.sv(svindex)).carrError(msStartTckVT);
    oldCarrNco(svindex)     = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT) - carrFreqBasis(svindex);
        
    svxyzr_tck_last(svindex,:) = zeros(1,3);    
    
    correction(svindex) = 0; % NLOS correction
%     transmitTimeVT(svindex) = eph(prn).TOW(1) + (sampleStartTckVT - sampleStart(svindex))/2/signal.Fs; % transmit time of the beginning of VTL
%     transmitTimeVT(svindex) = eph(prn).TOW(1) + file.skiptimeVT/solu.navSolPeriod/1000;
    transmitTimeVT(svindex) = navSolutionsCT.timeTransmit(1,svindex);
end

% PLL filter coefficientS
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

% 
amplitude = 0;
navi_data = 0;
navi_dataL035 = 0; 
deltaPr = zeros(1,length(Acquired.sv));
prRate = zeros(1,length(Acquired.sv));  
clkBias_last = clkBias;
estPos_last = estPos;


% 
localTime = min(transmitTimeVT);

% flag for position. When the file pointers in all tracking channels exceed
% the current measurement point (in samples), we do positioning
flg_pos = zeros(1,svlength);

%
posIndex = 0;

%% Start processing
h = waitbar(0,['Vector Tracking, Length: ',num2str(datalength),' ms,', '  Please wait...']);
tic
for msIndex = 1:datalength/track.pdi % 
    %waitbar(msIndex/(datalength/track.pdi),h)
    for svindex = 1:length(Acquired.sv)
        prn = sv(svindex);
        numSample(svindex) = ceil((signal.codelength*track.pdi-remChip(svindex))/(codeFreq(svindex)/signal.Fs));  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(file.fid, file_ptr(svindex)*1,'bof');  % file_ptr has already considered the data type of 'int16'.
        if file.dataPrecision == 2
            [rawsignal, ~] = fread(file.fid,numSample(svindex)*file.dataType,'int16') ;
            rawsignal = rawsignal';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            [rawsignal, ~]  = fread(file.fid,numSample(svindex)*file.dataType,'int8');
            rawsignal = rawsignal';
            if file.dataType == 2
                rawsignal = rawsignal(1:2:length(rawsignal)) + 1i*rawsignal(2:2:length(rawsignal));% For NSL STEREO LBand only
            end
        end
                
        % transmitTime of the last sample of this block
        transmitTimeVT(svindex) =  transmitTimeVT(svindex) + numSample(svindex)/signal.Fs;   
        tot_est_tck(svindex) = transmitTimeVT(svindex);
              
        % SV position and Velocity at the transmitTime
        [svxyz_tck(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel(prn)] = ...
            svPosVel(prn,eph, tot_est_tck(svindex),eph_idx(svindex));
        
        %% Iono, trop correction, follow the book of Paul: pp.268 and eq(7.34) (svxyz - estPos).
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt
            svenu           = xyz2enu(svxyz_tck(svindex,:), estPos);
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            %             az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az_rad(svindex) = atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex) * 180/pi;
            el(svindex)     = el_rad(svindex) * 180/pi;
            
            temp = xyz2llh(estPos);
            user_ll	= [temp(1:2).*180/pi temp(3)];
            ionodel(svindex)        = ionocorr(tot_est_tck(svindex), svxyz_tck(svindex,:), cnslxyz); %estPos(1:3)
            tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex))); %
            
            counter_corr(svindex)   = 0;
        end
        
        
        %% Predict code freq. based on the navigation solution
        r = sqrt(sum((svxyz_tck(svindex,:) - estPos).^2));
        % Apply corrections
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex);%
        % earth rotation correction
        svxyzr_tck(svindex,:) = erotcorr(svxyz_tck(svindex,:),predictedPr_tck(svindex));
        % Corrected range
        r = sqrt(sum((svxyzr_tck(svindex,:) - estPos).^2));
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex) ;%
        
        % Predicted code freq.
        if msIndex == 1
            codeFreq(svindex) = TckResultCT(Acquired.sv(svindex)).codeFreq(msStartTckVT);%fix(sampleStartTckVT/2/signal.Fs)
        else
            deltaPr(svindex) = (predictedPr_tck(svindex) - predictedPr_last(svindex))/(track.pdi*signal.ms);
            codeFreq(svindex) = signal.codeFreqBasis*(1 - deltaPr(svindex)/cmn.cSpeed); % 
        end 
        predictedPr_last(svindex) = predictedPr_tck(svindex); 
           
        % new code phase step
        codePhaseStep(svindex) = codeFreq(svindex)/signal.Fs;
        
        %% spacing = -0.6:0.05:0.6
%         t_CodeEarly       = (0 + Spacing(5) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(5) + remChip(svindex));
%         t_CodePrompt      = (0 + Spacing(15) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(15) + remChip(svindex));
%         t_CodeLate        = (0 + Spacing(25) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(25) + remChip(svindex));
%        
%         indx = 1;%%%
%         CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + indx));
%         CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + indx));
%         CodeLate       = Code(svindex,(ceil(t_CodeLate) + indx));
 
        t_CodeEarly_060   = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex)); 
        t_CodeEarly_055   = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));        
        t_CodeEarly       = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex)); 
        t_CodeEarly_045   = (0 + Spacing(4) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(4) + remChip(svindex));        
        t_CodeEarly_040   = (0 + Spacing(5) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex)+ Spacing(5) + remChip(svindex));       
        t_CodeEarly_035   = (0 + Spacing(6) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(6) + remChip(svindex));       
        t_CodeEarly_030   = (0 + Spacing(7) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(7) + remChip(svindex));       
        t_CodeEarly_025   = (0 + Spacing(8) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(8) + remChip(svindex));       
        t_CodeEarly_020   = (0 + Spacing(9) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(9) + remChip(svindex));       
        t_CodeEarly_015   = (0 + Spacing(10) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(10) + remChip(svindex));       
        t_CodeEarly_010   = (0 + Spacing(11) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(11) + remChip(svindex));       
        t_CodeEarly_005   = (0 + Spacing(12) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(12) + remChip(svindex));
        t_CodePrompt      = (0 + Spacing(13) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(13) + remChip(svindex));
        t_CodeLate005     = (0 + Spacing(14) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(14) + remChip(svindex));
        t_CodeLate010     = (0 + Spacing(15) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(15) + remChip(svindex));
        t_CodeLate015     = (0 + Spacing(16) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(16) + remChip(svindex));
        t_CodeLate020     = (0 + Spacing(17) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(17) + remChip(svindex));
        t_CodeLate025     = (0 + Spacing(18) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(18) + remChip(svindex));
        t_CodeLate030     = (0 + Spacing(19) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(19) + remChip(svindex));
        t_CodeLate035     = (0 + Spacing(20) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(20) + remChip(svindex));
        t_CodeLate040     = (0 + Spacing(21) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(21) + remChip(svindex));
        t_CodeLate045     = (0 + Spacing(22) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(22) + remChip(svindex));
        t_CodeLate        = (0 + Spacing(23) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(23) + remChip(svindex));
        t_CodeLate055     = (0 + Spacing(24) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(24) + remChip(svindex));
        t_CodeLate060     = (0 + Spacing(25) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(25) + remChip(svindex));

        index = 2;
        CodeEarly_060 = Code(svindex,(ceil(t_CodeEarly_060) + index));
        CodeEarly_055 = Code(svindex,(ceil(t_CodeEarly_055) + index));
        CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + index));
        CodeEarly_045 = Code(svindex,(ceil(t_CodeEarly_045) + index));
        CodeEarly_040 = Code(svindex,(ceil(t_CodeEarly_040) + index));
        CodeEarly_035 = Code(svindex,(ceil(t_CodeEarly_035) + index));
        CodeEarly_030 = Code(svindex,(ceil(t_CodeEarly_030) + index));
        CodeEarly_025 = Code(svindex,(ceil(t_CodeEarly_025) + index));
        CodeEarly_020 = Code(svindex,(ceil(t_CodeEarly_020) + index));
        CodeEarly_015 = Code(svindex,(ceil(t_CodeEarly_015) + index));
        CodeEarly_010 = Code(svindex,(ceil(t_CodeEarly_010) + index));
        CodeEarly_005 = Code(svindex,(ceil(t_CodeEarly_005) + index));
        CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + index));
        CodeLate005    = Code(svindex,(ceil(t_CodeLate005) + index));
        CodeLate010    = Code(svindex,(ceil(t_CodeLate010) + index));
        CodeLate015    = Code(svindex,(ceil(t_CodeLate015) + index));
        CodeLate020    = Code(svindex,(ceil(t_CodeLate020) + index));
        CodeLate025    = Code(svindex,(ceil(t_CodeLate025) + index));
        CodeLate030    = Code(svindex,(ceil(t_CodeLate030) + index));
        CodeLate035    = Code(svindex,(ceil(t_CodeLate035) + index));
        CodeLate040    = Code(svindex,(ceil(t_CodeLate040) + index));
        CodeLate045    = Code(svindex,(ceil(t_CodeLate045) + index));
        CodeLate       = Code(svindex,(ceil(t_CodeLate) + index));
        CodeLate055    = Code(svindex,(ceil(t_CodeLate055) + index));
        CodeLate060    = Code(svindex,(ceil(t_CodeLate060) + index));
      
        
        %% mixed with local carrier replica
        CarrTime = (0: numSample(svindex)) ./ signal.Fs;
        Wave = (2*pi*((carrFreq(svindex)) .* CarrTime)) + remCarrPhase(svindex);
        %         remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1),2*pi);
        
        carrsig = exp(1i.* Wave(1:numSample(svindex)));
        InphaseSignal    = imag(rawsignal .* carrsig); %
        QuadratureSignal = real(rawsignal .* carrsig); 
        
        %%
%         E_i  = sum(CodeEarly    .*InphaseSignal);
%         E_q = sum(CodeEarly    .*QuadratureSignal);
%         P_i  = sum(CodePrompt   .*InphaseSignal);
%         P_q = sum(CodePrompt   .*QuadratureSignal);
%         L_i  = sum(CodeLate     .*InphaseSignal);
%         L_q = sum(CodeLate     .*QuadratureSignal); 

        E_i_060  = sum(CodeEarly_060    .*InphaseSignal);
        E_q_060  = sum(CodeEarly_060    .*QuadratureSignal);
        E_i_055  = sum(CodeEarly_055    .*InphaseSignal);
        E_q_055  = sum(CodeEarly_055    .*QuadratureSignal);
        E_i  = sum(CodeEarly    .*InphaseSignal);  
        E_q = sum(CodeEarly    .*QuadratureSignal);
        E_i_045  = sum(CodeEarly_045    .*InphaseSignal);
        E_q_045  = sum(CodeEarly_045    .*QuadratureSignal);
        E_i_040  = sum(CodeEarly_040    .*InphaseSignal);
        E_q_040  = sum(CodeEarly_040    .*QuadratureSignal);
        E_i_035  = sum(CodeEarly_035    .*InphaseSignal);
        E_q_035  = sum(CodeEarly_035    .*QuadratureSignal);
        E_i_030  = sum(CodeEarly_030    .*InphaseSignal);
        E_q_030  = sum(CodeEarly_030    .*QuadratureSignal);
        E_i_025  = sum(CodeEarly_025    .*InphaseSignal);
        E_q_025  = sum(CodeEarly_025    .*QuadratureSignal);
        E_i_020  = sum(CodeEarly_020    .*InphaseSignal);
        E_q_020  = sum(CodeEarly_020    .*QuadratureSignal);
        E_i_015  = sum(CodeEarly_015    .*InphaseSignal);
        E_q_015  = sum(CodeEarly_015    .*QuadratureSignal);
        E_i_010  = sum(CodeEarly_010    .*InphaseSignal);
        E_q_010  = sum(CodeEarly_010    .*QuadratureSignal);
        E_i_005  = sum(CodeEarly_005    .*InphaseSignal);
        E_q_005  = sum(CodeEarly_005    .*QuadratureSignal);
        P_i  = sum(CodePrompt   .*InphaseSignal);  
        P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i005  = sum(CodeLate005     .*InphaseSignal); 
        L_q005  = sum(CodeLate005     .*QuadratureSignal); 
        L_i010  = sum(CodeLate010     .*InphaseSignal); 
        L_q010  = sum(CodeLate010     .*QuadratureSignal);
        L_i015  = sum(CodeLate015     .*InphaseSignal); 
        L_q015  = sum(CodeLate015     .*QuadratureSignal);
        L_i020  = sum(CodeLate020     .*InphaseSignal);  
        L_q020  = sum(CodeLate020     .*QuadratureSignal);
        L_i025  = sum(CodeLate025     .*InphaseSignal);  
        L_q025  = sum(CodeLate025     .*QuadratureSignal);
        L_i030  = sum(CodeLate030     .*InphaseSignal);  
        L_q030  = sum(CodeLate030     .*QuadratureSignal);
        L_i035  = sum(CodeLate035     .*InphaseSignal);  
        L_q035  = sum(CodeLate035     .*QuadratureSignal);
        L_i040  = sum(CodeLate040     .*InphaseSignal);  
        L_q040  = sum(CodeLate040     .*QuadratureSignal);
        L_i045  = sum(CodeLate045     .*InphaseSignal);  
        L_q045  = sum(CodeLate045     .*QuadratureSignal);
        L_i  = sum(CodeLate     .*InphaseSignal); 
        L_q = sum(CodeLate     .*QuadratureSignal); 
        L_i055  = sum(CodeLate055     .*InphaseSignal);  
        L_q055  = sum(CodeLate055     .*QuadratureSignal);
        L_i060  = sum(CodeLate060     .*InphaseSignal); 
        L_q060  = sum(CodeLate060     .*QuadratureSignal); 
   
        
        % remaining code and carrier phase
        remChip(svindex) = (t_CodePrompt(numSample(svindex)) + codePhaseStep(svindex)) - 1023*track.pdi;  
        remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1),2*pi);
        
        % some parameters at last epcoh
        clkBias_last = clkBias;
        estPos_last = estPos;
        svxyzr_tck_last(svindex,:) = svxyzr_tck(svindex,:); 
        
        %% CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_VT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*signal.ms*track.pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        %% PLL discriminator
        carrError(svindex)      = atan(P_q/P_i)/(2.0 * pi);
        carrNco(svindex)        = oldCarrNco(svindex) + (tau2carr/tau1carr) * (carrError(svindex) - oldCarrError(svindex)) + carrError(svindex)*(track.pdi*1e-3/tau1carr) ;
        oldCarrNco(svindex)     = carrNco(svindex);
        oldCarrError(svindex)   = carrError(svindex);
        carrFreq(svindex)       = carrFreqBasis(svindex) + carrNco(svindex);%
        oldCarrFreq(svindex)    = carrFreq(svindex);
        
        %% DLL discriminator
        E	= sqrt(E_i^2+E_q^2);
        L	= sqrt(L_i^2+L_q^2);
        codeError(svindex)	= -0.5*(E-L)/(E+L); 
        
        % pseudorange measurements; pseudorange error correction
        Z(msIndex,svindex) = codeError(svindex)*cmn.cSpeed/codeFreq(svindex);
        
        %% Result Record
%         TckResultVT(Acquired.sv(svindex)).E_i(msIndex)            = E_i;
%         TckResultVT(Acquired.sv(svindex)).E_q(msIndex)            = E_q;
%         TckResultVT(prn).P_i(msIndex)                             = P_i;
%         TckResultVT(prn).P_q(msIndex)                             = P_q;
%         TckResultVT(Acquired.sv(svindex)).L_i(msIndex)            = L_i;
%         TckResultVT(Acquired.sv(svindex)).L_q(msIndex)            = L_q;       

        TckResultVT(Acquired.sv(svindex)).E_i_060(msIndex)         = E_i_060;
        TckResultVT(Acquired.sv(svindex)).E_i_055(msIndex)         = E_i_055;
        TckResultVT(Acquired.sv(svindex)).E_i(msIndex)            = E_i;
        TckResultVT(Acquired.sv(svindex)).E_i_045(msIndex)         = E_i_045;
        TckResultVT(Acquired.sv(svindex)).E_i_040(msIndex)         = E_i_040;
        TckResultVT(Acquired.sv(svindex)).E_i_035(msIndex)         = E_i_035;
        TckResultVT(Acquired.sv(svindex)).E_i_030(msIndex)         = E_i_030;
        TckResultVT(Acquired.sv(svindex)).E_i_025(msIndex)         = E_i_025;
        TckResultVT(Acquired.sv(svindex)).E_i_020(msIndex)         = E_i_020;
        TckResultVT(Acquired.sv(svindex)).E_i_015(msIndex)         = E_i_015;
        TckResultVT(Acquired.sv(svindex)).E_i_010(msIndex)         = E_i_010;
        TckResultVT(Acquired.sv(svindex)).E_i_005(msIndex)         = E_i_005;
        TckResultVT(Acquired.sv(svindex)).E_q_060(msIndex)         = E_q_060;
        TckResultVT(Acquired.sv(svindex)).E_q_055(msIndex)         = E_q_055;
        TckResultVT(Acquired.sv(svindex)).E_q(msIndex)            = E_q;
        TckResultVT(Acquired.sv(svindex)).E_q_045(msIndex)         = E_q_045;
        TckResultVT(Acquired.sv(svindex)).E_q_040(msIndex)         = E_q_040;
        TckResultVT(Acquired.sv(svindex)).E_q_035(msIndex)         = E_q_035;
        TckResultVT(Acquired.sv(svindex)).E_q_030(msIndex)         = E_q_030;
        TckResultVT(Acquired.sv(svindex)).E_q_025(msIndex)         = E_q_025;
        TckResultVT(Acquired.sv(svindex)).E_q_020(msIndex)         = E_q_020;
        TckResultVT(Acquired.sv(svindex)).E_q_015(msIndex)         = E_q_015;
        TckResultVT(Acquired.sv(svindex)).E_q_010(msIndex)         = E_q_010;
        TckResultVT(Acquired.sv(svindex)).E_q_005(msIndex)         = E_q_005;
        TckResultVT(prn).P_i(msIndex)                             = P_i;
        TckResultVT(prn).P_q(msIndex)                             = P_q; 
        TckResultVT(Acquired.sv(svindex)).L_i005(msIndex)            = L_i005;
        TckResultVT(Acquired.sv(svindex)).L_i010(msIndex)            = L_i010;
        TckResultVT(Acquired.sv(svindex)).L_i015(msIndex)            = L_i015;
        TckResultVT(Acquired.sv(svindex)).L_i020(msIndex)            = L_i020;
        TckResultVT(Acquired.sv(svindex)).L_i025(msIndex)            = L_i025;
        TckResultVT(Acquired.sv(svindex)).L_i030(msIndex)            = L_i030;
        TckResultVT(Acquired.sv(svindex)).L_i035(msIndex)            = L_i035;
        TckResultVT(Acquired.sv(svindex)).L_i040(msIndex)            = L_i040;
        TckResultVT(Acquired.sv(svindex)).L_i045(msIndex)            = L_i045;
        TckResultVT(Acquired.sv(svindex)).L_i(msIndex)            = L_i;
        TckResultVT(Acquired.sv(svindex)).L_i055(msIndex)            = L_i055;
        TckResultVT(Acquired.sv(svindex)).L_i060(msIndex)            = L_i060;
        TckResultVT(Acquired.sv(svindex)).L_q005(msIndex)            = L_q005;
        TckResultVT(Acquired.sv(svindex)).L_q010(msIndex)            = L_q010;
        TckResultVT(Acquired.sv(svindex)).L_q015(msIndex)            = L_q015;
        TckResultVT(Acquired.sv(svindex)).L_q020(msIndex)            = L_q020;
        TckResultVT(Acquired.sv(svindex)).L_q025(msIndex)            = L_q025;
        TckResultVT(Acquired.sv(svindex)).L_q030(msIndex)            = L_q030;
        TckResultVT(Acquired.sv(svindex)).L_q035(msIndex)            = L_q035;
        TckResultVT(Acquired.sv(svindex)).L_q040(msIndex)            = L_q040;
        TckResultVT(Acquired.sv(svindex)).L_q045(msIndex)            = L_q045;
        TckResultVT(Acquired.sv(svindex)).L_q(msIndex)            = L_q;
        TckResultVT(Acquired.sv(svindex)).L_q055(msIndex)            = L_q055; 
        TckResultVT(Acquired.sv(svindex)).L_q060(msIndex)            = L_q060;  
 
        
        
        
        
        TckResultVT(Acquired.sv(svindex)).amplitude(msIndex)            = amplitude;
        TckResultVT(Acquired.sv(svindex)).navi_data(msIndex)            = sum(navi_data);
        TckResultVT(Acquired.sv(svindex)).navi_dataL035(msIndex)        = sum(navi_dataL035);
        
        
        %%
        TckResultVT(prn).carrError(msIndex)          = carrError(svindex);
        TckResultVT(prn).codeError(msIndex)          = codeError(svindex);
        TckResultVT(prn).remChip(msIndex)            = remChip(svindex);
        TckResultVT(prn).remCarrPhase(msIndex)       = remCarrPhase(svindex);
        TckResultVT(prn).codeFreq(msIndex)           = codeFreq(svindex);
        TckResultVT(prn).carrFreq(msIndex)           = carrFreq(svindex);
        TckResultVT(prn).carrNco(msIndex)           = carrNco(svindex);
        TckResultVT(prn).absoluteSample(msIndex)     = ftell(file.fid);
        file_ptr(svindex)                            = TckResultVT(prn).absoluteSample(msIndex);
        TckResultVT(prn).sv_vel(msIndex,:)           = sv_vel(svindex,:);
        TckResultVT(prn).codedelay(msIndex)          = mod(TckResultVT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),signal.Fs*signal.ms);
        
        codedelay_tck(svindex)                       = TckResultVT(prn).codedelay(msIndex);
        
        TckResultVT(prn).deltaPr(msIndex)           = deltaPr(svindex);
        TckResultVT(prn).prRate(msIndex)           = prRate(svindex);        
    end % end for svindex in Tracking
     
    
    %% Navigation solution solving  
    numSample_min = min(numSample)-1;
    
    
    for svindex = 1:svlength
        prn = sv(svindex);
        
        tot_est_pos(svindex) = tot_est_tck(svindex) - (numSample(svindex)-numSample_min)/signal.Fs;           
        localTime = min(tot_est_pos);
    
        [svxyz_pos(svindex,:), sv_vel_pos(svindex,:), sv_clk_pos(prn), sv_clk_vel(prn), grpdel(prn)] = ...
            svPosVel(prn,eph, tot_est_pos(svindex), eph_idx(svindex));
        
        r = sqrt(sum((svxyz_pos(svindex,:) - estPos).^2));
        predictedPr_pos(svindex) = r  + clkBias + sv_clk_pos(prn) - grpdel(prn)*cmn.cSpeed ...
            - tropodel_unb3(svindex) - ionodel(svindex);    %
        svxyzr_pos(svindex,:) = erotcorr(svxyz_pos(svindex,:),(predictedPr_pos(svindex)));%
        r = sqrt(sum((svxyzr_pos(svindex,:) - estPos).^2));
        a_pos(svindex,:) = (svxyzr_pos(svindex,:)-estPos) / r;
        H_pos(svindex,:) = [-a_pos(svindex,:) 0 0 0 1 0];
        H_pos(svindex+svlength,:) = [0 0 0 -a_pos(svindex,:) 0 1];
        
        % Find pseudorange rate error
        prr_measured(svindex)	= (carrFreq(svindex) - signal.IF)*cmn.cSpeed/signal.Fc; %% 22 Jun 2021, for Labsat3w
        prr_predicted(svindex)	= (estVel - sv_vel_pos(svindex,:))*a_pos(svindex,:)'; %
        Z(msIndex,svlength+svindex) = prr_predicted(svindex) - prr_measured(svindex) - clkDrift + sv_clk_vel(prn);
    end
    
    newZ = Z(msIndex,:); % complete measurements
    
    % Kalman filter
    error_state = zeros(num_state,1);
    error_state = Transistion_Matrix * error_state;
    
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;
    kalman_gain = state_cov * transpose(H_pos) * inv(H_pos * state_cov * transpose(H_pos) + mesurement_noise);
    
    counterUptR = counterUptR + 1;  % counter for update measurement noise variance, R
    recordR(counterUptR,:) = ((newZ' - H_pos * error_state));
    
    error_state = error_state + kalman_gain * (newZ' - H_pos * error_state);
    state_cov = (eye(num_state) - kalman_gain * H_pos) * state_cov;
        
    total_state =  total_state + error_state;
    estPos = total_state(1:3)';
    estVel = total_state(4:6)';  %
    clkBias = total_state(7);
    clkDrift = total_state(8);
    
    %% record results
    llh     = xyz2llh(cnslxyz);  %
    L_b     = llh(1);
    lamda_b = llh(2);
    C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
        -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
        -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];
    usrenuvel(msIndex,:) = C_e_n * estVel';
    
    usrenu(msIndex,:) = xyz2enu(estPos,cnslxyz);
    usrllh(msIndex,:) = xyz2llh(estPos);
    usrllh(msIndex,1:2)	= usrllh(msIndex,1:2)*180/pi;
    navSolutionsVT.localTime(msIndex,:)        = localTime;
    navSolutionsVT.usrPos(msIndex,:)         = estPos;
    navSolutionsVT.usrVel(msIndex,:)         = estVel;
    navSolutionsVT.usrPosENU(msIndex,:)      = usrenu(msIndex,:);
    navSolutionsVT.usrVelENU(msIndex,:)      = usrenuvel(msIndex,:);
    navSolutionsVT.usrPosLLH(msIndex,:)   	 = usrllh(msIndex,:);
    navSolutionsVT.clkDrift(msIndex,:)  	 = clkDrift;
    navSolutionsVT.clkBias(msIndex,:)        = clkBias;
    navSolutionsVT.satePos(msIndex,:)   	 = svxyzr_pos(svindex,:);
    navSolutionsVT.sateVel(msIndex,:)   	 = sv_vel_pos(svindex,:);
    navSolutionsVT.state(msIndex,:)          = error_state;
    navSolutionsVT.svxyz_pos(:,:,msIndex)    = svxyz_pos;
    navSolutionsVT.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsVT.state_cov(msIndex,:)      = diag(state_cov);
    navSolutionsVT.meas_inno(msIndex,:)      = ((newZ' - H_pos * error_state));
    navSolutionsVT.newZ(msIndex,:)           = newZ;
    navSolutionsVT.predicted_z(msIndex,:)    = H_pos * error_state;
    navSolutionsVT.satEA(msIndex,:)      = el;
    navSolutionsVT.satAZ(msIndex,:)      = az;
    %     navSolutionsVT.DOP(msIndex,:)      = [HDOP, TDOP]; 
    
    %% predict postion and clkBias at next epoch
    total_state  = Transistion_Matrix * total_state;
    estPos = total_state(1:3)';
    clkBias = total_state(7)'; 
    
    %% update Q and R by measurement variance
    if flag_corrCovEst2 == 1 && counterUptR == thresUptR
        tmpR = diag(1/counterUptR*(sum(recordR.^2)));
        mesurement_noise(1:svlength,1:svlength) = tmpR(1:svlength,1:svlength)*10; %
        mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = ...
            tmpR(svlength+1:2*svlength,svlength+1:2*svlength)*1;
        
        for idx = 1 : svlength
            if mesurement_noise(idx,idx) >= 12000
                mesurement_noise(idx,idx) = 12000;
            elseif mesurement_noise(idx,idx) <= 0.01
                mesurement_noise(idx,idx) = 0.01;
            end
            if mesurement_noise(idx+svlength,idx+svlength) >= 400
                mesurement_noise(idx+svlength,idx+svlength) = 400;
            elseif mesurement_noise(idx+svlength,idx+svlength) <= 0.01
                mesurement_noise(idx+svlength,idx+svlength) = 0.01;
            end
            
        end
        counterUptR = 0;
        counter_r = counter_r + 1;
        navSolutionsVT.R(counter_r,:) = diag(mesurement_noise);
    end
    
    navSolutionsVT.record_correction(msIndex,:)      = correction;
     
    %% 
    if mod(msIndex, 100) == 0
        disp([msIndex, localTime, usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3), clkBias, clkDrift]); 
    end
     
end % end for msIndex
close(h);
toc

TckResultVT_mltCorr = TckResultVT;
navSolutionsVT_mltCorr = navSolutionsVT;

%% Save results
save(['navSolVT_mltCorr_',file.fileName,'_updated'], 'navSolutionsVT_mltCorr','eph');
save(['tckRstVT_mltCorr_',file.fileName,'_updated'], 'TckResultVT_mltCorr','CN0_VT');


