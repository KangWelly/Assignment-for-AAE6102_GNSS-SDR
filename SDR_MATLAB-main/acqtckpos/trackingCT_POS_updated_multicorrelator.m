function [TckResultCT_mltCorr, navSolutionsCT_mltCorr] = trackingCT_POS_updated_multicorrelator(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph)
%Purpose:                                                                                                                         
%   Scalar tracking and positioning. Positioning and tracking are
%   implemented at the same time.
%   Conventional tracking and positioning using EKF and WLS
%Inputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%   TckResult_Eph - Tracking results that used for decoding eph., which
%   also contains information like indexes for first nav. bit transition, subframe,
%   absolute sample index in the IF file for each ms, etc
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%
%Outputs:
%	TckResultCT         - conventional tracking results
%	navSolutionsCT   	- navigation solutions in conventional tracking mode
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
%
% Written by B. XU and L. T. HSU


%%
sv_clk              = zeros(1,32);
clkBias_kf         	=  0;
usr_clk_wls         = 0;
clkDrift            = 0;
oldclkDrift         = 0;
estusr              = zeros(1,3);
estusr_wls          =cnslxyz;% zeros(1,3);
estusr_kf           = cnslxyz;%
estVel              = zeros(1,3);
oldestVel          	= estVel;
num_state           = 8;

Spacing = 0.6:-0.05:-0.6;

sv      = Acquired.sv ;
f0      = signal.codeFreqBasis;
fs      = signal.Fs ;
pdi     = track.pdi;
t       = signal.ms;
svlength    = length(sv);
datalength  = track.msPosCT;

% Kalman Filter Parameter
num_state   = 8;

% error state vector
error_state = zeros(num_state,1);
total_state = [cnslxyz(1:3),zeros(1,5)]';%zeros(num_state,1);

% system transition matrix
Dynamic_Model = diag(zeros(1,num_state));
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model*pdi*t;

% error covariance matrix
state_cov = 1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

% process (System) noise noise covariance matrix
process_noise(1:3,1:3) = diag(ones(1,3)*2e-1);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;

% measurement noise covariance matrix
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;

% parameters for measurement noise variance update
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 20/pdi;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);
%

% initialize tracking parameters using acquisition results
for svindex = 1:length(sv)
    prn                     = sv(svindex);
    codetemp                = generateCAcode(prn);
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)  codetemp(2)];
    AcqCodeDelay(svindex)   = Acquired.codedelay(svindex);
    %     file_ptr(svindex)       = signal.Sample - AcqCodeDelay(svindex) -1 + file.skip *fs*t + 0;
    %     file_ptr(svindex)       = signal.Sample - AcqCodeDelay(svindex) +1 + file.skip *fs*t + 0;   % 29/04/2020
    
    % 29/04/29, for each channel, move the file pointer to the first
    % subframe (not subframe 1), to ensure 20 ms T_coh and do positioning
    file_ptr(svindex)       = (signal.Sample - AcqCodeDelay(svindex) +1 ...
        + file.skip *fs*t ... 
        )*file.dataPrecision*file.dataType;
    
    carrFreq(svindex)       = Acquired.fineFreq(svindex);
    AcqFreq(svindex)        = Acquired.fineFreq(svindex);
    
    oldcodedelay_pos(svindex) = 0;
    oldabsoluteSample_pos(svindex) = 0;
end

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

eph_idx     = ones(1,svlength);

corrUpdateSec   = 0.01;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1);

% Tracking parameters
carrNco      = zeros(1,svlength);
oldCarrNco  = zeros(1,svlength);
oldCarrError       = zeros(1,svlength);
codeNco         = zeros(1,svlength);
code_outputLast     = zeros(1,svlength);
DLLdiscriLast       = zeros(1,svlength);
remChip             = zeros(1,svlength);
codeFreq            = ones(1,svlength)*f0;
remCarrPhase        = zeros(1,svlength);
carrError           = zeros(1,svlength);
codeError           = zeros(1,svlength);
delayValue          = zeros(svlength,datalength/pdi);

%%
localTime = inf;

% Find start and end of measurement point locations in IF signal stream with available
% measurements
sampleStart = zeros(1, svlength);
sampleEnd = inf(1, svlength);

for channelNr = 1:svlength
    prn = sv(channelNr);
    sampleStart(channelNr) = TckResult_Eph(prn).absoluteSample(sbf.nav1(prn)+eph(prn).sfb(1)*20); % first subframe, in unit of ms
    
    sampleEnd(channelNr) = TckResult_Eph(prn).absoluteSample(end);
end
sampleStartMea = max(sampleStart) + 1; % Now, in  unit of sample
sampleEndMea = min(sampleEnd) - 1;

%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(signal.Fs * solu.navSolPeriod/1000)*file.dataType;%file.dataType;

% flag for position. When the file pointers in all tracking channels exceed
% the current measurement point (in samples), we do positioning
flg_pos = zeros(1,svlength);

% Index for positioning
posIndex = 0;

%%
h = waitbar(0,['Conventional Tracking and Positioning, Length: ',num2str(datalength),' ms,', '  Please wait...']);
tic

%%
for msIndex = 1: datalength/pdi % Note that for pdi > 1ms, the index is still denoted as msIndex. 30/04/2020, BING XU
    %waitbar(msIndex/(datalength/pdi),h)
    for svindex = 1 :svlength
        prn = sv(svindex);
        
        % read raw data file
        codePhaseStep(svindex) = codeFreq(svindex)/signal.Fs;
        numSample = ceil((signal.codelength*pdi-remChip(svindex))/codePhaseStep(svindex));  %  /(1*2)
        
        delayValue(svindex,msIndex) = numSample - signal.Sample*pdi;
        
        fseek(file.fid, file_ptr(svindex),'bof');
        
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal = fread(file.fid,numSample*file.dataType,'int8')'; %
            if file.dataType == 2
                rawsignal = rawsignal(1:2:length(rawsignal)) + 1i*rawsignal(2:2:length(rawsignal)); 
            end
        end
        
        file_ptr(svindex)   = file_ptr(svindex) + numSample*file.dataType;  %%%%%%   
        
        %% spacing = -0.6:0.05:0.6
        % Updated by LZD on 2021.10.07
%         t_CodeEarly       = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex));
%         t_CodePrompt      = (0 + Spacing(13) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(13) + remChip(svindex));
%         t_CodeLate        = (0 + Spacing(23) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(23) + remChip(svindex));
%        
%         indx = 1;
%         CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + indx));
%         CodePrompt     = Code(svindex,(ceil(t_CodePrompt+0.05) + indx));
%         CodeLate       = Code(svindex,(ceil(t_CodeLate) + indx));
        t_CodeEarly_060   = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex); 
        t_CodeEarly_055   = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex);        
        t_CodeEarly       =   (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex); 
        t_CodeEarly_045   = (0 + Spacing(4) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(4) + remChip(svindex);        
        t_CodeEarly_040   = (0 + Spacing(5) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex)+ Spacing(5) + remChip(svindex);       
        t_CodeEarly_035   = (0 + Spacing(6) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(6) + remChip(svindex);       
        t_CodeEarly_030   = (0 + Spacing(7) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(7) + remChip(svindex);       
        t_CodeEarly_025   = (0 + Spacing(8) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(8) + remChip(svindex);       
        t_CodeEarly_020   = (0 + Spacing(9) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(9) + remChip(svindex);       
        t_CodeEarly_015   = (0 + Spacing(10) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(10) + remChip(svindex);       
        t_CodeEarly_010   = (0 + Spacing(11) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(11) + remChip(svindex);       
        t_CodeEarly_005   = (0 + Spacing(12) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(12) + remChip(svindex);
        t_CodePrompt      = (0 + Spacing(13) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(13) + remChip(svindex);
        t_CodeLate005     = (0 + Spacing(14) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(14) + remChip(svindex);
        t_CodeLate010     = (0 + Spacing(15) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(15) + remChip(svindex);
        t_CodeLate015     = (0 + Spacing(16) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(16) + remChip(svindex);
        t_CodeLate020     = (0 + Spacing(17) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(17) + remChip(svindex);
        t_CodeLate025     = (0 + Spacing(18) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(18) + remChip(svindex);
        t_CodeLate030     = (0 + Spacing(19) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(19) + remChip(svindex);
        t_CodeLate035     = (0 + Spacing(20) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(20) + remChip(svindex);
        t_CodeLate040     = (0 + Spacing(21) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(21) + remChip(svindex);
        t_CodeLate045     = (0 + Spacing(22) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(22) + remChip(svindex);
        t_CodeLate        = (0 + Spacing(23) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(23) + remChip(svindex);
        t_CodeLate055     = (0 + Spacing(24) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(24) + remChip(svindex);
        t_CodeLate060     = (0 + Spacing(25) + remChip(svindex)) : codePhaseStep(svindex) : (numSample -1) * codePhaseStep(svindex) + Spacing(25) + remChip(svindex);
        
        index = 2;
        CodeEarly_060 = Code(svindex,(ceil(t_CodeEarly_060) + index));
        CodeEarly_055 = Code(svindex,(ceil(t_CodeEarly_055) + index));
        CodeEarly         = Code(svindex,(ceil(t_CodeEarly) + index));
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

        
        %% Mixed with local carrier replica
        remChip(svindex) = t_CodePrompt(numSample) + codePhaseStep(svindex) - signal.codelength*pdi;
        
        CarrTime = (0:numSample)./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1), 2*pi);
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);        
        
        %%
        % Updated by LZD on 2021.10.07
%         E_i      = sum(CodeEarly    .*InphaseSignal);
%         E_q      = sum(CodeEarly    .*QuadratureSignal);
%         P_i      = sum(CodePrompt   .*InphaseSignal);
%         P_q      = sum(CodePrompt   .*QuadratureSignal);
%         L_i     = sum(CodeLate     .*InphaseSignal);
%         L_q     = sum(CodeLate     .*QuadratureSignal);
        
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


        %% Calculate CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_CT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        % Implement code loop filter and generate NCO command
        E = sqrt(E_i^2+E_q^2);
        L = sqrt(L_i^2+L_q^2);
        codeError(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator
        codeNco(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(codeError(svindex)...
            - DLLdiscriLast(svindex)) + codeError(svindex)* ((pdi*t)/tau1code);
        DLLdiscriLast(svindex) = codeError(svindex);
        code_outputLast(svindex) = codeNco(svindex);
        %         codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
        codeFreq(svindex) = signal.codeFreqBasis + codeNco(svindex);
        
        % PLL discriminator
        carrError(svindex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
        carrNco(svindex) = oldCarrNco(svindex) + (tau2carr/tau1carr)*(carrError(svindex) ...
            - oldCarrError(svindex)) + carrError(svindex) * ((pdi*t)/tau1carr);
        oldCarrNco(svindex) = carrNco(svindex);
        oldCarrError(svindex) = carrError(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrNco(svindex);  % Modify carrier freq
        
        %% Data Recording
        % Updated by LZD on 2021.10.07
%         TckResultCT(prn).E_i(msIndex) = E_i;
%         TckResultCT(prn).E_q(msIndex) = E_q;
%         TckResultCT(prn).P_i(msIndex) = P_i;
%         TckResultCT(prn).P_q(msIndex) = P_q;
%         TckResultCT(prn).L_i(msIndex) = L_i;
%         TckResultCT(prn).L_q(msIndex) = L_q;
        TckResultCT(prn).E_i_060(msIndex)         = E_i_060;
        TckResultCT(prn).E_i_055(msIndex)         = E_i_055;
        TckResultCT(prn).E_i(msIndex)            = E_i;
        TckResultCT(prn).E_i_045(msIndex)         = E_i_045;
        TckResultCT(prn).E_i_040(msIndex)         = E_i_040;
        TckResultCT(prn).E_i_035(msIndex)         = E_i_035;
        TckResultCT(prn).E_i_030(msIndex)         = E_i_030;
        TckResultCT(prn).E_i_025(msIndex)         = E_i_025;
        TckResultCT(prn).E_i_020(msIndex)         = E_i_020;
        TckResultCT(prn).E_i_015(msIndex)         = E_i_015;
        TckResultCT(prn).E_i_010(msIndex)         = E_i_010;
        TckResultCT(prn).E_i_005(msIndex)         = E_i_005;
        TckResultCT(prn).E_q_060(msIndex)         = E_q_060;
        TckResultCT(prn).E_q_055(msIndex)         = E_q_055;
        TckResultCT(prn).E_q(msIndex)            = E_q;
        TckResultCT(prn).E_q_045(msIndex)         = E_q_045;
        TckResultCT(prn).E_q_040(msIndex)         = E_q_040;
        TckResultCT(prn).E_q_035(msIndex)         = E_q_035;
        TckResultCT(prn).E_q_030(msIndex)         = E_q_030;
        TckResultCT(prn).E_q_025(msIndex)         = E_q_025;
        TckResultCT(prn).E_q_020(msIndex)         = E_q_020;
        TckResultCT(prn).E_q_015(msIndex)         = E_q_015;
        TckResultCT(prn).E_q_010(msIndex)         = E_q_010;
        TckResultCT(prn).E_q_005(msIndex)         = E_q_005;
        TckResultCT(prn).P_i(msIndex)                             = P_i;
        TckResultCT(prn).P_q(msIndex)                             = P_q; 
        TckResultCT(prn).L_i005(msIndex)            = L_i005;
        TckResultCT(prn).L_i010(msIndex)            = L_i010;
        TckResultCT(prn).L_i015(msIndex)            = L_i015;
        TckResultCT(prn).L_i020(msIndex)            = L_i020;
        TckResultCT(prn).L_i025(msIndex)            = L_i025;
        TckResultCT(prn).L_i030(msIndex)            = L_i030;
        TckResultCT(prn).L_i035(msIndex)            = L_i035;
        TckResultCT(prn).L_i040(msIndex)            = L_i040;
        TckResultCT(prn).L_i045(msIndex)            = L_i045;
        TckResultCT(prn).L_i(msIndex)            = L_i;
        TckResultCT(prn).L_i055(msIndex)            = L_i055;
        TckResultCT(prn).L_i060(msIndex)            = L_i060;
        TckResultCT(prn).L_q005(msIndex)            = L_q005;
        TckResultCT(prn).L_q010(msIndex)            = L_q010;
        TckResultCT(prn).L_q015(msIndex)            = L_q015;
        TckResultCT(prn).L_q020(msIndex)            = L_q020;
        TckResultCT(prn).L_q025(msIndex)            = L_q025;
        TckResultCT(prn).L_q030(msIndex)            = L_q030;
        TckResultCT(prn).L_q035(msIndex)            = L_q035;
        TckResultCT(prn).L_q040(msIndex)            = L_q040;
        TckResultCT(prn).L_q045(msIndex)            = L_q045;
        TckResultCT(prn).L_q(msIndex)            = L_q;
        TckResultCT(prn).L_q055(msIndex)            = L_q055; 
        TckResultCT(prn).L_q060(msIndex)            = L_q060;  

        
        
        
        TckResultCT(prn).carrError(msIndex)       = carrError(svindex);
        TckResultCT(prn).codeError(msIndex)       = codeError(svindex);
        TckResultCT(prn).codeFreq(msIndex)        = codeFreq(svindex);
        TckResultCT(prn).carrFreq(msIndex)        = carrFreq(svindex);
        TckResultCT(prn).numSample(msIndex)       = numSample;
        TckResultCT(prn).remChip(msIndex)         = remChip(svindex);
        TckResultCT(prn).remCarrPhase(msIndex)    = remCarrPhase(svindex);
        TckResultCT(prn).absoluteSample(msIndex)  = ftell(file.fid);
        TckResultCT(prn).absoluteSampleCodedelay(msIndex)  = mod(TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),fs*t );
        TckResultCT(prn).codedelay(msIndex)       = signal.Sample - AcqCodeDelay(svindex) +1 + sum(delayValue(svindex,(1:msIndex)));
        TckResultCT(prn).codedelay2(msIndex)      = mod( TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),fs*t );
        TckResultCT(prn).delayValue(msIndex)      = delayValue(svindex,msIndex);
        
        
        
    end % end for svindex in Tracking
    
    
    %%
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStartMea + measSampleStep*posIndex;
    
    for svIndex=1:svlength
        prn = sv(svIndex);
        if TckResultCT(prn).absoluteSample(msIndex) > currMeasSample
            flg_pos(svIndex) = 1;
        else
            flg_pos(svIndex) = 0;
        end
    end
    
    %%
    if sum(flg_pos) == svlength  
        posIndex = posIndex + 1;
        for svIndex=1:svlength
            prn = sv(svIndex);
            
            % Find index of I_P stream whose integration contains current
            % measurment point location
            for index = 1: length(TckResultCT(prn).absoluteSample)
                if(TckResultCT(prn).absoluteSample(index) > currMeasSample )
                    break
                end
            end
            index = index - 1;
            
            % Update the phasestep based on code freq and sampling frequency
            codePhaseStepX = TckResultCT(prn).codeFreq(index)/signal.Fs;
            
            codePhaseMeas(svIndex) = TckResultCT(prn).remChip(index) + ...
                codePhaseStepX*((currMeasSample - TckResultCT(prn).absoluteSample(index))/file.dataType);%file.dataType
            
            transmitTime(svIndex) = codePhaseMeas(svIndex)/signal.codelength/1000 + ... %TckResultCT(prn).codedelay(msIndex)/(signal.Fs/1000) 
                (index - (sbf.nav1(prn)+eph(prn).sfb(1)*20))/1000 + ...
                eph(prn).TOW(1);
        end
      
        % At first time of fix, local time is initialized by transmitTime and
        % settings.startOffset
        if (localTime == inf)
            maxTime   = max(transmitTime);
            localTime = maxTime + 75/1000; % 68 ms is an assumed travel time
        end
        pseudorange = (ones(1,svlength).*localTime - transmitTime)*cmn.cSpeed;
        
        %
        usr_clk = usr_clk_wls ;%%%
        estusr = estusr_wls;%%%
        
        for svindex = 1 : svlength
            prn = sv(svindex);
            
            tot_est_pos(svindex) = transmitTime(svindex);% ...
%                                                 + (1/cmn.cSpeed)*sv_clk(prn);
            
            % find the sv pos in ECEF at the time of transmision
            [svxyz(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel] = ...
                svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
            
            % C/A-code pseudorange corrected for satellite clock (in meters) and Tgd(in sec)
            prvec(svindex)      = pseudorange(svindex) + sv_clk(prn) - grpdel*cmn.cSpeed;% -sv_clk(prn)?
            
            % Adjust satellite position coordinates for earth rotation correction
            svxyzr(svindex,:)   = erotcorr(svxyz(svindex,:),prvec(svindex));
            
            % tropospheric and ionospheric delay correction
            counter_corr(svindex) = counter_corr(svindex) + 1;
            if counter_corr(svindex) ==  corrUpt
                svenu           = xyz2enu(svxyzr(svindex,:), estusr(1:3));%
                el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
                %             az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
                az_rad(svindex) =  atan2(svenu(1),svenu(2));
                az(svindex)     = az_rad(svindex)*180/pi;
                el(svindex)     = el_rad(svindex)*180/pi;
                temp            = xyz2llh(estusr(1:3));
                user_ll         = [temp(1:2).*180/pi temp(3)];
                ionodel(svindex)        = ionocorr(tot_est_pos(svindex),svxyzr(svindex,:), estusr(1:3));
                tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex)));
                counter_corr(svindex)   = 0;
            end
            
            prvec(svindex) = prvec(svindex) - ionodel(svindex) - tropodel_unb3(svindex); % sign of iono. and trop. error?
        end % for svindex=1:svlength
        
        
        %% Record Ppseudorange measurement 
        navSolutionsWLS.rawPseudorange(posIndex,:) = pseudorange ;        
        
        %% Position cal using LS method
        [estusr_wls, dop]       = olspos(prvec,svxyzr,estusr_wls); % ordinary least square
        [VR, dtRV, ~]     = ...
            LS_SA_code_Vel(estusr_wls(1:3)', svxyzr, sv_vel, carrFreq'-ones(length(carrFreq),1)*signal.IF, 0.190293672798365, sv_clk_vel(sv)); %carrFreq'-ones(length(carrFreq),1)*signal.IF
        
        usrenu_wls(posIndex,:)   = xyz2enu(estusr_wls(1:3),cnslxyz);
        usr_clk_wls             = estusr_wls(4);
        
        llh     = xyz2llh(estusr_wls(1:3));
        L_b     = llh(1);
        lamda_b = llh(2);
        C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
            -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
            -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];
        usr_velENU(posIndex,:) = C_e_n * VR ;
        
        usrenu_wls(posIndex,:)                 	= xyz2enu(estusr_wls(1:3),cnslxyz);
        usrllh_wls(posIndex,:)                   = xyz2llh(estusr_wls(1:3));
        usrllh_wls(posIndex,1:2)                 = usrllh_wls(posIndex,1:2)*180/pi;
        
        %     navSolutionsWLS.RxTime(msIndex)       = RxTime(msIndex);
        navSolutionsWLS.usrPos(posIndex,:)       = estusr_wls(1:3);
        navSolutionsWLS.usrVel(posIndex,:)       = VR;
        navSolutionsWLS.usrPosENU(posIndex,:)    = usrenu_wls(posIndex,:);
        navSolutionsWLS.usrPosLLH(posIndex,:)    = usrllh_wls(posIndex,:);
        navSolutionsWLS.clkBias(posIndex)        = usr_clk_wls;
        navSolutionsWLS.usrVelENU(posIndex,:)        = usr_velENU(posIndex,:);
        navSolutionsWLS.clkDrift(posIndex)   = dtRV; % m/s
        navSolutionsWLS.DOP(posIndex,:)       = dop;
        navSolutionsWLS.satEA(posIndex,:)      = el;
        navSolutionsWLS.satAZ(posIndex,:)      = az;
        navSolutionsWLS.timeTransmit(posIndex,:)      = transmitTime;
        
        navSolutionsWLS.codePhaseMeas(posIndex,:) = codePhaseMeas;
        %     navSolutionsWLS.test(msIndex,:)      = test;
        
        
        %=== Correct local time by clock error estimation =================
        localTime = localTime - navSolutionsWLS.clkBias(posIndex)/cmn.cSpeed;
        navSolutionsWLS.localTime(posIndex,:) = localTime;
        
        %=== Update local time by measurement  step  ====================================
        localTime = localTime + measSampleStep/signal.Fs;
        
       if posIndex == 20
           xx = 0;
       end
        navSolutionsCT = navSolutionsWLS;
        
        if mod(posIndex, 1) == 0
            fprintf('WLS: index = %4d; localTime: %f;  E = %f N = %f U = %f  VE = %f VN = %f VU = %f B = %f D = %f\n\n', ...
                posIndex, localTime, usrenu_wls(posIndex,1),usrenu_wls(posIndex,2),usrenu_wls(posIndex,3),usr_velENU(posIndex,1), usr_velENU(posIndex,2), usr_velENU(posIndex,3), usr_clk_wls, dtRV);
        end
    end % end for positioning at current measurement epoch
    
end % end for msIndex

close(h);

TckResultCT_mltCorr = TckResultCT;
navSolutionsCT_mltCorr = navSolutionsCT;

save(['navSolCT_',num2str(pdi),'ms_mltCorr_',file.fileName], 'navSolutionsCT_mltCorr' );
save(['tckRstCT_',num2str(pdi),'ms_mltCorr_',file.fileName], 'TckResultCT_mltCorr','CN0_CT');

