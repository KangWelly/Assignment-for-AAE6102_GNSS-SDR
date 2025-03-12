function [TckResultCT_pos, navSolutionsCT] = trackingCT_POS_updated(file,signal,track,cmn, Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu)
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
load(['countinx'],'countinx');
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
pdi     = 1;%track.pdi;
t       = signal.ms;
svlength    = length(sv);
datalength  = track.ctPOS;%3000;

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
for svIndex = 1:length(sv)
    prn                     = sv(svIndex);
    codetemp                = generateCAcode(prn);
    Code(svIndex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    AcqCodeDelay(svIndex)   = Acquired.codedelay(svIndex);
    %     file_ptr(svIndex)       = signal.Sample - AcqCodeDelay(svIndex) -1 + file.skip *fs*t + 0;
    %     file_ptr(svIndex)       = signal.Sample - AcqCodeDelay(svIndex) +1 + file.skip *fs*t + 0;   % 29/04/2020
    
    % 29/04/29, for each channel, move the file pointer to the first
    % subframe (not subframe 1), to ensure 20 ms T_coh and do positioning
    %countinx=TckResult_Eph.countbit(svIndex);
    %fseek(file.fid,(signal.Sample-AcqCodeDelay+1 + file.skip*signal.Sample)*file.dataPrecision*file.dataType,'bof');
       %@@@@bit transition form TracResultCT
%     file_ptr(svIndex)       = (signal.Sample - AcqCodeDelay(svIndex) +1 ...
%         + (file.skip+countinx(svIndex)) *fs*t ... 
%         )*file.dataPrecision*file.dataType;
    file_ptr(svIndex)       = (signal.Sample - AcqCodeDelay(svIndex) +1 ...
        + (file.skip) *fs*t ... 
        )*file.dataPrecision*file.dataType;
    %file_ptr(svIndex)= TckResult_Eph(prn).absoluteSample(countinx(svIndex)-1);
    
    carrFreq(svIndex)       = Acquired.fineFreq(svIndex);
    AcqFreq(svIndex)        = Acquired.fineFreq(svIndex);
    
    oldcodedelay_pos(svIndex) = 0;
    oldabsoluteSample_pos(svIndex) = 0;
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
delayValue          = zeros(svlength,datalength);

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
Index=0;
for msIndex = 1: datalength/1 % Note that for pdi > 1ms, the index is still denoted as msIndex. 30/04/2020, BING XU
    waitbar(msIndex/(datalength/1),h)
    Index=Index+1;
        for svIndex = 1 :svlength
            if (msIndex<=(track.msToProcessCT_1ms+countinx(svIndex)))
                pdi=1;
                prn = sv(svIndex);
                
                % read raw data file
                codePhaseStep(svIndex) = codeFreq(svIndex)/signal.Fs;
                numSample = ceil((signal.codelength*pdi-remChip(svIndex))/codePhaseStep(svIndex));
                
                delayValue(svIndex,Index) = numSample - signal.Sample*pdi;
                
                fseek(file.fid, file_ptr(svIndex),'bof');
                
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
                
                file_ptr(svIndex)   = file_ptr(svIndex) + numSample*file.dataType;  %%%%%%
                
                %% spacing = -0.6:0.05:0.6
                t_CodeEarly       = (0 + Spacing(3) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(3) + remChip(svIndex));
                t_CodePrompt      = (0 + Spacing(13) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(13) + remChip(svIndex));
                t_CodeLate        = (0 + Spacing(23) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(23) + remChip(svIndex));
                
                indx = 1;
                CodeEarly      = Code(svIndex,(ceil(t_CodeEarly) + indx));
                CodePrompt     = Code(svIndex,(ceil(t_CodePrompt+0.05) + indx));
                CodeLate       = Code(svIndex,(ceil(t_CodeLate) + indx));
                
                %%
                remChip(svIndex) = t_CodePrompt(numSample) + codePhaseStep(svIndex) - signal.codelength*pdi;
                
                CarrTime = (0:numSample)./signal.Fs;
                Wave = 2*pi*((carrFreq(svIndex)).*CarrTime) + remCarrPhase(svIndex);
                remCarrPhase(svIndex) = rem(Wave(numSample+1), 2*pi);
                carrsig = exp(1i.* Wave(1:numSample));
                InphaseSignal    = imag(rawsignal .* carrsig);
                QuadratureSignal = real(rawsignal .* carrsig);
                
                %%
                E_i      = sum(CodeEarly    .*InphaseSignal);
                E_q      = sum(CodeEarly    .*QuadratureSignal);
                P_i      = sum(CodePrompt   .*InphaseSignal);
                P_q      = sum(CodePrompt   .*QuadratureSignal);
                L_i     = sum(CodeLate     .*InphaseSignal);
                L_q     = sum(CodeLate     .*QuadratureSignal);
                
                % Calculate CN0
                if (flag_snr(svIndex) == 1)
                    index_int(svIndex) = index_int(svIndex) + 1;
                    Zk(svIndex,index_int(svIndex)) = P_i^2 + P_q^2;
                    if mod(index_int(svIndex),K) == 0
                        meanZk  = mean(Zk(svIndex,:));
                        varZk   = var(Zk(svIndex,:));
                        NA2     = sqrt(meanZk^2-varZk);
                        varIQ   = 0.5 * (meanZk - NA2);
                        CN0_CT(snrIndex(svIndex),svIndex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                        index_int(svIndex)  = 0;
                        snrIndex(svIndex)   = snrIndex(svIndex) + 1;
                    end
                end
                
                % Implement code loop filter and generate NCO command
                E = sqrt(E_i^2+E_q^2);
                L = sqrt(L_i^2+L_q^2);
                codeError(svIndex) = 0.5*(E-L)/(E+L);  % DLL discriminator
                %code_output     = code_outputLast + (tau2code/tau1code)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/tau1code);
                codeNco(svIndex) = code_outputLast(svIndex) + (tau2code/tau1code)*(codeError(svIndex)...
                    - DLLdiscriLast(svIndex)) + codeError(svIndex)* (t/tau1code);
                DLLdiscriLast(svIndex) = codeError(svIndex);
                code_outputLast(svIndex) = codeNco(svIndex);
                %         codeFreq(svIndex) = signal.codeFreqBasis - codeNco(svIndex);
                codeFreq(svIndex) = signal.codeFreqBasis + codeNco(svIndex);
                
                % PLL discriminator
                carrError(svIndex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
                carrNco(svIndex) = oldCarrNco(svIndex) + (tau2carr/tau1carr)*(carrError(svIndex) ...
                    - oldCarrError(svIndex)) + carrError(svIndex) * (t/tau1carr);
                oldCarrNco(svIndex) = carrNco(svIndex);
                oldCarrError(svIndex) = carrError(svIndex);
                carrFreq(svIndex)  = AcqFreq(svIndex) + carrNco(svIndex);  % Modify carrier freq
                
                %% Data Recording
                TckResultCT(prn).E_i(Index) = E_i;
                TckResultCT(prn).E_q(Index) = E_q;
                TckResultCT(prn).P_i(Index) = P_i;
                TckResultCT(prn).P_q(Index) = P_q;
                TckResultCT(prn).L_i(Index) = L_i;
                TckResultCT(prn).L_q(Index) = L_q;
                
                
                TckResultCT(prn).carrError(Index)       = carrError(svIndex);
                TckResultCT(prn).codeError(Index)       = codeError(svIndex);
                TckResultCT(prn).codeFreq(Index)        = codeFreq(svIndex);
                TckResultCT(prn).carrFreq(Index)        = carrFreq(svIndex);
                TckResultCT(prn).numSample(Index)       = numSample;
                TckResultCT(prn).remChip(Index)         = remChip(svIndex);
                TckResultCT(prn).remCarrPhase(Index)    = remCarrPhase(svIndex);
                TckResultCT(prn).absoluteSample(Index)  = ftell(file.fid);
                TckResultCT(prn).absoluteSampleCodedelay(Index)  = mod(TckResultCT(prn).absoluteSample(Index)/(file.dataPrecision*file.dataType),fs*t );
                TckResultCT(prn).codedelay(Index)       = signal.Sample - AcqCodeDelay(svIndex) +1 + sum(delayValue(svIndex,(1:Index)));
                TckResultCT(prn).codedelay2(Index)      = mod( TckResultCT(prn).absoluteSample(Index)/(file.dataPrecision*file.dataType),fs*t );
                TckResultCT(prn).delayValue(Index)      = delayValue(svIndex,Index);
                
            elseif  (msIndex>(track.msToProcessCT_1ms+countinx(svIndex)))
                pdi=10;
                prn = sv(svIndex);
                Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model*pdi*t;
                thresUptR = 20/pdi;
                codetemp                = generateCAcode(prn);
                Code2(svIndex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
                corrUpt         = corrUpdateSec/(pdi*t);
                counter_corr    = corrUpt-1 * ones(svlength,1);
                % read raw data file
                codePhaseStep(svIndex) = codeFreq(svIndex)/signal.Fs;
                numSample = ceil((signal.codelength*pdi-remChip(svIndex))/codePhaseStep(svIndex));
                
                delayValue(svIndex,Index) = numSample - signal.Sample*pdi;
                
                fseek(file.fid, file_ptr(svIndex),'bof');
                
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
                
                file_ptr(svIndex)   = file_ptr(svIndex) + numSample*file.dataType;  %%%%%%
                
                %% spacing = -0.6:0.05:0.6
                t_CodeEarly       = (0 + Spacing(3) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(3) + remChip(svIndex));
                t_CodePrompt      = (0 + Spacing(13) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(13) + remChip(svIndex));
                t_CodeLate        = (0 + Spacing(23) + remChip(svIndex)) : codePhaseStep(svIndex) : ((numSample -1) * codePhaseStep(svIndex) + Spacing(23) + remChip(svIndex));
                
                indx = 1;
                CodeEarly      = Code2(svIndex,(ceil(t_CodeEarly) + indx));
                CodePrompt     = Code2(svIndex,(ceil(t_CodePrompt+0.05) + indx));
                CodeLate       = Code2(svIndex,(ceil(t_CodeLate) + indx));
                
                %%
                remChip(svIndex) = t_CodePrompt(numSample) + codePhaseStep(svIndex) - signal.codelength*pdi;
                
                CarrTime = (0:numSample)./signal.Fs;
                Wave = 2*pi*((carrFreq(svIndex)).*CarrTime) + remCarrPhase(svIndex);
                remCarrPhase(svIndex) = rem(Wave(numSample+1), 2*pi);
                carrsig = exp(1i.* Wave(1:numSample));
                InphaseSignal    = imag(rawsignal .* carrsig);
                QuadratureSignal = real(rawsignal .* carrsig);
                
                %%
                E_i      = sum(CodeEarly    .*InphaseSignal);
                E_q      = sum(CodeEarly    .*QuadratureSignal);
                P_i      = sum(CodePrompt   .*InphaseSignal);
                P_q      = sum(CodePrompt   .*QuadratureSignal);
                L_i     = sum(CodeLate     .*InphaseSignal);
                L_q     = sum(CodeLate     .*QuadratureSignal);
                
                % Calculate CN0
                if (flag_snr(svIndex) == 1)
                    index_int(svIndex) = index_int(svIndex) + 1;
                    Zk(svIndex,index_int(svIndex)) = P_i^2 + P_q^2;
                    if mod(index_int(svIndex),K) == 0
                        meanZk  = mean(Zk(svIndex,:));
                        varZk   = var(Zk(svIndex,:));
                        NA2     = sqrt(meanZk^2-varZk);
                        varIQ   = 0.5 * (meanZk - NA2);
                        CN0_CT(snrIndex(svIndex),svIndex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                        index_int(svIndex)  = 0;
                        snrIndex(svIndex)   = snrIndex(svIndex) + 1;
                    end
                end
                
                % Implement code loop filter and generate NCO command
                E = sqrt(E_i^2+E_q^2);
                L = sqrt(L_i^2+L_q^2);
                codeError(svIndex) = 0.5*(E-L)/(E+L);  % DLL discriminator
                %code_output     = code_outputLast + (tau2code/tau1code)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/tau1code);
                codeNco(svIndex) = code_outputLast(svIndex) + (tau2code/tau1code)*(codeError(svIndex)...
                    - DLLdiscriLast(svIndex)) + codeError(svIndex)* (t/tau1code);
                DLLdiscriLast(svIndex) = codeError(svIndex);
                code_outputLast(svIndex) = codeNco(svIndex);
                %         codeFreq(svIndex) = signal.codeFreqBasis - codeNco(svIndex);
                codeFreq(svIndex) = signal.codeFreqBasis + codeNco(svIndex);
                
                % PLL discriminator
                carrError(svIndex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
                carrNco(svIndex) = oldCarrNco(svIndex) + (tau2carr/tau1carr)*(carrError(svIndex) ...
                    - oldCarrError(svIndex)) + carrError(svIndex) * (t/tau1carr);
                oldCarrNco(svIndex) = carrNco(svIndex);
                oldCarrError(svIndex) = carrError(svIndex);
                carrFreq(svIndex)  = AcqFreq(svIndex) + carrNco(svIndex);  % Modify carrier freq
                
                %% Data Recording
                TckResultCT(prn).E_i(Index) = E_i;
                TckResultCT(prn).E_q(Index) = E_q;
                TckResultCT(prn).P_i(Index) = P_i;
                TckResultCT(prn).P_q(Index) = P_q;
                TckResultCT(prn).L_i(Index) = L_i;
                TckResultCT(prn).L_q(Index) = L_q;
                
                
                TckResultCT(prn).carrError(Index)       = carrError(svIndex);
                TckResultCT(prn).codeError(Index)       = codeError(svIndex);
                TckResultCT(prn).codeFreq(Index)        = codeFreq(svIndex);
                TckResultCT(prn).carrFreq(Index)        = carrFreq(svIndex);
                TckResultCT(prn).numSample(Index)       = numSample;
                TckResultCT(prn).remChip(Index)         = remChip(svIndex);
                TckResultCT(prn).remCarrPhase(Index)    = remCarrPhase(svIndex);
                TckResultCT(prn).absoluteSample(Index)  = ftell(file.fid);
                TckResultCT(prn).absoluteSampleCodedelay(Index)  = mod(TckResultCT(prn).absoluteSample(Index)/(file.dataPrecision*file.dataType),fs*t );
                TckResultCT(prn).codedelay(Index)       = signal.Sample - AcqCodeDelay(svIndex) +1 + sum(delayValue(svIndex,(1:Index)));
                TckResultCT(prn).codedelay2(Index)      = mod( TckResultCT(prn).absoluteSample(Index)/(file.dataPrecision*file.dataType),fs*t );
                TckResultCT(prn).delayValue(Index)      = delayValue(svIndex,Index);
                
                
                %
                
            end % end for svIndex in Tracking
         end    
       %end % end for svIndex in Tracking
    
        
    
    
    %%
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStartMea + measSampleStep*posIndex;
    
    for svIndex=1:svlength
        prn = sv(svIndex);
        if TckResultCT(prn).absoluteSample(Index) > currMeasSample
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
                ( (index-track.msToProcessCT_1ms-countinx(svIndex))*pdi + (track.msToProcessCT_1ms+countinx(svIndex))  - (sbf.nav1(prn)+eph(prn).sfb(1)*20))/1000 + ...
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
        
        for svIndex = 1 : svlength
            prn = sv(svIndex);
            
            tot_est_pos(svIndex) = transmitTime(svIndex);% ...
%                                                 + (1/cmn.cSpeed)*sv_clk(prn);
            
            % find the sv pos in ECEF at the time of transmision
            [svxyz(svIndex,:), sv_vel(svIndex,:), sv_clk(prn), sv_clk_vel(prn), grpdel] = ...
                svPosVel(prn,eph,tot_est_pos(svIndex),eph_idx(svIndex));
            
            % C/A-code pseudorange corrected for satellite clock (in meters) and Tgd(in sec)
            prvec(svIndex)      = pseudorange(svIndex) + sv_clk(prn) - grpdel*cmn.cSpeed;% -sv_clk(prn)?
            
            % Adjust satellite position coordinates for earth rotation correction
            svxyzr(svIndex,:)   = erotcorr(svxyz(svIndex,:),prvec(svIndex));
            
            % tropospheric and ionospheric delay correction
            counter_corr(svIndex) = counter_corr(svIndex) + 1;
            if counter_corr(svIndex) ==  corrUpt
                svenu           = xyz2enu(svxyzr(svIndex,:), estusr(1:3));%
                el_rad(svIndex) = atan(svenu(3)/norm(svenu(1:2)));
                %             az_rad(svIndex) = (pi/2)-atan2(svenu(1),svenu(2));
                az_rad(svIndex) =  atan2(svenu(1),svenu(2));
                az(svIndex)     = az_rad(svIndex)*180/pi;
                el(svIndex)     = el_rad(svIndex)*180/pi;
                temp            = xyz2llh(estusr(1:3));
                user_ll         = [temp(1:2).*180/pi temp(3)];
                ionodel(svIndex)        = ionocorr(tot_est_pos(svIndex),svxyzr(svIndex,:), estusr(1:3));
                tropodel_unb3(svIndex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svIndex)));
                counter_corr(svIndex)   = 0;
            end
            
            prvec(svIndex) = prvec(svIndex) - ionodel(svIndex) - tropodel_unb3(svIndex); % sign of iono. and trop. error?
        end % for svIndex=1:svlength
        
        
        %% Record Ppseudorange measurement 
        navSolutionsWLS.rawPseudorange(posIndex,:) = pseudorange ;        
        
        %% Position cal using LS method
        [estusr_wls, dop]       = olspos(prvec,svxyzr,estusr_wls); % ordinary least square
        [VR, dtRV, ~]     = ...
            LS_SA_code_Vel(estusr_wls(1:3)', svxyzr, sv_vel, carrFreq'-ones(length(carrFreq),1)*signal.IF, 0.190293672798365, sv_clk_vel(sv)); %carrFreq'-ones(length(carrFreq),1)*signal.IF
        %LS_SA_code_Vel(estusr_wls(1:3)', svxyzr, sv_vel, carrFreq'-ones(length(carrFreq),1)*signal.IF, 0.190293672798365, sv_clk_vel(sv));
        
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

TckResultCT_pos = TckResultCT;

save(['navSolCT_',num2str(pdi),'ms_',file.fileName], 'navSolutionsCT' );
save(['tckRstCT_',num2str(pdi),'ms_',file.fileName], 'TckResultCT_pos','CN0_CT');