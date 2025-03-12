function  TckResultCT = trackingCT_multiCorr(file,signal,track,Acquired)
%Purpose:
%   Perform signal tracking using conventional DLL and PLL
%Inputs:
%	file        - parameters related to the data file to be processed
%	signal      - parameters related to signals,a structure
%	track       - parameters related to signal tracking 
%	Acquired    - acquisition results
%Outputs:
%	TckResultCT	- conventional tracking results, e.g. correlation values in 
%                   inphase prompt (P_i) channel and in qudrature prompt 
%                   channel (P_q), etc.
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
% 
% Copyright (C) X X
% Written by X X 

%%

 
                        
                        
                        
Spacing = -0.6:0.05:0.6;

[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

datalength = 50000 ;%track.msToProcessCT;
delayValue = zeros(length(Acquired.sv),datalength);

svlength    = length(Acquired.sv);
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);
pdi=1;
sv  = Acquired.sv;
for svindex = 1:length(Acquired.sv) 
    codetemp                = generateCAcode(Acquired.sv(svindex));
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    remChip(svindex) = 0;
    remPhase=0;
    remSample = 0;
    carrier_output=0;
    carrier_outputLast=0;
    PLLdiscriLast=0;
    code_output=0;
    code_outputLast=0;
    DLLdiscriLast=0;
    msIndex = 0;
    AcqDoppler = Acquired.fineFreq(svindex)-signal.IF;
    AcqCodeDelay = Acquired.codedelay(svindex);
    
    Codedelay = AcqCodeDelay;
    codeFreq(svindex) = signal.codeFreqBasis;
    carrierFreqBasis = Acquired.fineFreq(svindex);
    carrierFreq = Acquired.fineFreq(svindex);
    
    % set the file position indicator according to the acquired code delay
    fseek(file.fid,(signal.Sample-AcqCodeDelay-1+file.skip*signal.Sample)*file.dataPrecision*file.dataType,'bof');  %
%     fseek(file.fid,(AcqCodeDelay-1+file.skip*signal.Sample)*file.dataPrecision*file.dataType,'bof');  %
    
%     Code = generateCAcode(Acquired.sv(svindex));
%     Code = [Code(end) Code Code(1)];
    
    h = waitbar(0,['Ch:',num2str(svindex),' Please wait...']);
    
    for msIndex = 1: datalength
        
        waitbar(msIndex/datalength) 
        
        
        remSample = ((signal.codelength-remChip(svindex)) / (codeFreq(svindex)/signal.Fs));
        numSample = ceil((signal.codelength-remChip(svindex))/(codeFreq(svindex)/signal.Fs));%round
        delayValue(svindex,msIndex) = numSample - signal.Sample;
        
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')'; 
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal0DC = sin_rawsignal - mean(sin_rawsignal) + 1i*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal0DC = fread(file.fid,numSample*file.dataType,'int8')';  
            rawsignal0DC = rawsignal0DC(1:2:length(rawsignal0DC)) + 1i*rawsignal0DC(2:2:length(rawsignal0DC));% For NSL STEREO LBand only
        end
  
                
        %% spacing = -0.6:0.05:0.6
        codePhaseStep(svindex) = codeFreq(svindex)/signal.Fs;
        t_CodeEarly_060   = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex)); 
        t_CodeEarly_055   = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));        
        t_CodeEarly       = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex)); 
        t_CodeEarly_045   = (0 + Spacing(4) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(4) + remChip(svindex));        
        t_CodeEarly_040   = (0 + Spacing(5) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex)+ Spacing(5) + remChip(svindex));       
        t_CodeEarly_035   = (0 + Spacing(6) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(6) + remChip(svindex));       
        t_CodeEarly_030   = (0 + Spacing(7) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(7) + remChip(svindex));       
        t_CodeEarly_025   = (0 + Spacing(8) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(8) + remChip(svindex));       
        t_CodeEarly_020   = (0 + Spacing(9) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(9) + remChip(svindex));       
        t_CodeEarly_015   = (0 + Spacing(10) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(10) + remChip(svindex));       
        t_CodeEarly_010   = (0 + Spacing(11) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(11) + remChip(svindex));       
        t_CodeEarly_005   = (0 + Spacing(12) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(12) + remChip(svindex));
        t_CodePrompt      = (0 + Spacing(13) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(13) + remChip(svindex));
        t_CodeLate005     = (0 + Spacing(14) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(14) + remChip(svindex));
        t_CodeLate010     = (0 + Spacing(15) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(15) + remChip(svindex));
        t_CodeLate015     = (0 + Spacing(16) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(16) + remChip(svindex));
        t_CodeLate020     = (0 + Spacing(17) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(17) + remChip(svindex));
        t_CodeLate025     = (0 + Spacing(18) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(18) + remChip(svindex));
        t_CodeLate030     = (0 + Spacing(19) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(19) + remChip(svindex));
        t_CodeLate035     = (0 + Spacing(20) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(20) + remChip(svindex));
        t_CodeLate040     = (0 + Spacing(21) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(21) + remChip(svindex));
        t_CodeLate045     = (0 + Spacing(22) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(22) + remChip(svindex));
        t_CodeLate        = (0 + Spacing(23) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(23) + remChip(svindex));
        t_CodeLate055     = (0 + Spacing(24) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(24) + remChip(svindex));
        t_CodeLate060     = (0 + Spacing(25) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(25) + remChip(svindex));

        
        CodeEarly_060  = Code(svindex,(ceil(t_CodeEarly_060) + 1));
        CodeEarly_055  = Code(svindex,(ceil(t_CodeEarly_055) + 1));
        CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + 1));
        CodeEarly_045  = Code(svindex,(ceil(t_CodeEarly_045) + 1));
        CodeEarly_040  = Code(svindex,(ceil(t_CodeEarly_040) + 1));
        CodeEarly_035  = Code(svindex,(ceil(t_CodeEarly_035) + 1));
        CodeEarly_030  = Code(svindex,(ceil(t_CodeEarly_030) + 1));
        CodeEarly_025  = Code(svindex,(ceil(t_CodeEarly_025) + 1));
        CodeEarly_020  = Code(svindex,(ceil(t_CodeEarly_020) + 1));
        CodeEarly_015  = Code(svindex,(ceil(t_CodeEarly_015) + 1));
        CodeEarly_010  = Code(svindex,(ceil(t_CodeEarly_010) + 1));
        CodeEarly_005  = Code(svindex,(ceil(t_CodeEarly_005) + 1));
        CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate005    = Code(svindex,(ceil(t_CodeLate005) + 1));
        CodeLate010    = Code(svindex,(ceil(t_CodeLate010) + 1));
        CodeLate015    = Code(svindex,(ceil(t_CodeLate015) + 1));
        CodeLate020    = Code(svindex,(ceil(t_CodeLate020) + 1));
        CodeLate025    = Code(svindex,(ceil(t_CodeLate025) + 1));
        CodeLate030    = Code(svindex,(ceil(t_CodeLate030) + 1));
        CodeLate035    = Code(svindex,(ceil(t_CodeLate035) + 1));
        CodeLate040    = Code(svindex,(ceil(t_CodeLate040) + 1));
        CodeLate045    = Code(svindex,(ceil(t_CodeLate045) + 1));
        CodeLate       = Code(svindex,(ceil(t_CodeLate) + 1));
        CodeLate055    = Code(svindex,(ceil(t_CodeLate055) + 1));
        CodeLate060    = Code(svindex,(ceil(t_CodeLate060) + 1));
        
        
%         t_CodeEarly    = (0 + Spacing(1) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(1) + remChip);
%         t_CodePrompt   = (0 + Spacing(2) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(2) + remChip);
%         t_CodeLate     = (0 + Spacing(3) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(3) + remChip);
%         CodeEarly      = Code(ceil(t_CodeEarly) + 1);
%         CodePrompt     = Code(ceil(t_CodePrompt) + 1);
%         CodeLate       = Code(ceil(t_CodeLate) + 1);
        remChip(svindex)   = (t_CodePrompt(numSample) + codeFreq(svindex)/signal.Fs) - signal.codeFreqBasis*signal.ms;
        
        CarrTime = (0 : numSample)./signal.Fs;
        Wave     = (2*pi*(carrierFreq .* CarrTime)) + remPhase ;  
        remPhase =  rem( Wave(numSample+1), 2*pi); 
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal0DC .* carrsig);
        QuadratureSignal = real(rawsignal0DC .* carrsig);
        
%         E_i  = sum(CodeEarly    .*InphaseSignal);  E_q = sum(CodeEarly    .*QuadratureSignal);
%         P_i  = sum(CodePrompt   .*InphaseSignal);  P_q = sum(CodePrompt   .*QuadratureSignal);
%         L_i  = sum(CodeLate     .*InphaseSignal);  L_q = sum(CodeLate     .*QuadratureSignal);
     
        %%
        E_i_060  = sum(CodeEarly_060    .*InphaseSignal);
        E_q_060  = sum(CodeEarly_060    .*QuadratureSignal);
        E_i_055  = sum(CodeEarly_055    .*InphaseSignal);
        E_q_055  = sum(CodeEarly_055    .*QuadratureSignal);
        E_i      = sum(CodeEarly    .*InphaseSignal);  
        E_q      = sum(CodeEarly    .*QuadratureSignal);
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
        P_i      = sum(CodePrompt   .*InphaseSignal);  
        P_q      = sum(CodePrompt   .*QuadratureSignal);
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
        L_i     = sum(CodeLate     .*InphaseSignal); 
        L_q     = sum(CodeLate     .*QuadratureSignal); 
        L_i055  = sum(CodeLate055     .*InphaseSignal);  
        L_q055  = sum(CodeLate055     .*QuadratureSignal);
        L_i060  = sum(CodeLate060     .*InphaseSignal); 
        L_q060  = sum(CodeLate060     .*QuadratureSignal); 
  
        
        % Calculate CN0
        flag_snr(svindex)=1;
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_CT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*signal.ms*pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        
        % DLL
        E               = sqrt(E_i^2+E_q^2);
        L               = sqrt(L_i^2+L_q^2);
        DLLdiscri       = 0.5 * (E-L)/(E+L);
        code_output     = code_outputLast + (tau2code/tau1code)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/tau1code);
        DLLdiscriLast   = DLLdiscri;
        code_outputLast = code_output;
        codeFreq(svindex)        = signal.codeFreqBasis - code_output;
        
        % PLL
        PLLdiscri           = atan(P_q/P_i) / (2*pi);
        carrier_output      = carrier_outputLast + (tau2carr/tau1carr)*(PLLdiscri - PLLdiscriLast) + PLLdiscri * (0.001/tau1carr);
        carrier_outputLast  = carrier_output;  
        PLLdiscriLast       = PLLdiscri;
        carrierFreq         = carrierFreqBasis + carrier_output;  % Modify carrier freq based on NCO command
        
        % Data Record
        
        TckResultCT(Acquired.sv(svindex)).E_i_060(msIndex)         = E_i_060;
        TckResultCT(Acquired.sv(svindex)).E_i_055(msIndex)         = E_i_055;
        TckResultCT(Acquired.sv(svindex)).E_i(msIndex)            = E_i;
        TckResultCT(Acquired.sv(svindex)).E_i_045(msIndex)         = E_i_045;
        TckResultCT(Acquired.sv(svindex)).E_i_040(msIndex)         = E_i_040;
        TckResultCT(Acquired.sv(svindex)).E_i_035(msIndex)         = E_i_035;
        TckResultCT(Acquired.sv(svindex)).E_i_030(msIndex)         = E_i_030;
        TckResultCT(Acquired.sv(svindex)).E_i_025(msIndex)         = E_i_025;
        TckResultCT(Acquired.sv(svindex)).E_i_020(msIndex)         = E_i_020;
        TckResultCT(Acquired.sv(svindex)).E_i_015(msIndex)         = E_i_015;
        TckResultCT(Acquired.sv(svindex)).E_i_010(msIndex)         = E_i_010;
        TckResultCT(Acquired.sv(svindex)).E_i_005(msIndex)         = E_i_005;
        TckResultCT(Acquired.sv(svindex)).E_q_060(msIndex)         = E_q_060;
        TckResultCT(Acquired.sv(svindex)).E_q_055(msIndex)         = E_q_055;
        TckResultCT(Acquired.sv(svindex)).E_q(msIndex)            = E_q;
        TckResultCT(Acquired.sv(svindex)).E_q_045(msIndex)         = E_q_045;
        TckResultCT(Acquired.sv(svindex)).E_q_040(msIndex)         = E_q_040;
        TckResultCT(Acquired.sv(svindex)).E_q_035(msIndex)         = E_q_035;
        TckResultCT(Acquired.sv(svindex)).E_q_030(msIndex)         = E_q_030;
        TckResultCT(Acquired.sv(svindex)).E_q_025(msIndex)         = E_q_025;
        TckResultCT(Acquired.sv(svindex)).E_q_020(msIndex)         = E_q_020;
        TckResultCT(Acquired.sv(svindex)).E_q_015(msIndex)         = E_q_015;
        TckResultCT(Acquired.sv(svindex)).E_q_010(msIndex)         = E_q_010;
        TckResultCT(Acquired.sv(svindex)).E_q_005(msIndex)         = E_q_005;
        TckResultCT(Acquired.sv(svindex)).P_i(msIndex)                             = P_i;
        TckResultCT(Acquired.sv(svindex)).P_q(msIndex)                             = P_q; 
        TckResultCT(Acquired.sv(svindex)).L_i005(msIndex)            = L_i005;
        TckResultCT(Acquired.sv(svindex)).L_i010(msIndex)            = L_i010;
        TckResultCT(Acquired.sv(svindex)).L_i015(msIndex)            = L_i015;
        TckResultCT(Acquired.sv(svindex)).L_i020(msIndex)            = L_i020;
        TckResultCT(Acquired.sv(svindex)).L_i025(msIndex)            = L_i025;
        TckResultCT(Acquired.sv(svindex)).L_i030(msIndex)            = L_i030;
        TckResultCT(Acquired.sv(svindex)).L_i035(msIndex)            = L_i035;
        TckResultCT(Acquired.sv(svindex)).L_i040(msIndex)            = L_i040;
        TckResultCT(Acquired.sv(svindex)).L_i045(msIndex)            = L_i045;
        TckResultCT(Acquired.sv(svindex)).L_i(msIndex)            = L_i;
        TckResultCT(Acquired.sv(svindex)).L_i055(msIndex)            = L_i055;
        TckResultCT(Acquired.sv(svindex)).L_i060(msIndex)            = L_i060;
        TckResultCT(Acquired.sv(svindex)).L_q005(msIndex)            = L_q005;
        TckResultCT(Acquired.sv(svindex)).L_q010(msIndex)            = L_q010;
        TckResultCT(Acquired.sv(svindex)).L_q015(msIndex)            = L_q015;
        TckResultCT(Acquired.sv(svindex)).L_q020(msIndex)            = L_q020;
        TckResultCT(Acquired.sv(svindex)).L_q025(msIndex)            = L_q025;
        TckResultCT(Acquired.sv(svindex)).L_q030(msIndex)            = L_q030;
        TckResultCT(Acquired.sv(svindex)).L_q035(msIndex)            = L_q035;
        TckResultCT(Acquired.sv(svindex)).L_q040(msIndex)            = L_q040;
        TckResultCT(Acquired.sv(svindex)).L_q045(msIndex)            = L_q045;
        TckResultCT(Acquired.sv(svindex)).L_q(msIndex)            = L_q;
        TckResultCT(Acquired.sv(svindex)).L_q055(msIndex)            = L_q055; 
        TckResultCT(Acquired.sv(svindex)).L_q060(msIndex)            = L_q060;  
        
        
        TckResultCT(Acquired.sv(svindex)).PLLdiscri(msIndex)      = PLLdiscri;
        TckResultCT(Acquired.sv(svindex)).DLLdiscri(msIndex)      = DLLdiscri;
        TckResultCT(Acquired.sv(svindex)).codedelay(msIndex)      = Codedelay + sum(delayValue(1:msIndex));
        TckResultCT(Acquired.sv(svindex)).remChip(msIndex)        = remChip(svindex);
        TckResultCT(Acquired.sv(svindex)).codeFreq(msIndex)       = codeFreq(svindex);  
        TckResultCT(Acquired.sv(svindex)).carrierFreq(msIndex)    = carrierFreq;  
        TckResultCT(Acquired.sv(svindex)).remPhase(msIndex)       = remPhase;
        TckResultCT(Acquired.sv(svindex)).remSample(msIndex)      = remSample;
        TckResultCT(Acquired.sv(svindex)).numSample(msIndex)      = numSample;
        TckResultCT(Acquired.sv(svindex)).delayValue(msIndex)     = delayValue(svindex,msIndex);
    end
    close(h);
end % end for

save(['TckResultCT_multiCorr_',file.fileName], 'TckResultCT','CN0_CT');