function [pseudorange, relative_time]= pr_est_interp(eph, Acquired,for_prest,signal, codedelay_pos, pdi, msIndex)
%Purpose
%   pseudorange estimation
%Inputs: 
%	eph             - ephemeris 
%	Acquired        - acquisition results
%	for_prest     	- parameters used for pseudorange estimation
%	signal          - parameters related to signals
%	codedelay_pos 	- code delay used for positioning
%	pdi             - predetection integration interval
%	msIndex         - index of millisecond 
%Outputs:
%	pseudorange         - estimated pseudorange
%	relative_time   	- time difference between different channels
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
% 
% Copyright (C) X X
% Written by X X

%% 
sv      = Acquired.sv;
c       = 299792458;
bSec    = 75*1e-3; % base second: the distance between sv and usr is between 67 to 86 sec.
fs      = signal.Fs; 
ms      = 1e-3; 
finefreq = 1/fs;

for svindex = 1 : length(Acquired.sv)
    prn = sv(svindex);
    if length(for_prest.sfb1(prn)) > 0
        sfb1            = for_prest.sfb1(prn); % unit: 20ms, the beginning of subframe 1
        nav1            = for_prest.nav1(prn); % unit: ms, the first navigation data point     
%         bca             = TckResult(prn).codedelay(ms_pr) - 1; % unit: finefreq, beginning of the C/A code
%         dat(svindex)    = 20 * 16367.667 * sfb1 + 16367.667 * nav1 + bca; % unit: finefreq.
%         bca             = TckResult(prn).codedelay(ms_pr) - 1; %%%%%%%%%%%%%% unit: finefreq, beginning of the C/A code
        bca             = codedelay_pos(svindex) -1;  
        dat(svindex)    = 20*fs*ms*sfb1 + fs*ms*nav1 + bca; % unit: finefreq.
    end
end % end of svindex

base = min(dat);
for svindex = 1 : length(Acquired.sv)
    if length(for_prest.sfb1(Acquired.sv(svindex))) > 0
        diff_of_dat(svindex)   = dat(svindex) - base;
        pseudorange(svindex)   = c * (bSec + diff_of_dat(svindex) * finefreq);
%         pseudorange(svindex)   = c*(bSec + (dat(svindex)* finefreq));
        relative_time(svindex) = diff_of_dat(svindex) * finefreq;
    end
end

% end % end of function