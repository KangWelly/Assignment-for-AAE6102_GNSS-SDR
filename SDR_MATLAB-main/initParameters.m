function [file, signal, acq, track, solu, cmn] = initParameters()
%Purpose:
%   Parameter initialization
%Inputs: 
%	None
%Outputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%	acq         - parameters related to signal acquisition,a structure
%	track       - parameters related to signal tracking,a structure
%	solu        - parameters related to navigation solution,a structure
%	cmn         - parameters commmonly used,a structure
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.1
% 
% Copyright (C) X X
% Written by X X

%% File parameters
file.fileName       = 'Opensky';%'File_130'; 
file.fileRoute      = ['C:\Users\welly\Downloads\',file.fileName,'.bin'];  
file.skip        	= 5000; % in unit of ms
solu.iniPos	= [22.328444770087565/180 * pi, 114.1713630049711/180 * pi, 4]; 
%22.304015/180 * pi, 114.178904/180 * pi, 11 171-173
%22.303851/180 * pi, 114.180020/180 * pi, 12 174-176
%22.299856/180 * pi, 114.1800424/180 * pi, 4 178-180
% Ground truth location 

global ALPHA BETA % Parameters for iono. and trop. correction; From RINEX file
ALPHA = [9.3132E-09  1.4901e-08 -5.9605e-08 -1.1921e-07 ];
BETA  = [8.8064e+04  4.9152e+04 -1.3107e+05 -3.2768e+05];
cmn.doy = 171; % Day of year

%% File parameters
file.fid           	= fopen(file.fileRoute,'r','ieee-le');
file.skiptimeVT     = 100; % skip time from the first measurement epoch of CT, in uint of msec
file.dataType       = 2;    %1:I; 2:IQ
file.dataPrecision  = 1;    %1:int8 or byte; 2; int16 

%% Signal parameters
signal.IF               =4.58e6;%1580e6-1575.42e6;%0;   % unit: Hz 
signal.Fs               =58e6;%58e6;  	% unit: Hz
signal.Fc               = 1575.42e6; % unit: Hz	
signal.codeFreqBasis	= 1.023e6; % unit: Hz 	
signal.ms               = 1e-3; % unit: s
signal.Sample           = ceil(signal.Fs*signal.ms);	
signal.codelength       = signal.codeFreqBasis * signal.ms;

%% Acquisition parameters
acq.prnList     = 1:32;	% PRN list
acq.freqStep    = 500;	% unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.datalen     = 20; % msec
acq.L           = 10; % number of ms to perform FFT

%% Tracking parameters
% track.mode                  = 0;    % 0:conventional tracking; 1:vector tracking
track.CorrelatorSpacing  	= 0.5;  % unit: chip
track.DLLBW               	= 2;	% unit: Hz
track.DLLDamp           	= 0.707; 
track.DLLGain            	= 0.1;	
track.PLLBW              	= 15;
track.PLLDamp             	= 0.707;
track.PLLGain              	= 0.25; 	
track.msToProcessCT_1ms     = 1000; % unit: ms 
track.msToProcessCT_10ms    =40000;%unit: ms
track.ctPOS               = 3000; % unit: number of index
track.msToProcessVT         = 5000; %track.msPosCT - file.skiptimeVT; %
track.pdi                   = 1; % By now, not supporting other coherent integration time like 20 msec


%% Navigation solution parameters
solu.navSolPeriod = 20;     % unit: ms 
solu.mode  	= 2;    % 0:conventional LS/WLS; 1:conventionalKF; 2:VT


%% commonly used parameters
cmn.vtEnable  	=1;%   % 0: disable vector tracking; 1:enable vector tracking
cmn.mltCorrON = [1 0]; %% 0: disable multicorr; 1: enable mltCorr; cmn.mltCorrON(1): CT;    %cmn.mltCorrON(1): CT   ;  cmn.mltCorrON(2): VT  
                                            % Suggested STEPS:
                                            %"cmn.vtEnable = 0"-> [0 0]->[1 0]
                                            %"cmn.vtEnable = 1"-> [1 0]->[1 1]
cmn.cSpeed      = 299792458;    % speed of light, [m/s]
cmn.equip = 'stereo';   % labsat/stereo