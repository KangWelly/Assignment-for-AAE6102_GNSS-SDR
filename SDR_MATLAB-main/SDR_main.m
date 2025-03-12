%Purpose:
%   Main function of the software-defined radio (SDR) receiver platform
%
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.1
% 
% Copyright (C) X X  
% Written by X X
 
clear; 
clc;
format long g;
addpath geo             %  
addpath acqtckpos       % Acquisition, tracking, and postiong calculation functions

%% Parameter initialization 
[file, signal, acq, track, solu, cmn] = initParameters();
 

%% Acquisition 
if ~exist(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat'])
    Acquired = acquisition(file,signal,acq); %
    save(['Acquired_',file.fileName,'_',num2str(file.skip)],'Acquired');    
else
    load(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat']);
    showAcqResult(Acquired, signal);
end 
if isempty(Acquired.sv)
    fprintf('No satellites acquired. Check parameter settings. \n ');
    return;
end

%% Do conventional signal tracking and obtain satellites ephemeris
if ~exist(['eph_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat'])
    % tracking using conventional DLL and PLL
    if ~exist(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']) %
        fprintf('Tracking for navigation data decoding ... \n\n');
        [TckResultCT, CN0_Eph,countinx] =  trackingCT(file,signal,track,Acquired); 
        TckResult_Eph = TckResultCT;
%           [TckResultCT_20ms, CN0_Eph] =  trackingCT_20ms(file,signal,track,Acquired);
%           TckResult_Eph = TckResultCT_20ms;
        save(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000)], 'TckResult_Eph','CN0_Eph','countinx');        
    else   
        load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']);
        %%%%%%%%%
    end 
    if isempty(TckResult_Eph)
        fprintf('Not enough raw data for navigation data decoding.  \n\n ');
        return
    end
    
    % navigaion data decoding
    fprintf('Navigation data decoding ... \n\n');
    [eph, ~, sbf] = naviDecode_updated(Acquired, TckResult_Eph);
    save(['eph_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000)], 'eph');
    save(['sbf_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000)], 'sbf');
else
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']);
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']);
    load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']);
end 
 
  
%% Find satellites that can be used to calculate user position
posSV  = findPosSV(file,Acquired,eph);


%% Do positiong in conventional or vector tracking mode
fprintf('Positioning ...\n\n');
cnslxyz = llh2xyz(solu.iniPos); % initial position in ECEF coordinate
 
disp(['Trying to load: navSolCT_1ms_', file.fileName, '.mat']);
if cmn.vtEnable == 1    
    fprintf('Positioning (VTL) ... \n\n');
  
    % load data to initilize VT
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;  
    
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']); % load eph
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT_10ms/1000),'.mat']); % 
    
    load(['tckRstCT_10ms_',file.fileName,'.mat']);%,'_Grid'
    load(['navSolCT_10ms_',file.fileName,'.mat']); 

    disp(['cmn.vtEnable = ', num2str(cmn.vtEnable)]);
    disp(['cmn.mltCorrON(2) = ', num2str(cmn.mltCorrON(2))]);


%%~~~~~VT/Multi~~~~~%%
    if cmn.mltCorrON(2) == 1 
%         load(['tckRstCT_1ms_',file.fileName,'.mat']);
%         load(['navSolCT_1ms',file.fileName,'.mat']);
        [TckResultVT_mltCorr, navSolutionsVT_mltCorr] =trackingVT_POS_updated_multicorrelator(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);                
    else
        if cmn.mltCorrON(2) == 0 
%                 load(['tckRstCT_1ms_',file.fileName,'_updated.mat']);
%                 load(['navSolCT_1ms',file.fileName,'_updated.mat']);
                [TckResultVT, navSolutionsVT] =trackingVT_POS_updated(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);  
        end
    end   
%%~~~~~~~~~~~~~~~%%

else 
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;
    
%%~~~~~CT/Multi~~~~~%%
    if cmn.mltCorrON(1) == 1
        [TckResultCT_mltCorr, navSolutionsCT_mltCorr] = trackingCT_POS_updated_multicorrelator(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph); %trackingCT_POS_multiCorr_1ms    
        load(['tckRstCT_10ms_mltCorr_',file.fileName,'.mat']);
        load(['navSolCT_10ms_mltCorr_',file.fileName,'.mat']);
    else
         if cmn.mltCorrON(1) == 0
            [TckResunavSolutionsCTltCT_pos, navSolutionsCT] =trackingCT_POS_updated(file,signal,track,cmn,Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu); %trackingCT_POS_multiCorr_1ms
            load(['tckRstCT_10ms_',file.fileName,'.mat']);
            load(['navSolCT_10ms_',file.fileName,'.mat']);
         end
    end
%%~~~~~~~~~~~~~~%%
end 

fprintf('Tracking and Positioning Completed.\n\n');


% Updated by LZD on 10.08
if cmn.vtEnable ==1 && cmn.mltCorrON(2) ==1
    load(['tckRstVT_',file.fileName,'_updated.mat']);
    load(['navSolVT_',file.fileName,'_updated.mat']);
    load(['tckRstVT_mltCorr_',file.fileName,'_updated.mat']);
    load(['navSolVT_mltCorr_',file.fileName,'_updated.mat']);
    
    load(['tckRstCT_10ms_',file.fileName,'.mat']);
    load(['navSolCT_10ms_',file.fileName,'.mat']);
    load(['tckRstCT_10ms_mltCorr_',file.fileName,'.mat']);
    load(['navSolCT_10ms_mltCorr_',file.fileName,'.mat']);
    
    save([file.fileName, '_', cmn.equip]);
    fprintf('The final results have been saved \n\n');
end
 



%     [TckResultVT_mltCorr, navSolutionsVT_mltCorr] =trackingVT_POS_updated_multicorrelator(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);                
%     [TckResultVT, navSolutionsVT] =trackingVT_POS_updated(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);  
%     [TckResultCT_mltCorr, navSolutionsCT_mltCorr] = trackingCT_POS_updated_multicorrelator(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph); %trackingCT_POS_multiCorr_1ms    
%     [TckResultCT_pos, navSolutionsCT] =trackingCT_POS_updated(file,signal,track,cmn,Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu); %trackingCT_POS_multiCorr_1ms


