
% clear

%% init parameters
% load 'TckResultCT_multiCorr_Data_20190116_data1.mat'
% TckResults_multiCorr = TckResultCT; 
% prnList = [3, 8, 14, 16, 22, 23, 26, 27, 31, 32];
% eleList = [24.96, 8.78, 50.45, 62.68, 25.81, 17.09, 62.89, 33.56, 42.7, 31.0];
% Data_20190116_data1_labeled = zeros(1,9); 

% load 'TckResultCT_multiCorr_Data_20190116_data4.mat'
% TckResults_multiCorr = TckResultCT; 
% prnList = [3, 8, 9, 16, 22, 23, 26, 27, 31, 32];
% eleList = [21.5, 21.26, 37.33, 62.38, 19.93, 28.32, 52.95, 48.22, 39.68, 18.53];
% Data_20190116_data4_labeled = zeros(1,8); 
% 
% 
% load 'TckResultCT_multiCorr_Data_20190116_data11.mat'
% TckResults_multiCorr = TckResultCT; 
% prnList = [3, 8, 14, 16, 22, 23, 26, 27, 31];
% eleList = [8.36, 54.23, 24.23, 47.1, 6.39, 50.04, 32.49, 83.8, 27.88];
% Data_20190116_data11_labeled = zeros(1,8); 
 


% load 'TckResultCT_multiCorr_Data_CENTRAL_static_c.mat'
% TckResults_multiCorr = TckResultCT; 
% % prnList = [3, 22, 23, 26, 32];
% % eleList = [22.503, 21.538, 25.865, 53.744, 19.663];
% prnList = [ 32];
% eleList = [ 19.663];
% Data_CENTRAL_static_c_labeled = zeros(1,11); %%%%%%%%%%%%%%%


% load 'TckRstct_data_static_TSTE_c_40.mat'
% TckResults_multiCorr = TckResultCT; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [5];
% eleList = [58.69];
% data_static_TSTE_c_40_labeled = zeros(1,8); 

% load 'TckRstct_data_static_TSTE_c_40.mat'
% TckResults_multiCorr = TckResultCT_multiCorr; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [10];
% eleList = [58.69];
% data_static_TSTE_c_40_lfd = zeros(1,8); 


% load 'TckRstct_data_static_TSTE_c_40.mat'
% TckResults_multiCorr = TckResultCT_multiCorr; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [4];
% eleList = [58.69];
% data_mp = zeros(1,8); 


% % load 'TckRstct_austinTower_3_50_30.mat'
% TckResults_multiCorr = TckResultCT ; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [23];
% eleList = [52];
% austinTower_3_50_30_labeled = zeros(1,8); 


% % load 'TckRstct_austinTower_3_50_30.mat'
% TckResults_multiCorr = TckResultCT ; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [22];%14, 22, 31, 32];
% eleList = [31.55];%67.96, 31.55, 46.61, 52.87];
% data_20190227_01_00PM_labeled = zeros(1,8);


% % load 'TckRstct_austinTower_3_50_30.mat'
% TckResults_multiCorr = TckResultCT ; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [14];%14, 22, 31, 32];
% eleList = [64.45];%67.96, 31.55, 46.61, 52.87];
% data_mtcgps_12_12_labeled = zeros(1,8); 

% % load 'TckRstct_austinTower_3_50_30.mat'
% TckResults_multiCorr = TckResultCT ; 
% % prnList = [5, 13, 17, 19];
% % eleList = [58.69, 32.0, 16.82, 43.12];
% prnList = [14];%14, 22, 31, 32];
% eleList = [67.99];%67.96, 31.55, 46.61, 52.87];
% data_mtcgps_12_30_labeled = zeros(1,8); 

% load 'TckRstct_austinTower_3_50_30.mat'


%D_20190421_MegBOx_GPS_IQ_1550
% prnList = [1 7 9 11 22]; 
% eleList = [79.23 63.57 6.92 65.07 27.31];


% D_20190422_TSTE2_GPS_IQ_1330
% prnList = [7 8 9 18 27]; 
% eleList = [36.5 71.63 34.77 42.42 44.2];

% D_20190422_TSTE2_GPS_IQ_1335
% prnList = [1 7 8 9 18 23 27]; 
% eleList = [25.02 38.6 69.725 34.235 44.575 34.075 42.065];

% D_20190422_TSTE2_GPS_IQ_1335


% D_20190422_TSTE2_GPS_IQ_1335
% prnList = [11 16 23]; 
% eleList = [33.41 33.1 44.83];


% D_20190424_AustinTower_1400
% prnList = [8 11 18 30]; 
% eleList = [57.55 65.54 58.67 19.48];


% D_20190424_AustinTower_1410
% prnList = [7 11 18  ]; 
% eleList = [57.33 70.53 62.44];








% load 'tckResut_TckResultVT_PNT.mat'
% TckResults_multiCorr = TckResultVT; 
% prnList = [14, 22, 26, 31, 32];
% eleList = [68.8690070180876, 29.8899167453969, 70.2660966952284, 48.2252845727556, 56.9475219430146];
% PNT_VT_labeled = zeros(1,8); 


% load 'TckResultCT_multiCorr_GPS_lband_60sec.mat'
% TckResults_multiCorr = TckResultCT; 
% prnList = [10 12 14 20 25 31 32];
% eleList = [8.36, 54.23, 24.23, 47.1, 6.39, 50.04, 32.49];
% GPS_lband_60sec_labeled = zeros(1,8); 

% load 'TckResultCT_multiCorr_Data_20190116_data5.mat'
% TckResults_multiCorr = TckResultCT; 
% % prnList = [3 8 14 16 22 23 26 27 31];
% % eleList = [8.36, 54.23, 24.23, 47.1, 6.39, 50.04, 32.49, 50.04, 32.49];
% prnList = [3 8 14 16 ];
% eleList = [8.36, 54.23, 24.23, 47.1 ];
% Data_20190116_data5_labeled = zeros(1,8); 


% load 'PNT_Data_John_dyn_c_VT_125s_multiCorr_exclude3_29_BigR31_PNTpaper.mat'
% clear sum
% TckResults_multiCorr = TckResultCT; 
% prnList = [14 22 26 31 32];
% eleList = [8.36, 54.23, 24.23, 47.1, 6.39];
% Data_John_dyn_c_VT_125s_multiCorr_labeled = zeros(1,8); 

%% Initial settings (by LZD)
data_labeled = zeros(1,8); 
TckResults_multiCorr = TckResultVT ;  
prnList = [1 7 8 9 11 18]; 
eleList = [31.54 45 64.42 32.25 56.6 51.04 ];

%% Parameters (by LZD)

pdi = 1;
corrSpacing = 0.05;
maxDelay = 0.6;
numDelay = 2*maxDelay/corrSpacing + 1;
cenDelay = maxDelay/corrSpacing + 1;
spacing = -0.6:corrSpacing:0.6;
ind_start = 1; % start point, ms
numOverLap = 100; % ms

%%
for ind = 1: length(prnList)    
    prn = prnList(ind);     
    ele = eleList(ind); 
    
    a =[4092.9779845217          340.423503277404         -2.99026922880033        0.0251763660254827];  % 
    
    expectedCorr = a(1) + a(2)*ele^1 + a(3)*ele^2 + a(4)*ele^3;
    
    
%     expectedCN0 = a2(1) + a2(2)*ele^1 + a2(3)*ele^2 + a2(4)*ele^3;
%     expectedCorr = 1.8*10^(-5)*ele^4-0.0064*ele^3+0.73*ele^2+48*ele+1.2e3; %%

    dataLength = length(TckResults_multiCorr(prn).P_i); % in unit of ms
    
    % calculate the correlation value
    for i= 1:dataLength
        corr(i,:) = [
            sqrt((TckResults_multiCorr(prn).E_i_060(i))^2 + (TckResults_multiCorr(prn).E_q_060(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_055(i))^2 + (TckResults_multiCorr(prn).E_q_055(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i(i))^2 + (TckResults_multiCorr(prn).E_q(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_045(i))^2 + (TckResults_multiCorr(prn).E_q_045(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_040(i))^2 + (TckResults_multiCorr(prn).E_q_040(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_035(i))^2 + (TckResults_multiCorr(prn).E_q_035(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_030(i))^2 + (TckResults_multiCorr(prn).E_q_030(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_025(i))^2 + (TckResults_multiCorr(prn).E_q_025(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_020(i))^2 + (TckResults_multiCorr(prn).E_q_020(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_015(i))^2 + (TckResults_multiCorr(prn).E_q_015(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_010(i))^2 + (TckResults_multiCorr(prn).E_q_010(i))^2) ...
            sqrt((TckResults_multiCorr(prn).E_i_005(i))^2 + (TckResults_multiCorr(prn).E_q_005(i))^2) ...
            sqrt((TckResults_multiCorr(prn).P_i(i))^2 + (TckResults_multiCorr(prn).P_q(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i005(i))^2 + (TckResults_multiCorr(prn).L_q005(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i010(i))^2 + (TckResults_multiCorr(prn).L_q010(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i015(i))^2 + (TckResults_multiCorr(prn).L_q015(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i020(i))^2 + (TckResults_multiCorr(prn).L_q020(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i025(i))^2 + (TckResults_multiCorr(prn).L_q025(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i030(i))^2 + (TckResults_multiCorr(prn).L_q030(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i035(i))^2 + (TckResults_multiCorr(prn).L_q035(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i040(i))^2 + (TckResults_multiCorr(prn).L_q040(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i045(i))^2 + (TckResults_multiCorr(prn).L_q045(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i(i))^2 + (TckResults_multiCorr(prn).L_q(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i055(i))^2 + (TckResults_multiCorr(prn).L_q055(i))^2) ...
            sqrt((TckResults_multiCorr(prn).L_i060(i))^2 + (TckResults_multiCorr(prn).L_q060(i))^2)
            ]; 
    end
    
    codeDiscOut = TckResults_multiCorr(prn).DLLdiscri;%;    % Code disc output    codeError
    
    numMLC_1ms = 0;
    
%     meanMaxForFFT = zeros(1,300);
        
    % calculate the mean and variance of the maximum correlation value, in unit of chips
    for m=1: (dataLength-ind_start)/numOverLap 
        for n=1:numOverLap
            [maxCorr(n) maxCorrInd(m,n)] = max(corr(ind_start-1+n+(m-1)*numOverLap,:)); 
            maxCorr(n)  = corr(ind_start-1+n+(m-1)*numOverLap,13);
            tmpDelay(n) = (maxCorrInd(m,n)- mean(maxCorrInd(m,:)))*corrSpacing;%
            
            tempCodeDisc(n) = codeDiscOut(ind_start-1+n+(m-1)*numOverLap);
            
        
%             tempcorr_1ms = corr(ind_start-1+n+(m-1)*numOverLap,:);  % For 1 ms
% %             if tempcorr_1ms(1) > tempcorr_1ms(2)
% %                 numMLC_1ms = numMLC_1ms + 1;
% %             end
%             for ind_j=2:numDelay-1
%                 if tempcorr_1ms(ind_j)>tempcorr_1ms(ind_j-1) && tempcorr_1ms(ind_j)>tempcorr_1ms(ind_j+1)
%                     numMLC_1ms = numMLC_1ms + 1;
%                 end
%             end
% %             if tempcorr_1ms(end) > tempcorr_1ms(end-1)
% %                 numMLC_1ms = numMLC_1ms + 1;
% %             end
        end
%         numMLC(m) = numMLC_1ms/numOverLap;
%         numMLC_1ms = 0;
            
        
        meanMax(m) = mean(maxCorr);
        
%         meanMaxForFFT = [meanMaxForFFT(2:end), meanMax(m)];  
%         fs = 10;
%         N = 1024;  
%         mag = abs(fft(meanMaxForFFT-mean(meanMaxForFFT), N));
%         [M, I] = max(mag);
%         fftFreq(m) = (I-1)*fs/N;
%         fftMag(m) = M/(N/2);
        
        meanDelay(m) = mean((maxCorrInd(m,:)-cenDelay)*corrSpacing);
        varDelay(m) = sum(tmpDelay.^2)/numOverLap;
%         varDelay2(m) = sum((maxCorrInd(m,:)-cenDelay)*corrSpacing)/numOverLap;
        meanCodeDisc(m) = mean(tempCodeDisc);
        varCodeDisc(m) = var(tempCodeDisc);
        
        data_labeled = [data_labeled; ...
                                            [prn, ele, meanMax(m), ...
                                            meanMax(m)/expectedCorr, ...    % F1: mean maxmum corr vs. ele 
                                            -(meanDelay(m)+0.00)*1, ...               % F2: mean delay of the max corr
                                            varDelay(m), ...                % F3: var delay of the max corr
                                            meanCodeDisc(m), ...            % F4: mean code discr.
                                            varCodeDisc(m)%, ...             % F5: var code discr.
%                                             meanMax(m)/expectedCorr;    %  
%                                             numMLC(m), ...                  % F6: no. max local max
%                                             fftFreq(m), ...                 % F7: fft freq. of the max corr
%                                             fftMag(m)                       % F8: fft mag. of the max corr
                                            ]
                                         ]; 
    end
     
    
   
%    	meanMaxForFFT = zeros(1,300);
    clear   maxCorr maxCorrInd tmpDelay tempCodeDisc numMLC
    
end

%% 
prn = 9;
figure;hold on;grid on;
if flag==1
    advt = 1000;
else
    advt=0;
end
for i= 1:20
%     plot(spacing+TckResults_multiCorr(prn).codePhase(ind_start+i+69000-1 +advt),(corr(ind_start+i+69000-1 +advt,:)),'linewidth',1 ) %,'Marker','s', 'MarkerSize',8
%     plot(spacing+TckResults_multiCorr(prn).remChip(ind_start+i+1000-1 +advt),(corr(ind_start+i+1000-1 +advt,:)),'linewidth',1 ) %,'Marker','s', 'MarkerSize',8  remChip
    plot(spacing +TckResults_multiCorr(prn).remChip(ind_start+i+25000-1 +advt),(corr(ind_start+i+25000-1 +advt,:)),'linewidth',1 ) %,'Marker','s', 'MarkerSize',8  remChip
    frame = getframe(gcf);
    imind = frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if i==1
        imwrite(imind,cm,'correlation_outputs.gif','gif', 'Loopcount',inf,'DelayTime', 1e0);
    else
        imwrite(imind,cm,'correlation_outputs.gif','gif','WriteMode','append','DelayTime', 1e0);
    end

end
title(['PRN #', int2str(prn)] ,'fontsize',14);
xlabel('Time delay (Chip)','fontsize',14);
ylabel('Correlation value','fontsize',14);
set(gca,'FontSize',14);
