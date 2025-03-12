%%%%%%%%%%多径信号仿真
clear all;clc
%mul_tao=[97 146 195 244];
%mul_tao=[97 146 195 244 293 342 390 439 488 610 732 855 977 1099 1343 1465 1587 1832 1954 2076 2320 2443 2565 2809 2931];
%mul_tao=[1343 1465 1587 1832 1954 2076 2320 2443 2565 2809 2931];
mul_tao=linspace(0,2500,26);
% taom=(1:10);
% mul_tao=[];
% mul_tao=taom*0.01*5000/1.023;
deg=3.1415926/180;
%mul_num=1;
mul_num=1;
for x=mul_tao
    for w=1:1
        multao=x;
        losphase=0/180*3.1415;
        mulphase=10/180*3.1415;
        mulamp = 0.5;
        
        settings.fs=5e9;
        fs=settings.fs;
        settings.code = 1.023e6;
        settings.fre=1575.42e6;
        
        nn=[1:fs*0.001];
        cosin = cos(2*pi*settings.fre*nn/fs+losphase);
        cosin2 = cos(2*pi*settings.fre*nn/fs+mulphase);
        CAcodesTable1 = generateCAcode(1);
        CAcodesTable2 = generateCAcode(1);
        test = mod(ceil([1:fs*0.001]/fs*(1.023e6)),1023);
        test (4995113:end) = 1023;%%%
        gold1 = CAcodesTable1(test);
        
        test2=mod(ceil([1:fs*0.001]/fs*(1.023e6)),1023);
        
        test2 (4995113:end) = 1023;%%%
        goldoffset = CAcodesTable1(test2);
        
        gold2 = CAcodesTable1(test);%modn(ceil([1:fs*0.001]/fs*(1.023e6)),1023)
        
        long = length(gold1);
        gold1 = [gold1(long/2+1:end) gold1(1:long/2)].*cosin;
        goldoffset =  mulamp*[goldoffset(long/2+multao+1:end) goldoffset(1:long/2+multao)].*cosin2;
        codefreq=conj(fft(gold2));
        noise = randn(1,long)*5;
        %load randnoise.mat
        goldsum=[];
        goldsum = gold1+goldoffset;
        goldsum=goldsum+noise;
        %     gold3=gold1+noise;
        %     for i =1:long
        %     if(gold3(i))>=0
        %         gold3(i) = 1;
        %     else
        %         gold3(i) = -1;
        %     end
    end
    for i =1:long
        if(goldsum(i))>=0
            goldsum(i) = 1;
        else
            goldsum(i) = -1;
        end
    end
    BandWidth = 2*settings.code;
    %验证是否符合带通采样 或者低通采样
%     n=round(2*settings.fi/settings.fs)+1;
%     for i=1:n
%         if ((2*settings.fi-BandWidth)/settings.fs>=i && (2*settings.fi+BandWidth)/settings.fs<=i+1 )|| (settings.fs>2*BandWidth)
%             break;
%         else
%             disp (['Sampling rate is not valid ']);
%         end
%     end
end 
    %低通滤波
    band_mul = 3;
    data_f = fft(data);
    data_f(1:round(settings.fi*settings.PIT-BandWidth*settings.PIT*band_mul))=0;
    data_f(round(settings.fi*settings.PIT+BandWidth*settings.PIT*band_mul):round(settings.fs*settings.PIT-settings.fi*settings.PIT-BandWidth*settings.PIT*band_mul))=0;
    data_f(round(settings.fs*settings.PIT-settings.fi*settings.PIT+BandWidth*settings.PIT*band_mul):settings.fs*settings.PIT)=0;
    data = real(ifft(data_f));
    clear data_f;
    %%去噪声
    %     band_mul = 20;
    %     data_f = fft(goldsum);
    %     data_f(1:(settings.fre/1000-settings.code*band_mul/1000))=0;
    %     data_f((settings.fre/1000+settings.code*band_mul/1000):(5e6-settings.fre/1000-settings.code*band_mul/1000))=0;
    %     data_f((5e6-(settings.fre/1000-settings.code*band_mul/1000)):5e6)=0;
    %     goldsum = ifft(data_f);
    %     clear data_f;
    % fc1  = settings.fre+500*(32-101)+50*(2-11);
    %%%%%%%%%%%捕获
    phase=0;
    % for i=1:21
    % fc(i) = fc1+50*(i-11);
    expfreq=exp(1i*(2*pi*settings.fre*nn/fs+phase/180*3.1415));
    result= ifft(fft(real(expfreq).*goldsum+1i*imag(expfreq).*goldsum) .* codefreq);
    
    result_I = single(real(result));
    result_Q = single(imag(result));
    acqabsresult = single(abs(result).^2);
    
    
    result_I = fliplr(single(real(result)));
    result_Q = fliplr(single(imag(result)));
    acqabsresult = fliplr(acqabsresult);
    
    [peak c2] = max(acqabsresult);
    
    figure(101)
    plot(result_I(c2-5000:c2+5000),'r')
    hold on
    plot(result_Q(c2-5000:c2+5000),'b')
    hold on
    fprintf('******************************\n');
    fprintf('%f\n',  peak/mean(acqabsresult));