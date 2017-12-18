% Orlando S. Hoilett
% Unit 7 Lab Projecct
% Thursday, December 7, 2017
% Question 3: Oversampling A/D conversion (Sigma-Delta A/D) Effects
%             understood through dithering applied of visual images

clc;
clear;


%% Read Serial Data
s = serial('COM4', 'BaudRate', 115200); %choose serial port
samples = 1000;
fopen(s);

ecg = zeros(1,4);
t=0:5:5000;
for n = 1:1:samples+1
    x = fscanf(s);
    %x = str2num(x);
    if(~isnan(str2double(x)))
        ecg(n) = str2double(x);
    end
    %arr(n) = str2double(x);
%      arr(n) = fscanf(s,'%d');
end

%figure(1)
%plot(t(11:end),ecg(10:end));

%!!!!Remember to close serial port before trying to run the code again
fclose(s)



%%
Fs = 200;
ecg = ecg(2:end);
%samples = length(ecg);
t = t(1:(end-1));


%%
ecg_fft = abs(fft(ecg));
ecg_psd = (ecg_fft).^2;
ecg_freq_vec = 0:(Fs/samples):Fs-(Fs/samples);


figure(1);
%subplot(311);
subplot(211);
plot(t, ecg);
xlabel('Time(ms)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Filtered and Amplified ECG Signal','FontSize',14);

%subplot(312);
subplot(212);
%plot(ecg_freq_vec, 10*log10(ecg_fft));
plot(ecg_freq_vec, ecg_fft);
xlabel('Frequency(Hz)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('FFT of Filtered and Amplified ECG Signal','FontSize',14);
xlim([1 100]);

% subplot(313);
% %plot(ecg_freq_vec, 10*log10(ecg_psd));
% plot(ecg_freq_vec, ecg_psd);
% xlabel('Frequency(Hz)','FontSize',18);
% ylabel('Amplitude','FontSize',18);
% title('PSD of Raw ECG Signal','FontSize',18);
% xlim([1 100]);

print -dtiff ECGFFT14font


%% Design simple FIR (MA and/or DIFF) filters
Acoefs=1; % set to 1 for FIR filters

% %% 3-pt Hanning filter 
%Bcoefs=[.25 .5 .25];  % 3-pt Hanning

% %% 8-pt MA filter 
N=18;   %length moving-average filter
Bcoefs=[1/N*ones(1,N)];

% %% 2-pt difference (derivative) filter 
% Bcoefs=[1 -1];  % 

% %% 3pt central difference  (2-pt MA+2-ptdifference) filter 
% Bcoefs=[1 0 -1];  % 

%%% IIR filter needs both A and B coefficients
%% 2-pt difference (derivative) filter with extra POLE to sharpen
% Bcoefs= 
% Acoefs=


%% Z-plane Analysis of Filter
figure(12); clf
zplane(Bcoefs,Acoefs)

% use tf2zp to get Poles /Zeros
[B,A]=eqtflength(Bcoefs,Acoefs);  % tf2zp needs num&den same length (append zeros where needed)
[B,A]=eqtflength(Bcoefs,Acoefs);

%% Frequency domain response of Filter
figure(13); clf
freqz(Bcoefs,Acoefs,512,'half',Fs)   % Use help freqz to see what all the options mean.

%% Pass signal through Filter and plot result
output = filter(Bcoefs,Acoefs,ecg);  % filter input signal into output

figure(14); clf
plot(t,ecg)
hold on 
plot(t,output,'r','LineWidth',1)
hold off
axis tight;
xlabel('Time in seconds');
ylabel('Signal (blue); Filtered Signal (red)');
xlim([500 5000]);


%%
max_ecg = max(ecg);
n = 1;
beats = 0;
%loc = zeros;
HR_temp = 0;
left_index = 0;
right_index = 0;
buffer = 0;
while(n < length(ecg))
    loc = find(ecg(n:end) > max_ecg*.95,1);
    if(~isempty(loc))
        beats = beats + 1;
        if(beats == 2)
            HR_temp = 60000/(loc*5+300);
            left_index = round((loc+n)-((loc+60)*.375));
            right_index = round((loc+n)+((loc+60)*.625));
            %segment = ecg((loc+n)-((loc+60)*.375):(loc+n)+((loc+60)*.625));
            segment = ecg(left_index:right_index);
            buffer = segment;
        end
        if(beats > 2)
%             HR_temp = 60000/(loc*5+300);
%             left_index = round((loc+n)-((loc+60)*.375));
%             right_index = round((loc+n)+((loc+60)*.625));
            buffer = buffer + ecg(left_index:right_index);
        end
        n = n + loc + 300/5; %needs to be 300ms
    else
        n = n + 300/5;
    end
end
disp(2)
HR = beats*(60000/((1000/Fs)*length(t)));

figure(2)
subplot(211);
plot(segment);

buffer = buffer ./ beats; %average all segments

figure(3);
plot(buffer);


%% Calculating Normalized Correlation Coefficient Based
%  on Equation 3 from Instructional Material

segmentbar = mean(segment);
ecgbar = mean(ecg);
segmentSamples = length(segment);
ecgSamples = length(ecg);

n = 1:1:segmentSamples;
gamma = zeros(ecgSamples,1);

%equation for normalized correlation coefficient
for k = 1:(ecgSamples-segmentSamples)
    gamma(k) = sum((segment(n) - segmentbar).*(ecg(n+k)-ecgbar)) ./ sqrt(sum((segment(n) - segmentbar).^2).*sum((ecg(n+k)-ecgbar).^2));
end

figure(2);
subplot(212);
plot(t,gamma);
title('Normalized Correlation Coefficient of Raw ECG Data','FontSize',16);
ylabel('Cross-correlation with Template','FontSize',16);
xlabel('Time (ms)','FontSize',16);

%find where the correlation factor is greater than a specified threshold
peaks = find(gamma >= 0.9);

%print -dtiff ECGnCorr


%% Segment and Average
% lenb = length(segment);
% buffer = zeros(lenb, 1);
% num_peaks = length(peaks);
% num_segments = 1;
% 
% %"fencepost" solution so add the first segment
% %outside of the loop and begin counter at 1
% buffer(1:lenb) = buffer(1:lenb) + ecg((peaks(2)-left_index+1):(peaks(2)+right_index));
% %buffer(1:lenb) = buffer(1:lenb) + ecg((peaks(1)-left_index):(peaks(1)+right_index));
% counter = 1; %counter to keep track of how many segments are identified
% target_peaks = zeros(12,1); %keep track of index of peaks for each segment
% %target_peaks(counter) = peaks(1);
% target_peaks(counter) = peaks(2);
% 
% %for n=3:num_peaks
% for n=3:num_peaks
%     %the correlation factor of the ECG wavform wil be greater than
%     %threshold for several points along the segment, adding this
%     %conditional ensures that one QRS peak is not counted twice
%     if (peaks(n) - peaks(n-1)) > 60
%         %summing each segment of ECG that has correlation
%         %factor greater than threshold
%         buffer(1:lenb) = buffer(1:lenb) + ecg((peaks(n)-left_index):(peaks(n)+right_index));
%         counter = counter + 1;
%         target_peaks(counter) = peaks(n); %keep track of index of peaks selected
%     end
% end
% buffer = buffer ./ counter; %average all segments
% 
% figure(3);
% subplot(313);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(buffer);
% xlabel('Time from Start of Signal (ms)','FontSize',16);
% ylabel('Relative Amplitude','FontSize',16);
% title('Average of 12 ECG Cycles','FontSize',16);
% 
% 
% %% Plotting segments 3 and 5 of raw ECG signal
% %The instructional material handed out for this lab indicates that we
% %need to pick out segments 3 and 5, however, we are pretty sure (from looking
% %at the shape of the segments) that it should be segments of 3 and 4. The
% %segment 5 does not segment 5 in the raw ECG signal. We are guessing this
% %doesn't really matter in the grand scheme of things so we decided to plot
% %segment 5 nonetheless.
% %figure(3) -- same figure as average of all cycles
% subplot(311);
% plot(1:650,ecg(target_peaks(3)-220:target_peaks(3)+429));
% xlabel('Time from Start of Signal (ms)','FontSize',16);
% ylabel('Relative Amplitude','FontSize',16);
% title('ECG Cycle Number 3','FontSize',16);
% 
% subplot(312);
% plot(1:650,ecg(target_peaks(5)-220:target_peaks(5)+429));
% xlabel('Time from Start of Signal (ms)','FontSize',16);
% ylabel('Relative Amplitude','FontSize',16);
% title('ECG Cycle Number 5','FontSize',16);
% 
% print -dtiff Segments
% %saveas(gcf, 'SegmentsAll.png');
% 
% 
% %% All synchronized realizations superimposed on one figure
% figure(4)
% %set(gcf, 'Position', get(0, 'Screensize'));
% for n=1:length(target_peaks)
%    plot(1:650,ecg(target_peaks(n)-220:target_peaks(n)+429));
%    hold on
% end
% 
% xlabel('Normalized Time to Start of Each Cycle (ms)','FontSize',14);
% ylabel('Relative Amplitude','FontSize',14);
% title('Superimposed Synchronized Realizations','FontSize',14);
% 
% print -dtiff Realizations
% %saveas(gcf, 'AllRealizations.png');
% 
% 
% %% Evaluating the effect of averaging on the RMS of the ST segment
% %  Isolate each ST segment and assign to a holding matrix
% %figure(6)
% for n=1:length(target_peaks)
%    st = ecg(target_peaks(n)+110:target_peaks(n)+179);
%    mean_adj = mean(st);     %used to adjust the mean of each segment to zero
%    st = st - mean_adj;      %adjust mean to zero
%    %fprintf('\n%d',mean(st))
%    st_hold(n, 1:70) = st;   %assign st segment to holding matrix
%    %plot(st)
%    %hold on
%    if n==1
%        root_mean_sq(n) = rms(st_hold(n, 1:70));
%    else
%        av_st = mean(st_hold);   %average all current ST segments
%        root_mean_sq(n) = rms(av_st);    %calculate RMS of averaged
%    end
% end
% 
% figure(5)
% plot(1:length(target_peaks), root_mean_sq)
% %set(gcf, 'Position', get(0, 'Screensize'));
% title('Zero-mean RMS v. Sample Size','FontSize',14);
% xlabel('Sample size, M','FontSize',14);
% ylabel('RMS of zero-mean ST segments','FontSize',14);
% 
% print -dtiff RMS_plot
% 
% %saveas(gcf, 'RMSPlot.png');
