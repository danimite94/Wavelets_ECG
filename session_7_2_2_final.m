% session_7_2_2
% October 2015


%waverec wrong
%why 6 approximation coefficients?
clear; close all; clc; 

load ECG
dwtmode('per','nodisp'); %periodic extension for decomposition
fs=500;
%% Filtering of the signal (copy from previous exercise or load the result from a saved variable)
% clean up the ECG signal
fs = 500; %sampling rate 1000Hz
s_len = length(data);
t = [1:s_len]/fs;

figure('name','ECG Original Signal')
plot(t,data)
xlabel('Time (s)')
ylabel('Amplitude')
title('ECG Original Signal')

fft_data=fft(data);

figure('name','Frequency domain')
plot(linspace(-250,250,length(fft_data)),fftshift(fft_data))
xlabel('Frequency (HZ)')
ylabel('Amplitude')
title('ECG Original Signal')


[B,A] = butter(4,25/250,'low');

BB=[1 -1]/(1-exp(-2j*pi*0.5));
AA=[1 -0.99]/(1-0.99*exp(-2j*pi*0.5));

wo1=2*pi*50/1000;

polinome1_1=[1,-(cos(wo1)+sin(wo1)*1j)];
polinome1_2=[1,-(cos(wo1)-sin(wo1)*1j)];

polinome1_1_and_2=conv(polinome1_1, polinome1_2);

%150 Hz

wo2=2*pi*150/1000;

polinome2_1=[1,-(cos(wo2)+sin(wo2)*1j)];
polinome2_2=[1,-(cos(wo2)-sin(wo2)*1j)];

polinome2_1_and_2=conv(polinome2_1, polinome2_2);

%250 Hz

wo3=2*pi*250/1000;

polinome3_1=[1,-(cos(wo3)+sin(wo3)*1j)];
polinome3_2=[1,-(cos(wo3)-sin(wo3)*1j)];

polinome3_1_and_2=conv(polinome3_1, polinome3_2);

%350 Hz

wo4=2*pi*350/1000;

polinome4_1=[1,-(cos(wo4)+sin(wo4)*1j)];
polinome4_2=[1,-(cos(wo4)-sin(wo4)*1j)];

polinome4_1_and_2=conv(polinome4_1, polinome4_2);

%450 Hz

wo5=2*pi*450/1000;

polinome5_1=[1,-(cos(wo5)+sin(wo5)*1j)];
polinome5_2=[1,-(cos(wo5)-sin(wo5)*1j)];

polinome5_1_and_2=conv(polinome5_1, polinome5_2);

polinome_1_and_2=conv(polinome1_1_and_2,polinome2_1_and_2);
polinome_1_and_2=conv(polinome_1_and_2,polinome3_1_and_2);
polinome_1_and_2=conv(polinome_1_and_2,polinome4_1_and_2);
polinome_1_and_2=conv(polinome_1_and_2,polinome5_1_and_2);

polinome_1_and_2=conv(polinome_1_and_2,[1 1]);

normalization_factor=sum(polinome_1_and_2);

BBB=polinome_1_and_2/normalization_factor;
AAA=[1];

filters_together_num=conv(B,BB);
filters_togeter_denom=conv(A,AA);

filters_together_num=conv(filters_together_num,BBB);
filters_togeter_denom=conv(filters_togeter_denom,AAA);

data_filtered=filter(filters_together_num,filters_togeter_denom,data);

figure('name','Filtering with the three filters in cascade')
title('Plot of the original signal and the filtered one with Cascade Filter on time and frequency domain')

subplot(2,1,1)
plot(t,data)
title('Original Signal - Time Domain')
xlabel('Time')
ylabel('ECG without filtering')
ylim([-600,600])

subplot(2,1,2)
plot(t,data_filtered);
title('Filtered Signal - Time Domain')
xlabel('Time')
ylabel('ECG with the Filters')
ylim([-600,600])
%% Wavelet decomposition
% decompose the signal using Daubechies wavelet 4 (db4)
% use functions wavedec, waverec, appcoef, detcoef

% ..
n=5;

[c,l] = wavedec(data_filtered,n,'db4');

app_coef = appcoef(c,l,'db4',n);

%detail coeeficients{5,4,3,2,1}

det_coef_5 = detcoef(c,l,5);
det_coef_4 = detcoef(c,l,4);
det_coef_3 = detcoef(c,l,3);
det_coef_2 = detcoef(c,l,2);
det_coef_1 = detcoef(c,l,1);

figure('name', 'Approximation and detail coefficients');

subplot(length(l)-1,1,1)
stem(1:l(1),app_coef,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Approximation coefficient 5');

subplot(length(l)-1,1,2)
stem(l(2):l(3),det_coef_5,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Detail coefficient 5');

subplot(length(l)-1,1,3)
stem(l(3)+1:l(4),det_coef_4,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Detail coefficient 4');

subplot(length(l)-1,1,4)
stem(l(4)+1:l(5),det_coef_3,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Detail coefficient 3');

subplot(length(l)-1,1,5)
stem(l(5)+1:l(6),det_coef_2,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Detail coefficient 2');

subplot(length(l)-1,1,6)
stem(l(6)+1:l(7),det_coef_1,'r*');
xlabel('Samples');
ylabel('Amplitude');
title('Detail coefficient 1');

%clearer graphs. how?

%% reconstruction of the signal
indexes=1:length(c);
temp=[abs(c)',indexes'];
prob_sort = sortrows(temp,-1);

M=10;
rest=[zeros(length(c)-M,1),prob_sort(M+1:end,2)];
new_c=[c(prob_sort(1:M,2))',prob_sort(1:M,2)];

n_c=[new_c;rest];

ordered_new_c=sortrows(n_c,2);
recon_sig10 = waverec(ordered_new_c,l,'db4');

figure('name','Kept coefficients with the largest absolute amplitudes');

subplot 411
stem(new_c(:,2),new_c(:,1)); axis tight
xlabel('Samples');
ylabel('Amplitude');
title('10 wavelet coefficients plot');

% Compression for M={25,50,100}
M=25;
rest=[zeros(length(c)-M,1),prob_sort(M+1:end,2)];
new_c=[c(prob_sort(1:M,2))',prob_sort(1:M,2)];

n_c=[new_c;rest];

ordered_new_c=sortrows(n_c,2);
recon_sig25 = waverec(ordered_new_c,l,'db4');

subplot 412
stem(new_c(:,2),new_c(:,1)); axis tight
xlabel('Samples');
ylabel('Amplitude');
title('25 wavelet coefficients plot');

M=50;
rest=[zeros(length(c)-M,1),prob_sort(M+1:end,2)];
new_c=[c(prob_sort(1:M,2))',prob_sort(1:M,2)];

n_c=[new_c;rest];

ordered_new_c=sortrows(n_c,2);
recon_sig50 = waverec(ordered_new_c,l,'db4');

subplot 413
stem(new_c(:,2),new_c(:,1)); axis tight
xlabel('Samples');
ylabel('Amplitude');
title('50 wavelet coefficients plot');

M=100;
rest=[zeros(length(c)-M,1),prob_sort(M+1:end,2)];
new_c=[c(prob_sort(1:M,2))',prob_sort(1:M,2)];

n_c=[new_c;rest];

ordered_new_c=sortrows(n_c,2);
recon_sig100 = waverec(ordered_new_c,l,'db4');

subplot 414
stem(new_c(:,2),new_c(:,1)); axis tight
xlabel('Samples');
ylabel('Amplitude');
title('100 wavelet coefficients plot');

%% Reconstructed signals plot

N1=1;
N2=10;

new_data=resample(data,N1,N2);
new_recon_sig10=resample(recon_sig10,N1,N2);
new_recon_sig25=resample(recon_sig25,N1,N2);
new_recon_sig50=resample(recon_sig50,N1,N2);
new_recon_sig100=resample(recon_sig100,N1,N2);

figure('name', 'Comparison between original ECG and the reconstructed signal for M=10');
plot(500:700,new_data(500:700),'b-');axis tight
hold on;
plot(499:699,new_recon_sig10(500:700),'r-'); axis tight
legend('Original ECG data','Signal reconstruction using 10 wavelet coefficients');

figure('name', 'Comparison between original ECG and the reconstructed signal for M=25');
stem(1:length(new_data),new_data);axis tight
hold on;
stem(1:length(new_recon_sig25),new_recon_sig25,'r*'); axis tight
legend('Original ECG data','Signal reconstruction using 25 wavelet coefficients');

figure('name', 'Comparison between original ECG and the reconstructed signal for M=50');
stem(1:length(new_data),new_data);axis tight
hold on;
stem(1:length(new_recon_sig50),new_recon_sig50,'r*'); axis tight
legend('Original ECG data','Signal reconstruction using 50 wavelet coefficients');

figure('name', 'Comparison between original ECG and the reconstructed signal for M=100');
stem(1:length(new_data),new_data);axis tight
hold on;
stem(1:length(new_recon_sig100),new_recon_sig100,'r*'); axis tight
legend('Original ECG data','Signal reconstruction using 100 wavelet coefficients');
%%
err10 = norm(data-recon_sig10);
err25 = norm(data-recon_sig25);
err50 = norm(data-recon_sig50);
err100 = norm(data-recon_sig100);
errors= [err10,err25,err50,err100];

label={'Normalized error'};

tabling=table(err10,err25,err50,err100,'VariableNames',{'err10','err25','err50','err100'}, 'RowNames',label)
