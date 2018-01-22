% session_7_2_1
% October 2015

clear; close all; clc; 

load ECG.mat

% clean up the ECG signal
fs = 500; %sampling rate 1000Hz
s_len = length(data);
t = [1:s_len]/fs;

figure('name','ECG Original Signal')
plot(t,data)  %data contains the ECG data, single channel, sampled at 500Hz
xlabel('Time (s)')
ylabel('Amplitude')
title('ECG Original Signal')


% Hanning Filter
% H(z) = (1/4) * ( 1 + 2z^(-1) + z^(-2) )

B=[0.25 0.5 0.25];
A=[1 0 0];

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
title('Plot of the original signal and the filtered one ')

plot(t,data,'r-')
xlabel('Time')
ylabel('ECG (amplitude)')
ylim([-200,200])
hold on;
plot(t,data_filtered);
ylim([-200,200])

legend('Original Signal - Time Domain', 'Filtered Signal - Time Domain');
% ..

%% Adapt the code from session_7_2_0_demo to compress the clean ECG signal using DCT
% compress using 70, 90 and 99% of the energy. Compute the percentage of
% coefficients that are kept and the associated RMSE of the reconstructions
X = dct(data_filtered);
[XX, ind] = sort(X.^2,'descend'); %squaring for energy measure

totalE = sum(XX);
partialE = cumsum(XX);

%perc = [0.7, 0.9,0.99]; %percentage of energy to be kept
rmse=[];
featfrac=[];
perc=0.8;
    %find the index where the remaining energy fraction is greater than perc
    numcoef = find(partialE/totalE >= perc,1); 

    %set lowest coefficients to zero and reconstruct the signal
    X(ind(numcoef+1:end)) = 0;
    y2 = idct(X);
    rmse = [rmse; sqrt(mean((data_filtered-y2).^2))];
    featfrac = [featfrac; round(numcoef/length(data_filtered)*1000)/10]; %in new matlab versions, you can use 100*round(numcoef/length(y),3)

    % plot results
    figure;

    subplot 311
    plot(data_filtered); axis tight
    xlabel('Time (s)');
    title('Original signal');

    subplot 312
    stem(X); axis tight
    title([num2str(numcoef) ' DCT coefficients (' num2str(featfrac) '%) kept for E/E_{tot} = ' num2str(i)]);

    subplot 313
    plot(data_filtered); hold on;
    plot(y2); axis tight
    xlabel('Time (s)');
    title(['Reconstruction: rmse = ' num2str(rmse)]);


labels={'RMSE (compressed signal={70%,90%,99%} of the signal energy)','Fraction of coefficients'};
table=table(rmse, featfrac,'RowNames',labels)
