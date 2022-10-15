%%Use case 1: ECG
close all
clear 

load ('data.mat');
fs=double(fs);
Dim_dat=length(data1);
tMax=Dim_dat/fs;
figure
subplot(2,1,1);
plot(time,data1);
xlabel('time[sec]');
ylabel('amplitude'); 
title('Raw Data');
subplot(2,1,2);
plot(time,data2);
xlabel('time[sec]');
ylabel('amplitude'); 
title('Raw Data');

N_finestre=1;

Matrix1 = reshape(data1,[Dim_dat/(tMax/N_finestre),tMax/N_finestre]);
data_dtr1=detrend(Matrix1);
data_dtr1=data_dtr1(:);
[values1,MaxInd1] = find_ecg_peak(data_dtr1);

RR_Distance1(length(MaxInd1)-1)=0;
for i=1:length(MaxInd1)-1
    RR_Distance1(i)= (MaxInd1(i+1)- MaxInd1(i))/fs;
end
RR_min1=mean(RR_Distance1)-std(RR_Distance1);
RR_medio1=mean(RR_Distance1);


Matrix2 = reshape(data2,[Dim_dat/(tMax/N_finestre),tMax/N_finestre]);
data_dtr2=detrend(Matrix2);
data_dtr2=data_dtr2(:);
[values2,MaxInd2] = find_ecg_peak(data_dtr2);

RR_Distance2(length(MaxInd2)-1)=0;
for i=1:length(MaxInd2)-1
    RR_Distance2(i)= (MaxInd2(i+1)- MaxInd2(i))/fs;
end
RR_min2=mean(RR_Distance2)-std(RR_Distance2);
RR_medio2=mean(RR_Distance2);


figure
subplot(2,1,1);
plot(time,data_dtr1);
xlabel('time[sec]');
ylabel('amplitude'); 
title('Detrended Data');
hold on
scatter((MaxInd1/fs),values1,'filled','r');
hold off

subplot(2,1,2);
plot(time,data_dtr2);
xlabel('time[sec]');
ylabel('amplitude'); 
title('Detrended Data');
hold on
scatter((MaxInd2/fs),values2,'filled','r');
hold off

bpm1=compute_heart_rate(length(MaxInd1),tMax);


bpm2=compute_heart_rate(length(MaxInd2),tMax);

M=1000;
Matrix1 = reshape(data_dtr1,[M,(Dim_dat/M)]);
Matrix1=Matrix1.';



