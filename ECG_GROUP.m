%% USE CASE 1: ECG
% plotto le due distribuzoni di dati forniti
clear
close all
load('data.mat');
N=length(data1); %trovo quanti dati sono
fs=double(fs);
tmax=N/fs; %trovo il tempo di durata della ECG

% plotto i dati raw
figure('Name','Data') % apro unna finestra grafica
subplot(2,1,1)
plot(time,data1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Raw Data1");

subplot(2,1,2)
plot(time ,data2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Raw Data2");


% detrendo

k =0.5; % [s]
n_fin=tmax/k;
data_dtr1=lineardetrend(data1,N,n_fin);
data_dtr2=lineardetrend(data2,N,n_fin);

%plotto i dati detrendati 

figure('Name','Detrend')
subplot(2,1,1)
plot(time,data_dtr1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data1");

subplot(2,1,2)
plot(time,data_dtr2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data2");

%% Find Peaks
% identifichiamo i picchi significativi del segnale

multiplier = 6;  % valore del moltiplicando
Soglia1=find_thr(data_dtr1,multiplier); %%  calcolo valore di soglia
Soglia2=find_thr(data_dtr2,multiplier);

[values2, Maxindex2] = find_ecg_peak(data_dtr2,Soglia2);  
[values1, Maxindex1] = find_ecg_peak(data_dtr1,Soglia1);  %% vettore dei picchi

% plotto il risultato

figure('Name','Peaks');
subplot(2,1,1)
plot(time,data_dtr1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData1");
hold on;
scatter(Maxindex1/fs,values1,'filled');


subplot(2,1,2)
plot(time,data_dtr2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData2");
hold on;
scatter(Maxindex2/fs,values2,'filled');

%% Troviamo la frequenza cardiaca
bpm1=compute_heart_rate(Maxindex1,tmax);
bpm2=compute_heart_rate(Maxindex2,tmax);

%% intervallo RR
RR1 = zeros([1 (length(Maxindex1)-1)]);
for i=1:length(Maxindex1)-1
    RR1(i) = time(Maxindex1(i+1)) - time(Maxindex1(i));
end
Media_RR1 = mean(RR1);

RR2 = zeros([1 (length(Maxindex2)-1)]);
for i=1:length(Maxindex2)-1
    RR2(i) = time(Maxindex2(i + 1)) - time(Maxindex2(i));
end
Media_RR2 = mean(RR2);

disp('intervallo RR medio 1 [sec]')
disp(Media_RR1)

disp('intervallo RR medio 2 [sec]')
disp(Media_RR2)

diff_RR = Media_RR2 - Media_RR1;
disp('differenza medie intervalli RR [sec]')
disp(diff_RR)
Media_RR = (Media_RR2 + Media_RR1)/2;
%% QRS Complex
% funzione per il calcolo della durata del complesso QRS e per la
% visualizzazione del complesso QRS medio.
M1 = Media_RR1;   % intervallo temporale in sec finestra 
Numfin1 = length(values1);  % numero finestre totale
Numfin2 = length(values2);
N_camp1 = round(tmax/M1); % numero campioni per ogni finestra 
Matrix_QRS1 = zeros([Numfin1 N_camp1]);
grayColor = [.7 .7 .7];
for i=1:Numfin1
    for j=1:N_camp1
            Matrix_QRS1(i,j) = data1(Maxindex1(i) - N_camp1/2 + j);
    end
end
Media_QRS1 = mean(Matrix_QRS1);
time_QRS1 = linspace((-M1/2),(M1/2),N_camp1);
figure
subplot(2,2,1)
plot(time_QRS1, Matrix_QRS1, ':', 'Color', grayColor);
xlabel("time[sec]");
ylabel("amplitude"); 
title("QRS Complex 1");
hold on
legend('Singles QRS Complexes 1');
plot(time_QRS1,(Media_QRS1)', 'k','DisplayName', 'Mean QRS Complex 1', 'LineWidth', 3);
hold off  % plotto il complesso medio QRS 
% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
M2 = Media_RR2;
N_camp2 = round(tmax/M2);
Matrix_QRS2 = zeros([Numfin2 N_camp2]);
for i=1:Numfin2
    for j=1:N_camp2
            Matrix_QRS2(i,j) = data2(Maxindex2(i) - N_camp2/2 + j);
    end
end
Media_QRS2 = mean(Matrix_QRS2);
time_QRS2 = linspace((-M2/2),(M2/2),N_camp2);
subplot(2,2,2)
plot(time_QRS2,Matrix_QRS2, ':', 'Color', grayColor);
xlabel("time[sec]");
ylabel("amplitude"); 
title("QRS Complex 2");
hold on
legend('Singles QRS Complexes 2');
plot(time_QRS2,(Media_QRS2)', 'k','DisplayName', 'Mean QRS Complex 2', 'LineWidth', 3);
hold off
%% QRS Complex Duration
 % trovo intervallo tra R e S e lo stimo pari a metà della durata QRS
 % QRS = RS*2; stima durata QRS
Max_QRS1 = max(Media_QRS1);
Min_QRS1 = min(Media_QRS1);
for i=1:N_camp1
    if(Media_QRS1(i) == Max_QRS1)
        Imax1 = i;
    end
end
for i=1:N_camp1
    if(Media_QRS1(i) == Min_QRS1)
        Imin1 = i;
    end
end
durQRS1 = (time_QRS1(Imin1) - time_QRS1(Imax1))*2;
disp('Mean Duration QRS Complex 1 [sec]')
disp(durQRS1)
Max_QRS2 = max(Media_QRS2);
Min_QRS2 = min(Media_QRS2);
for i=1:N_camp2
    if(Media_QRS2(i) == Max_QRS2)
        Imax2 = i;
    end
end
for i=1:N_camp2
    if(Media_QRS2(i) == Min_QRS2)
        Imin2 = i;
    end
end
durQRS2 = (time_QRS2(Imin2) - time_QRS2(Imax2))*2;
disp('Mean Duration QRS Complex 2 [sec]')
disp(durQRS2)
% single QRS durations

 %% Histogram
 % istogramma intervallo RR
n_bins = 100; % decido io il numero di bins
min_minRR = min([RR1,RR2]);
max_maxRR = max([RR1,RR2]);

bin_edgesRR = linspace(min_minRR, max_maxRR, n_bins);  % devo prendere il minimo globale tra i due e il massimo globale 
               % cerco i limiti dei vari bins per avere degli intervalli in
               % cui inserire i vari campioni X i-esimi

bin_centersRR = bin_edgesRR(1:end-1) + mean(diff((bin_edgesRR))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
figure('Name','RR interval Histogram');
histogram(RR1,bin_edgesRR);
hold on
histogram(RR2, bin_edgesRR,'FaceColor','r');
% istogramma QRS duration
% min_minQRS = min([x1,y2]);
% max_maxQRS = max([x1,y2]);
% 
% bin_edgesQRS = linspace(min_minQRS, max_maxQRS, n_bins);  % devo prendere il minimo globale tra i due e il massimo globale 
%                % cerco i limiti dei vari bins per avere degli intervalli in
%                % cui inserire i vari campioni X i-esimi
% 
% bin_centersQRS = bin_edgesQRS(1:end-1) + mean(diff((bin_edgesQRS))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
% figure('Name','QRS interval Histogram');
% histogram(x1,bin_edgesQRS);
% hold on
% histogram(y2, bin_edgesQRS,'FaceColor','r');