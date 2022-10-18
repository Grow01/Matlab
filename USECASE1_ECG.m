%% USE CASE 1: ECG
% plotto le due distribuzoni di dati forniti
close all
load('data.mat');
N = length(data1);  % numero campioni totale
T = length(time)/fs;  % periodo di campionamento
f1 = figure('Name', 'Heart Rate 1');
plot(time,data1);
xlabel('time[sec]');
ylabel('amplitude[]'); 
title('Raw data 1');
f2 = figure('Name', 'Heart Rate 2');
plot(time,data2);
xlabel('time[sec]');
ylabel('amplitude[]'); 
title('Raw data 2');
%% Detrend
% pulizia segnale
% definisco Nk intervalli temporali di durata t sec in cui applicare il
% detrend
t = 0.5;  % [sec]
k = t*fs; %% numero campioni in 0,5 sec 
Nk = N/k; %% numero finestre totale
newmatrix1 = reshape(data1,[k,Nk]);   % creo matrice con i dati forniti di k righe e Nk colonne
newmatrix2 = reshape(data2,[k,Nk]);   % per utilizzare il comando detrend
D1 = detrend(newmatrix1);  %% applico il detrend (lavora per colonne)
D2 = detrend(newmatrix2);
newdata1 = reshape(D1,[1,N]);   % vettore di dati ripuliti (1)
newdata2 = reshape(D2,[1,N]);   % " (2)
%% Find Peaks
% identifichiamo i picchi significativi del segnale
multiplier = 6;  % valore del moltiplicando
S1 = mean(newdata1) + multiplier*std(newdata1); %%  calcolo valore di soglia
S2 = mean(newdata2) + multiplier*std(newdata2);
% normdata2 = fitgmdist((newdata2)',1);
[pks2, index2] = findpeaks(newdata2,'MinPeakHeight',S2);  
[pks1, index1] = findpeaks(newdata1,'MinPeakHeight',S1);  %% vettore dei picchi
% plotto il risultato
f3 = figure('Name','Heart Rate Cleaned 1');
plot(time,newdata1);
xlabel('time[sec]');
ylabel('amplitude[]'); 
title('Detrended Data 1');
findpeaks(newdata1,'MinPeakHeight',S1);
f4 = figure('Name','Heart Rate Cleaned 2');
plot(time,newdata2);
xlabel('time[sec]');
ylabel('amplitude[]'); 
title('Detrended Data 2');
findpeaks(newdata2,'MinPeakHeight',S2);
%% Troviamo la frequenza cardiaca
fc = 1; % 1 battito al secondo
fc1 = (length(pks1)/double(T))*60; % battiti al minuto
fc2 = (length(pks2)/double(T))*60;
if ((length(pks1)/double(T)) >= fc)
    f5 = figure;
    plot(time,data1);
    xlabel('time[sec]');
    ylabel('amplitude[]'); 
    title('Raw data 1, physiological');
else
    f5 = figure;
    plot(time,data1);
    xlabel('time[sec]');
    ylabel('amplitude[]'); 
    title('Raw data 1, pathological');
end
if ((length(pks2)/double(T)) >= fc)
    f6 = figure;
    plot(time,data2);
    xlabel('time[sec]');
    ylabel('amplitude[]'); 
    title('Raw data 2, physiological');
else
    f6 = figure;
    plot(time,data2);
    xlabel('time[sec]');
    ylabel('amplitude[]'); 
    title('Raw data 2, pathological');
end
%% intervallo RR
RR1 = rand([1 length(index1)]);
for i=1:length(index1)
    if (i < length(index1))
       RR1(i) = time(index1(i+1)) - time(index1(i));
    else
        RR1(i) = 0;
    end
end
meanRR1 = mean(RR1);
RR2 = rand([1 length(index2)]);
for i=1:length(index2)
    if (i < length(index2))
       RR2(i) = time(index2(i + 1)) - time(index2(i));
    else
        RR2(i) = 0;
    end
end
meanRR2 = mean(RR2);
disp('intervallo RR medio 1')
disp(meanRR1)
disp('intervallo RR medio 2')
disp(meanRR2)
diffRR = meanRR2 - meanRR1;
disp('differenza medie intervalli RR')
disp(diffRR)
meanRR = (meanRR2 + meanRR1)/2;
%% QRS Complex
% funzione per il calcolo della durata del complesso QRS e per la
% visualizzazione del complesso QRS medio.
M = meanRR;   % intervallo temporale in sec finestra per QRS (0.08-0.12 sec)
Numfin1 = length(pks1);  % numero finestre totale
Numfin2 = length(pks2);
Nm = T/M; % numero campioni per ogni finestra 
matrixQRS1 = rand([Numfin1 Nm]);
for i=1:Numfin1
    for j=1:Nm
            matrixQRS1(i,j) = data1(index1(i) - Nm/2 + j);
    end
end
meanQRS1 = mean(matrixQRS1);
timeQRS = linspace((-M/2),(M/2),Nm);
f7 = figure('Name','QRS Complex 1');
plot(timeQRS, matrixQRS1, ':');
hold on
legend('Singles QRS Complexes 1');
plot(timeQRS,(meanQRS1)', 'k','DisplayName', 'Mean QRS Complex 1');
hold off  % plotto il complesso medio QRS 
% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
matrixQRS2 = rand([Numfin2 Nm]);
for i=1:Numfin2
    for j=1:Nm
            matrixQRS2(i,j) = data2(index2(i) - Nm/2 + j);
    end
end
meanQRS2 = mean(matrixQRS2);
f8 = figure('Name','QRS Complex 2');
plot(timeQRS,matrixQRS2, ':');
hold on
legend('Singles QRS Complexes 2');
plot(timeQRS,(meanQRS2)', 'k','DisplayName', 'Mean QRS Complex 2');
hold off
%% QRS Complex Duration
 % trovo intervallo tra R e S e lo stimo pari a metà della durata QRS
 % QRS = RS*2; stima durata QRS
 %% Histogram
 % istogramma
 
