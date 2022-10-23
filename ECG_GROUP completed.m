% USE CASE 1: ECG
% plotto le due distribuzoni di dati forniti
close all
load('data.mat');
N=length(data1); %trovo quanti dati sono
figure % apro unna finestra grafica
tmax=N/fs; %trovo il tempo di durata della ECG

% plotto i dati raw

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
nrfin=tmax/k;
newdata1=lineardetrend(data1,N,nrfin);
newdata2=lineardetrend(data2,N,nrfin);

%plotto i dati detrendati 

figure
subplot(2,1,1)
plot(time,newdata1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data1");

subplot(2,1,2)
plot(time,newdata2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data2");

%% Find Peaks
% identifichiamo i picchi significativi del segnale

m = 6;  % valore del moltiplicando
S=find_thr(newdata1,m); %%  calcolo valore di soglia

[values2, Maxindex2] = find_ecg_peak(newdata2,S);  
[values1, Maxindex1] = find_ecg_peak(newdata1,S);  %% vettore dei picchi

% plotto il risultato

figure;
subplot(2,1,1)
plot(time,newdata1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData");
hold on;
scatter(Maxindex1/double(fs),values1,'filled');


subplot(2,1,2)
plot(time,newdata2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData");
hold on;
scatter(Maxindex2/double(fs),values2,'filled');

%% Troviamo la frequenza cardiaca
bpm1=compute_heart_rate(Maxindex1,tmax);
bpm2=compute_heart_rate(Maxindex2,tmax);

%% intervallo RR
RR1 = rand([1 length(Maxindex1)]);
for i=1:length(Maxindex1)
    if (i < length(Maxindex1))
       RR1(i) = time(Maxindex1(i+1)) - time(Maxindex1(i));
    else
        RR1(i) = 0;
    end
end
meanRR1 = mean(RR1);
RR2 = rand([1 length(Maxindex2)]);
for i=1:length(Maxindex2)
    if (i < length(Maxindex2))
       RR2(i) = time(Maxindex2(i + 1)) - time(Maxindex2(i));
    else
        RR2(i) = 0;
    end
end
meanRR2 = mean(RR2);
disp('intervallo RR medio 1 [sec]')
disp(meanRR1)
disp('intervallo RR medio 2 [sec]')
disp(meanRR2)
diffRR = meanRR2 - meanRR1;
disp('differenza medie intervalli RR [sec]')
disp(diffRR)
meanRR = (meanRR2 + meanRR1)/2;
%% QRS Complex
% funzione per il calcolo della durata del complesso QRS e per la
% visualizzazione del complesso QRS medio.
M1 = meanRR1;   % intervallo temporale in sec finestra 
Numfin1 = length(values1);  % numero finestre totale
Numfin2 = length(values2);
Nm1 = tmax/M1; % numero campioni per ogni finestra 
matrixQRS1 = rand([Numfin1 Nm1]);
grayColor = [.7 .7 .7];
for i=1:Numfin1
    for j=1:Nm1
            matrixQRS1(i,j) = data1(Maxindex1(i) - Nm1/2 + j);
    end
end
meanQRS1 = mean(matrixQRS1);
timeQRS1 = linspace((-M1/2),(M1/2),Nm1);
figure('Name','QRS Complex 1');
plot(timeQRS1, matrixQRS1, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 1');
plot(timeQRS1,(meanQRS1)', 'k','DisplayName', 'Mean QRS Complex 1', 'LineWidth', 3);
hold off  % plotto il complesso medio QRS 
% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
M2 = meanRR2;
Nm2 = tmax/M2;
matrixQRS2 = rand([Numfin2 Nm2]);
for i=1:Numfin2
    for j=1:Nm2
            matrixQRS2(i,j) = data2(Maxindex2(i) - Nm2/2 + j);
    end
end
meanQRS2 = mean(matrixQRS2);
timeQRS2 = linspace((-M2/2),(M2/2),Nm2);
figure('Name','QRS Complex 2');
plot(timeQRS2,matrixQRS2, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 2');
plot(timeQRS2,(meanQRS2)', 'k','DisplayName', 'Mean QRS Complex 2', 'LineWidth', 3);
hold off
%% QRS Complex Duration
% trovo i picchi con findpeaks (senza ribaltare) per ogni finestra di campioni
 % calcolo la durata temporale tra il secondo picco prima e dopo rispetto al picco R e
 % lo stimo circa pari alla QRS duration
% ne prendo due per confronto con il QRS medio
durQRS1 = zeros([1,Numfin1]);
durQRS2 = zeros([1,Numfin2]);

% detrend
newmatrixQRS1 = zeros([Numfin1,Nm1]);
newmatrixQRS2 = zeros([Numfin2,Nm2]);
k1 = 5; % campioni per finestra detrend
k2 = 3;

for i=1:Numfin1
    Nk = round(numel(matrixQRS1(i,:))/k1);
    D1 = detrend(reshape(matrixQRS1(i,:),Nk,[]));
    newmatrixQRS1(i,:) = reshape(D1,[1,numel(matrixQRS1(i,:))]);
end
for i=1:Numfin2
    Nk = round(numel(matrixQRS2(i,:))/k2);
    D2 = detrend(reshape(matrixQRS2(i,:),Nk,[]));
    newmatrixQRS2(i,:) = reshape(D2,[1,numel(matrixQRS2(i,:))]);
end

% abs

ABS1 = zeros([Numfin1,Nm1]);
ABS2 = zeros([Numfin2,Nm2]);

for i=1:Numfin1
    ABS1(i,:) = abs(newmatrixQRS1(i,:));
end
for i=1:Numfin2
    ABS2(i,:) = abs(newmatrixQRS2(i,:));
end

% findpeaks (selezione solo picchi R e minimi S)

sogliaQRS = 0.013;
K1 = 1.73;   % rapporto tra QRS medio e RS medio nel primo e secondo caso ottenuto da grafici QRS medi 
K2 = 1.82;   % dalla durata RS posso stimare quella del complesso totale 

for i=1:Numfin1
    [picchi,indici] = findpeaks(ABS1(i,:), 'MinPeakHeight',sogliaQRS,'NPeaks',2,'MinPeakDistance',8);
    durQRS1(1,i) = (timeQRS1(indici(2)) - timeQRS1(indici(1)))*K1;
end
for i=1:Numfin2
    [picchi,indici] = findpeaks(ABS2(i,:), 'MinPeakHeight',sogliaQRS,'NPeaks',2, 'MinPeakDistance',8);
    durQRS2(1,i) = (timeQRS2(indici(2)) - timeQRS2(indici(1)))*K2;
end
% volendo si può aggiungere un controllo in più con un ciclo for che prende
% solo le durate minori di 0,12 sec e ne fa la media, in modo da eliminare
% gli effetti del rumore ed eventuali errori 

% meandurQRS1 = mean(durQRS1);
% meandurQRS2 = mean(durQRS2);

sumQRS = 0;
cont = 0;
for i=1:Numfin1
    if(durQRS1(1,i) <= 0.12)
        tmt = sumQRS;
        sumQRS = tmt + durQRS1(1,i);
        contt = cont;
        cont = contt + 1;
    end
end
meandurQRS1 = sumQRS/cont;
sumQRS = 0;
cont = 0;
for i=1:Numfin2
    if(durQRS2(1,i) <= 0.12)
        tmt = sumQRS;
        sumQRS = tmt + durQRS2(1,i);
        contt = cont;
        cont = contt + 1;
    end
end
meandurQRS2 = sumQRS/cont;

disp('Mean QRS Complex Duration 1 [sec] 1')
disp(meandurQRS1)
disp('Mean QRS Complex Duration 2 [sec] 1')
disp(meandurQRS2)

 %% Histogram
 % istogramma intervallo RR
n_bins = 100; 
min_minRR = min([RR1,RR2]);
max_maxRR = max([RR1,RR2]);

bin_edgesRR = linspace(min_minRR, max_maxRR, n_bins);  

% bin_centersRR = bin_edgesRR(1:end-1) + mean(diff((bin_edgesRR))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
figure('Name','RR interval Histogram');
histogram(RR1,bin_edgesRR);
hold on
histogram(RR2, bin_edgesRR,'FaceColor','r');
legend('PDF RR interval 1', 'PDF RR interval 2');
 % istogramma QRS duration
 n_bins = 25; 
min_minQRS = min([durQRS1,durQRS2]);
max_maxQRS = max([durQRS1,durQRS2]);

bin_edgesQRS = linspace(min_minQRS, max_maxQRS, n_bins); 

% bin_centersQRS = bin_edgesQRS(1:end-1) + mean(diff((bin_edgesQRS))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
figure('Name','QRS duration Histogram');
histogram(durQRS1,bin_edgesQRS);
hold on
histogram(durQRS2, bin_edgesQRS,'FaceColor','r');
legend('PDF QRS interval 1', 'PDF QRS interval 2');
% written by Davide, Edo, Anna, Matte
