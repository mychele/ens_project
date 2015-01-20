% Questo codice acquisisce, elabora riproduce un file .wav in tempo reale.

% Inizializzazione
clear all;
clc;
close all;

fileName = 'segnale_107.wav';
saveTo = 'Polese_Michele.wav';
% Numero di campioni che vengono acquisiti
FrameSize = 2000;
% Inizializza l'oggetto che importa i campioni del file .wav
AR = dsp.AudioFileReader('Filename', fileName, 'SamplesPerFrame', FrameSize);
Fs = AR.SampleRate; % frequenza di campionamento

% Inizializza l'oggetto che riproduce i campioni
AP = dsp.AudioPlayer('SampleRate', Fs);

% Inizializza l'oggetto che  riconosce le componenti sinusoidali e filtra il segnale
% I parametri sono commentati a parte nella relazione
iterativeSinsDel = deleteSins('Threshold', 120, 'Fs', Fs, 'SamplesPerFrame', FrameSize, 'notch_band', 200);

% Inizializza l'oggetto che registra il file audio su file
toFile = dsp.AudioFileWriter('Filename', saveTo, 'SampleRate', Fs);

% Vettore che registra il numero di componenti sinusoidali presenti
% contemporanemente
sinStartTimes = ones(1, Fs/FrameSize); % Fs/FrameSize � il numero massimo di iterazioni
% indice
count = 1;

while ~isDone(AR)
    % acquisizione di FrameSize campioni del segnale
    audioIn = step(AR);
    % filtraggio
    [audioOut, nSins] = step(iterativeSinsDel, audioIn);
    sinStartTimes(count) = nSins;
    % salva su file
    step(toFile, audioOut);
    % riproduzione
    step(AP, audioOut);

    count = count + 1;
end

% rilascia le risorse
release(AR);
release(AP);
release(toFile);
release(iterativeSinsDel);
