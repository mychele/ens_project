classdef deleteSins < matlab.System & matlab.system.mixin.Propagates ...
        & matlab.system.mixin.CustomIcon
    % deleteSins Questo System Object è usato per rilevare componenti sinusoidali
    % in un segnale, con una procedura iterativa.
    
    properties (PositiveInteger, Nontunable)
        %Threshold
        % Specifica una soglia: se max(segnale) > threshold c'è un seno
        Threshold = 200
        
        %SampleRate
        % Frequenza di campionamento
        Fs = 16000
        
        %SamplesPerFrame
        % Campioni analizzati ad ogni iterazione
        SamplesPerFrame = 4000
        
        %Notch band
        % Banda del filtro ai 3 dB
        notch_band = 200;
    end
    
    properties (Access = private)
        % Alcune variabili che saranno inizializzate in seguito
        NFFT = 0;
        f;
        window;
        numbSins;
        r;
    end
    
    properties (DiscreteState)
        % Stato dell'oggetto, sarà inizializzato in seguito
        State;
    end
    
    methods
        % Costruttore
        function obj = deleteSins(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        %% Funzioni comuni
        function setupImpl(obj,~)
            % Azioni che vengono compiute una sola volta, quando l'oggetto
            % viene inizializzato
            obj.State = zeros(obj.SamplesPerFrame, 1);
            obj.NFFT = 128*2^nextpow2(obj.SamplesPerFrame); % zero packing
            obj.f = obj.Fs*linspace(0, 1, obj.NFFT); % vettore delle frequenze
            obj.window = triang(obj.SamplesPerFrame); % finestra triangolare
            % calcola i parametri del notch che non dipendono dalla
            % frequenza da eliminare
            delta = pi*obj.notch_band/(2*obj.Fs);
            obj.r=1-delta;
        end
        
        function [out, nSins, amps] = stepImpl(obj,x)
            % Funzione che viene eseguita ad ogni iterazione. Riconosce le
            % componenti sinusoidali e le filtra. Restituisce il segnale
            % filtrato, il numero di seni trovati e la loro ampiezza
            
            y = obj.window.*x; % applica la finestra triangolare
            % calcola la FFT del segnale originale e di quello finestrato
            X = fft(x, obj.NFFT);
            Y = fft(y, obj.NFFT);
            
            amps = zeros(1, 3);
%             figure
%             plot(obj.f/1000, 20*log10(abs(X)));
%             title('DFT - modulo in dB');%,'Interpreter','latex')
%             xlabel('f [kHz]');%, 'Interpreter','latex');
%             ylabel('|X| [dB]');%, 'Interpreter','latex');
%             
%             set(gcf, 'PaperUnits', 'points');
%             %set(gcf, 'PaperPosition', [0 0 1920 1080]);
%             set(gcf, 'PaperSize', [1200, 700]);
%             set(gcf, 'Color', 'w');
%             %format_ticks(gca);
%             fig=gcf;
%             set(findall(fig,'-property','FontSize'),'FontSize',33)
%             %export_fig instThr_c2_200.png -q101 -nocrop
%             
%             pause();
            
            % Calcola il massimo dello spettro del segnale
            [R, index] = max(abs(Y(1:length(Y)/2)));
            obj.numbSins = 0;
            while(R > obj.Threshold) % se supera la soglia data
                % trova la frequenza corrispondente
                fm = obj.f(index);
                obj.numbSins = obj.numbSins + 1;
                amps(obj.numbSins) = R*2/obj.NFFT;
                
                % calcola i parametri del filtro notch
                b1 = -2*cos(2*pi*fm/obj.Fs);
                b = [1 b1 1];
                a2=obj.r^2; a1 = obj.r*b1;
                a= [1 a1 a2];
                % e il filtro notch
                [N1, ~]=freqz(b, a, obj.NFFT, 'whole', obj.Fs);
                
                % filtra X, Y
                X = N1 .* X;
                Y = N1 .* Y;
                
                % ricalcola R
                [R, index] = max(abs(Y(1:length(Y)/2)));
            end
            
%             figure
%             plot(obj.f/1000, 20*log10(abs(X)));
%             title('DFT - modulo in dB del segnale filtrato');%,'Interpreter','latex')
%             xlabel('f [kHz]');%, 'Interpreter','latex');
%             ylabel('|X| [dB]');%, 'Interpreter','latex');
%             
%             set(gcf, 'PaperUnits', 'points');
%             %set(gcf, 'PaperPosition', [0 0 1920 1080]);
%             set(gcf, 'PaperSize', [1200, 700]);
%             set(gcf, 'Color', 'w');
%             %format_ticks(gca);
%             fig=gcf;
%             set(findall(fig,'-property','FontSize'),'FontSize',33)
%             %export_fig instThr_c2_200.png -q101 -nocrop
%             
%             pause();
            
            % ottieni il segnale filtrato con ifft
            out_compl = ifft(X);
            % e aggiungi il residuo dell'iterazione precedente (dovuto a
            % zero packing e finestratura)
            out = out_compl(1:length(x)) + obj.State;
            % salva l'attuale residuo nel vettore di stato
            obj.State = out_compl(length(x) + 1:2*length(x));
            % e ritorna il numero di seni trovati
            nSins = obj.numbSins;
        end
        
        function resetImpl(obj)
            % Resetta lo stato
            obj.State = zeros(obj.SamplesPerFrame, 1);
        end
        
        %% Funzioni di backup e restore, implementate di default
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.Threshold = obj.Threshold;
                s.Fs = obj.Fs;
                s.State = obj.State;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            if wasLocked
                obj.Threshold = s.Threshold;
                obj.Fs = s.Fs;
                s.State = obj.State;
            end
            loadObjectImpl@matlab.System(obj, s, wasLocked);
        end
        
        %% Simulink functions
        function z = getDiscreteStateImpl(obj)
            % Return structure of states with field names as
            % DiscreteState properties.
            z = struct([]);
        end
        
        function flag = isInputSizeLockedImpl(obj,index)
            % Set true when the input size is allowed to change while the
            % system is running.
            flag = false;
        end
        
        function sz = getOutputSizeImpl(obj)
            % Implement if input size does not match with output size.
            sz = [1 1];
        end
        
        function icon = getIconImpl(obj)
            % Define a string as the icon for the System block in Simulink.
            icon = mfilename('class');
        end
    end
    
    methods(Static, Access = protected)
        %% Simulink customization functions
        function header = getHeaderImpl(obj)
            % Define header for the System block dialog box.
            header = matlab.system.display.Header(mfilename('class'));
        end
        
        function group = getPropertyGroupsImpl(obj)
            % Define section for properties in System block dialog box.
            group = matlab.system.display.Section(mfilename('class'));
        end
    end
end
