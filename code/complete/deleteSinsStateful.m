classdef deleteSinsStateful < matlab.System & matlab.system.mixin.Propagates ...
        & matlab.system.mixin.CustomIcon
    % deleteSins This System object is used to detect sinusoidal components
    % in a signal, in an iterative fashion, according to a certain
    % threshold which has to be tuned accoding to the kind of signal we are
    % considering.
    % It computes R = max(sig)/median(sig) and if R is above this threshold
    % a sin is detected and then filtered out.
    
    
    properties
        % Public, tunable properties.
    end
    
    properties (PositiveInteger, Nontunable)
        %Threshold
        % Specify a threshold such that if R > threshold there is a sin
        Threshold = 200
        
        %SampleRate
        % Fs of the signal
        Fs = 16000
        
        %SamplesPerFrame
        % Samples in a frame
        SamplesPerFrame = 4000
    end
    
    properties (Nontunable)
        %StopBand
        % StopBand of the notch filters
        StopBand = 1
    end
    
    
    properties (Access = private)
        NFFTx = 0;
        NFFTy = 0;
        fx;
        fy;
        window;
        numbSins;
    end
    
    properties (DiscreteState)
        State;
        prev_samples;
    end
    
    methods
        % Constructor
        function obj = deleteSinsStateful(varargin)
            % Support name-value pair arguments when constructing the object.
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        %% Common functions
        function setupImpl(obj,~)
            % Implement tasks that need to be performed only once,
            % such as pre-computed constants.
            obj.State = zeros(obj.SamplesPerFrame, 1);
            obj.NFFTx = 4*2^nextpow2(obj.SamplesPerFrame);
            obj.NFFTy = 4*2^nextpow2(2*obj.SamplesPerFrame);
            obj.fx = obj.Fs*linspace(0, 1, obj.NFFTx);
            obj.fy = obj.Fs*linspace(0, 1, obj.NFFTy);
            obj.window = triang(2*obj.SamplesPerFrame);
            obj.prev_samples = zeros(obj.SamplesPerFrame, 1);
        end
        
        function [out, nSins] = stepImpl(obj,x)
            z = vertcat(obj.prev_samples, x);
            y = obj.window.*z;
            
            % compute the DFT
            X = fft(x, obj.NFFTx);
            Y = fft(y, obj.NFFTy);
            Y = Y - mean(abs(Y));
                        
            % Initial R
            R = max(abs(Y(1:length(Y)/2))) / median(abs(Y(1:length(Y)/2)));
            obj.numbSins = 0;
            while(R > obj.Threshold)
                [~, index] = max(Y(1:length(Y)/2));
                %compute the related frequency
                fm = obj.fy(index);
                obj.numbSins = obj.numbSins + 1;
                %disp(numbSins);
                %disp(fm);
                %compute notch filter parameters
                beta1 = cos(2*pi/obj.Fs*fm);
                alpha1 = 0.8;
                K1 = (1 + alpha1)/2;
                a1 = [K1, -2*K1*beta1, K1];
                b1 = [1, -beta1*(alpha1 + 1), alpha1];
                
                %compute notch filter
                [N1, ~]=freqz(a1, b1, obj.NFFTx, 'whole', obj.Fs);
                [N2, ~]=freqz(a1, b1, obj.NFFTy, 'whole', obj.Fs);

                
                %filter X and Y
                Y = N2 .* Y;
                X = N1 .* X;
                                
                %compute R
                R = max(abs(Y(1:length(Y)/2))) / median(abs(Y(1:length(Y)/2)));
            end
            
            obj.prev_samples = x;
            out_compl = ifft(X);
            out = out_compl(1:length(x)) + obj.State;
            obj.State = out_compl(length(x) + 1:2*length(x));
            nSins = obj.numbSins;
            
        end
         
        function resetImpl(obj)
            % Initialize discrete-state properties.
            obj.State = zeros(obj.SamplesPerFrame, 1);
        end
        
        %% Backup/restore functions
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
