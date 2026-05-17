classdef PAM
    %% Class PAM
    properties
        M % constellation size
        Rb % bit rate
        level_spacing % 'equally-spaced' or 'optimized'
        pulse_shape % struct containing the following fields {type: {'rectangular', 'root raised cosine', 'raised cosine'}, h: pulse shape impulse response i.e., coefficients of FIR filter,sbs: samples per symbol, and other parameters such as rolloff factor, etc}
        a % levels
        b % decision threshold
    end
    
    properties (Dependent)
        Rs % symbol rate
        optimize_level_spacing % logical variable
    end

       
    methods
        function obj = PAM(M, Rb, level_spacing, pulse_shape)
            %% Class constructor
            % Inputs
            % - M = constellation size
            % - Rb = bit rate
            % - level_spacing = 'equally-spaced' or 'optimized'
            % - puse_shape = % struct containing the following fields {type:
            % {'rectangular', 'root raised cosine', 'raised cosine'}, 
            % h: pulse shape impulse response i.e., coefficients of FIR filter,
            % Mct: oversampling ratio of pulse shapping filter, and
            % other parameters such as rolloff factor, etc}
            
            obj.M = M;
            obj.Rb = Rb;
            
            if exist('level_spacing', 'var')
                obj.level_spacing = level_spacing;
            else
                obj.level_spacing = 'equally-spaced';
            end
            
            if exist('pulse_shape', 'var')
                assert(~strcmp(class(pulse_shape), 'function_handle'), 'PAM: behaviour of PAM has changed. Now PAM expects a struct as opposed to the function handle.')
                obj.pulse_shape = pulse_shape;
                obj.pulse_shape.h = obj.norm_filter_coefficients(obj.pulse_shape.h);
            else
                obj.pulse_shape = select_pulse_shape('rect', 1);
            end
            
            obj = obj.reset_levels();
        end
    end
    
    methods 
        %% Get and set methods
        function Rs = get.Rs(self)
            %% Symbol-rate assuming rectangular pulse
            Rs = self.Rb/log2(self.M);
        end
        
        function optimize_level_spacing = get.optimize_level_spacing(self)
            %% True if level_spacing == 'optimized'
            optimize_level_spacing = strcmp(self.level_spacing, 'optimized');
        end
        
        function H = Hpshape(self, f)
            %% Frequency response of PAM pulse shape
            fs = self.Rs*self.pulse_shape.sps;
            delay = grpdelay(self.pulse_shape.h, 1, 1);
            H = freqz(self.pulse_shape.h/abs(sum(self.pulse_shape.h)), 1, f, fs)...
                .*exp(1j*2*pi*f/fs.*delay); % remove group delay
        end        
    end
       
    methods
        %% Levels and decision thresholds     
        function self = set_levels(self, levels, thresholds)
            %% Set levels to desired values
            % Levels and decision thresholds are normalized that last level is unit
            assert(length(levels) == self.M, 'mpam/set_levels: invalid number of levels');
            assert(length(thresholds) == self.M-1, 'mpam/set_levels: invalid number of decision thresholds');
            if size(levels, 1) < size(levels, 2) % ensure levels are M x 1 vector
                levels = levels.'; 
            end
            if size(thresholds, 1) < size(thresholds, 2)  % ensure thresholds are M x 1 vector
                thresholds = thresholds.';
            end
                
            self.a = levels;
            self.b = thresholds;
        end
        
        function self = norm_levels(self)
            %% Normalize levels so that last level is 1
            self.b = self.b/self.a(end);
            self.a = self.a/self.a(end);
        end
        
        function self = unbias(self)
            %% Remove DC bias from levels and normalize to have excusion from -1 to 1
            self.b = self.b - mean(self.a);
            self.a = self.a - mean(self.a);
            self = self.norm_levels;
        end
            
        function self = reset_levels(self)
            %% Reset levels and decision thresholds to original configuration
            self.a = ((0:2:2*(self.M-1))/(2*(self.M-1))).';
            self.b = ((1:2:(2*(self.M-1)-1))/(2*(self.M-1))).';
        end

    end

    %% Auxiliary functions
    methods
        function h = norm_filter_coefficients(~, h)
            %% Normalize coefficients of FIR filter h so that impulse response of h at t = 0 is 1
            n = length(h);
            if mod(n, 2) == 0 % even
                h = 2*h/(h(n/2) + h(n/2+1));
            else
                h = h/h((n+1)/2);
            end
        end
    end
end