classdef laser 
    %% Class laser defining semiconductor laser characteristics and source impairments
    properties
        lambda 
        PdBm   
        RIN    
        linewidth 
        freqOffset = 0;  
        alpha = 0; 
        H = @(f) ones(size(f)); 
    end
    
    properties(Dependent)
        PW 
        wavelength 
    end      
    
    methods
        function obj = laser(lambda, PdBm, RIN, linewidth, freqOffset)
            %% Class constructor
            obj.lambda = lambda;
            obj.PdBm = PdBm;
            
            if exist('RIN', 'var') 
                obj.RIN = RIN;
            end
            
            if exist('linewidth', 'var')
                obj.linewidth = linewidth;
            end
            
            if exist('freqOffset', 'var')
                obj.freqOffset = freqOffset;
            end         
        end
        
        function LaserTable = summary(self)
            %% Generate a summary table of the laser class operational parameters
            disp('-- Laser class parameters summary:')
            rows = {'Wavelength'; sprintf('Maximum frequency offset from %.2f nm', self.lambda*1e9);...
                'Power'; 'Relative intensity noise'; 'Linewidth'};
            Parameters = {'lambda'; 'freqOffset'; 'PdBm'; 'RIN'; 'linewidth'};
            Values = [self.lambda*1e9; max(self.freqOffset)/1e9; self.PdBm; self.RIN; self.linewidth/1e3];
            Units = {'nm'; 'GHz'; 'dBm'; 'dB/Hz'; 'kHz'};
            LaserTable = table(Parameters, Values, Units, 'RowNames', rows);
        end
        
        %% Main Impairment Methods
        function sigma2 = varRIN(self, P, Df)
            %% Compute Relative Intensity Noise (RIN) variance over given bandwidth
            % Inputs:
            % - P:  Optical power vector (W)
            % - Df: One-sided evaluation bandwidth (Hz)
            sigma2 = [];
            if ~isempty(self.RIN)  
                sigma2 = 10^(self.RIN/10)*P.^2*Df; 
            end
        end
        
        function sigma2 = varPN(self, fs)
            %% Compute white frequency noise variance per sample (Wiener phase noise model)
            % Inputs:
            % - fs: Sampling frequency (Hz)
            sigma2 = [];
            if ~isempty(self.linewidth)
                sigma2 = 2*pi*self.linewidth/fs; 
            end
        end
        
        function Pout = addIntensityNoise(self, Pout, fs)
            %% Induce Gaussian Relative Intensity Noise onto the time-domain power vector
            Pout = Pout + sqrt(self.varRIN(Pout, fs/2)).*randn(size(Pout));
        end
        
        function [Eout, phase_noise] = addPhaseNosie(self, Eout, fs)
            %% Apply random-walk phase noise matching the designated laser linewidth
            initial_phase = pi*(2*rand(1)-1); 
            dtheta = [0, sqrt(self.varPN(fs))*randn(1, length(Eout)-1)]; 
            phase_noise = initial_phase + cumsum(dtheta, 2);
            Eout = Eout.*exp(1j*phase_noise); 
        end
        
        function Eout = addTransientChirp(self, Ein)
            %% Model transient chirp phase variations driven by the output optical power profiles
            Eout = Ein.*exp(1j*self.alpha/2*log(abs(Ein).^2));
        end
        
        %% Signal Generation Methods
        function Eout = cw(self, sim)
            %% Generates a continuous-wave (CW) laser field sequence with impairment tracking
            Pout = self.PW*ones(size(sim.t));
            
            if isfield(sim, 'RIN') && sim.RIN
                Pout = self.addIntensityNoise(Pout, sim.fs);
            end
            
            Eout = sqrt(Pout);
            
            if isfield(sim, 'phase_noise') && sim.phase_noise
                Eout = self.addPhaseNosie(Eout, sim.fs);
            end
            
            if length(self.freqOffset) == 1
                if self.freqOffset ~= 0
                    Eout = Eout.*exp(1j*2*pi*self.freqOffset*sim.t);
                end
            else
                Eout = Eout.*exp(1j*2*pi*self.freqOffset.*sim.t);
            end            
        end    
        
        function Eout = modulate(self, x, sim)
            %% Apply direct modulation onto the source field driven by signal vector x
            Pout = self.PW*x/mean(x); 
            Pout = real(ifft(fft(Pout).*ifftshift(self.H(sim.f))));
                                  
            if isfield(sim, 'RIN') && sim.RIN
                Pout = self.addIntensityNoise(Pout, sim.fs);
            end
            
            if self.alpha ~= 0
                Pout = self.addTransientChirp(Pout);
            end
            
            Pout(Pout < 0) = 0; 
            Eout = sqrt(Pout); 
        end          
    end
    
    %% Dependent Property Accessors
    methods 
        function PW = get.PW(self)
            %% Convert log power parameter (dBm) to linear scale (W)
            PW = dBm2Watt(self.PdBm);
        end
        
        function self = set.PW(self, Ptx)
            %% Standardize input linear power (W) back to operational dBm values
            self.PdBm = Watt2dBm(Ptx);
        end
        
        function self = setPower(self, PtxW)
            %% Alternative function wrapper to explicitly define laser launch power in Watts
            self.PdBm = Watt2dBm(PtxW);
        end
        
        function lamb = get.wavelength(self)
            %% Map wavelength query back to core internal lambda property
            lamb = self.lambda;
        end
        
        function self = set.wavelength(self, lamb)
            %% Standardize external user-defined wavelength values symmetrically
            self.lambda = lamb;
        end
    end
end