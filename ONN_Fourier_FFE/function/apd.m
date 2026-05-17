%% Avalanche Photodiode (APD) Class
classdef apd
    properties
        Gain    % Linear gain
        ka      % Impact ionization factor
        BW0     % Low-gain bandwidth
        GainBW  % Gain-bandwidth product
        R       % Responsivity (A/W)
        Id      % Dark current (A)
    end
    
    properties (Dependent)
        Fa      % Excess noise factor
        Geff    % Effective gain = Gain * Responsivity
        BW      % Operational bandwidth
        GaindB  % Gain in dB 
    end
    
    properties (Constant, Hidden)
        h = 6.62606957e-34; % Planck constant
        q = 1.60217657e-19; % Electron charge
        c = 299792458;      % Speed of light
    end
    
    properties (Constant, GetAccess=private, Hidden)
        cdf_accuracy = 1-1e-4; % Required accuracy of the CDF for PMF truncation
        Niterations = 1e6;     % Maximum number of iterations in while loops
        Ptail = 1e-6;          % Probability of clipped tail
    end
    
    properties (Dependent, GetAccess=private, Hidden)
        a % Auxiliary variable: G = 1/(1-ab)
        b % Auxiliary variable: b = 1/(1-ka)
    end
   
    methods
        function this = apd(GaindB, ka, BW, R, Id)
            %% Class Constructor
            % Inputs:
            %   - GaindB : Gain in dB
            %   - ka     : Impact ionization factor
            %   - BW     : (Optional, default = Inf) If scalar, specifies bandwidth.
            %              If 2x1 vector, [low-gain BW, gain-bandwidth product].
            %   - R      : (Optional, default = 1 A/W) Responsivity
            %   - Id     : (Optional, default = 0 A) Dark current
            
            this.Gain = 10^(GaindB/10);
            this.ka = ka;
            
            if exist('BW', 'var')
                if length(BW) == 1
                    this.BW0 = BW;
                    this.GainBW = Inf;
                else
                    this.BW0 = BW(1);
                    this.GainBW = BW(2);
                end
            else
                this.BW0 = Inf;
                this.GainBW = Inf;
            end
                        
            if exist('R', 'var')
                this.R = R;
            else 
                this.R = 1;
            end
            
            if exist('Id', 'var')
                this.Id = Id;
            else 
                this.Id = 0;
            end
        end
        
        function APDtable = summary(self)
            %% Generate Parameter Summary Table
            disp('-- APD class parameters summary:')
            rows = {'Gain'; 'Impact ionization factor'; 'Low-gain bandwidth';...
                'Gain-bandwidth product'; 'Responsivity'; 'Dark current'};
            Parameters = {'Gain'; 'ka'; 'BW0'; 'GainBW'; 'R'; 'Id'};
            Values = [self.Gain; self.ka; self.BW0/1e9; self.GainBW/1e9; self.R; self.Id*1e9];
            Units = {''; ''; 'GHz'; 'GHz'; 'A/W'; 'nA'};
            APDtable = table(Parameters, Values, Units, 'RowNames', rows);
        end
        
        %% Main Methods
        function l = lambda(this, P, dt)
            %% Poisson Process Rate
            % Calculates the rate for a given power P in an interval dt.
            % Ref: Personick, "Statistics of a General Class of Avalanche Detectors..."
            l = (this.R*P + this.Id)*dt/this.q; 
        end
        
        function Hapd = H(this, f)
            %% Frequency Response
            % Normalized to have unit gain at DC.
            if isinf(this.BW)
                Hapd = ones(size(f));
            else
                % Hapd = 1./sqrt(1 + (f/this.BW).^2);
                Hapd = 1./(1 + 1j*(f/this.BW))...
                    .*exp(1j*2*pi*f.*(this.BW./(2*pi*(this.BW^2 + f.^2)))); % Remove group delay
            end
        end
        
        function hapd = ht(this, t)
            %% Impulse Response
            if isinf(this.BW)
                hapd = double(t == 0);
            else
                % Hapd = 1./sqrt(1 + (f/this.BW).^2);
                hapd = this.BW*exp(-t*this.BW);
                hapd(t < 0) = 0;
            end            
        end
        
        function [Hw, y] = Hwhitening(self, f, P, N0, x)
            %% Whitening Filter Design
            % Inputs:
            %   - f  : Frequency vector
            %   - P  : Power at the APD input (average power used for shot noise PSD)
            %   - N0 : Power spectral density (PSD) of thermal noise
            %   - x  : (Optional) Input signal to be filtered by Hw
            if ~isinf(self.BW)
                r = N0/self.varShot(P, 1); % Ratio between thermal and shot noise PSD
                Hw = sqrt((1 + r)./(r + abs(self.H(f)).^2));
                [Hw, groupdelay] = Hgrpdelay(Hw, f);
                if exist('x', 'var')
                    y = ifft(fft(x).*ifftshift(Hw));
                end
            else
                Hw = 1;
                y = [];
                if exist('x', 'var')
                    y = x;
                end
            end
        end
        
        function bw = noisebw(this) 
            %% One-Sided Noise Bandwidth
            % Factor of pi/2 arises because this.Hapd energy is pi*BW
            bw = pi/2*this.BW; 
        end
        
        function sig2 = varShot(this, Pin, Df)
            %% Shot Noise Variance
            % Ref: Agrawal 4.4.17 (4th edition)
            % Inputs:
            %   - Pin : Input optical power (W)
            %   - Df  : Noise bandwidth (Hz)
            sig2 = 2*this.q*this.Gain^2*this.Fa*(this.R*Pin + this.Id)*Df; 
        end
        
        function output = detect(this, Ein, fs, noise_stats, N0)
            %% Direct Detection Simulation
            % Inputs:
            %   - Ein         : Input electric field
            %   - fs          : Sampling frequency (Hz)
            %   - noise_stats : 'gaussian', 'doubly-stochastic' (deprecated), or 'no noise'
            %   - N0          : (Optional) Thermal noise PSD added after detection
            
            if any(size(Ein) == 2) % Dual polarization support
                if size(Ein, 1) ~= 2
                    Ein = Ein.'; % Convert to 2 x N format
                end
                Pin = sum(abs(Ein).^2, 1);
            else
                Pin = abs(Ein).^2;
            end     
            
            switch noise_stats 
                case 'gaussian'
                    output = this.R*this.Gain*Pin + sqrt(this.varShot(Pin, fs/2)).*randn(size(Pin));
                    
                case 'doubly-stochastic'
                    error('apd/detect: noise_stats = doubly-stochastic is deprecated.');
                    % % Uses saddlepoint approximation to obtain PMF for random generation
                    % Plevels = unique(Pin);
                    % output = zeros(size(Pin));
                    % for k = 1:length(Plevels)
                    %     [px, x] = this.output_pdf_saddlepoint(Plevels(k), fs, 0);
                    %     cdf = cumtrapz(x, px);
                    %     pos = (Pin == Plevels(k));
                    %     u = rand(sum(pos), 1);
                    %     dist = abs(bsxfun(@minus, u, cdf));
                    %     [~, ix] = min(dist, [], 2);
                    %     output(pos) = x(ix);
                    % end
                    
                case 'no noise'
                    % Only amplifies and filters the signal (for BER estimation and debugging)
                    output = this.R*this.Gain*Pin;
                    
                otherwise 
                    error('apd/detect: Invalid Option!')
            end
            
            % Frequency vector setup
            df = fs/length(Pin); % Calculate frequency resolution
            f = (-fs/2:df:fs/2-df);
            if size(output, 1) > size(output, 2)
                f = f.';
            end
            
            % Apply APD frequency response (H has unit gain at DC)
            if ~isinf(this.BW)
                output = real(ifft(fft(output).*ifftshift(this.H(f))));
            end
                        
            % Add thermal noise if N0 is provided
            if exist('N0', 'var')
                output = output + sqrt(N0*fs/2).*randn(size(Pin)); 
            end 
        end
              
        function noise_std = stdNoise(this, Hrx, Hff, N0, RIN, sim)
            %% Noise Standard Deviation Calculator
            % Returns a function handle to calculate noise std for a given power level P.
            % Note: Power level P is assumed to be referred to AFTER the APD.
            % Thermal, shot, and RIN are assumed to be white AWGN.
            % Inputs:
            %   - Hrx : Receiver filter response (e.g., whitening + matched filter)
            %   - Hff : Equalizer frequency response
            %   - N0  : Thermal noise PSD
            %   - RIN : RIN in dB/Hz (omitted if empty)
            %   - sim : Simulation configuration struct
            
            Htot = Hrx.*Hff; 
            Df = 1/2*trapz(sim.f, abs(Htot).^2); % Filter noise BW (includes noise enhancement)
            
            if ~isinf(this.BW)
                Dfshot = 1/2*trapz(sim.f, abs(this.H(sim.f).*Htot).^2);
            else
                Dfshot = Df;
            end
            DfRIN = Dfshot; % Approximated: needs to include fiber response
            
            % 1. Thermal noise variance
            varTherm = N0*Df; 
            
            % 2. RIN variance
            if ~isempty(RIN) && isfield(sim, 'RIN') && sim.RIN
                varRIN = @(Plevel) 10^(RIN/10)*Plevel.^2*DfRIN;
            else
                varRIN = @(Plevel) 0;
            end
            
            % 3. Shot noise variance 
            % Plevel is divided by Geff to get power at APD input
            varShot = @(Plevel) this.varShot(Plevel/(this.Gain*this.R), Dfshot);
            
            % Combined Noise Standard Deviation Function Handle
            noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + varShot(Plevel));
        end
        
        function [Gopt, mpam] = optGain(this, mpam, tx, fiber, rx, sim)
            %% Optimize APD Gain
            % Finds the APD gain that minimizes required optical power for a target BER.
            disp('Optimizing APD gain for sensitivity');
                   
            if mpam.optimize_level_spacing
                % Level spacing optimization ensures target BER is met for a given gain
                [Gopt, ~, exitflag] = fminbnd(@(Gapd) ...
                    this.optimize_PAM_levels(Gapd, mpam, tx, fiber, rx, sim), eps, maxGain(this, mpam.Rs/5));  
                [~, mpam] = this.optimize_PAM_levels(Gopt, mpam, tx, fiber, rx, sim);
            else
                % Adjust power to reach target BER
                [Gopt, ~, exitflag] = fminbnd(@(Gapd) fzero(@(PtxdBm)...
                    log10(this.calc_apd_ber(PtxdBm, Gapd, mpam, tx, fiber, rx, sim)) - log10(sim.BERtarget), -20), 1, maxGain(this, mpam.Rs/5));
            end
            
            if exitflag ~= 1
                warning('apd/optGain: APD gain optimization did not converge (exitflag = %d)\n', exitflag);
            end 
            assert(Gopt >= 0, 'apd/optGain: Negative gain found while optimizing APD gain')
                    
            % Helper function for maximum allowed gain
            function Gmax = maxGain(apd, minBW)
                if isinf(apd.GainBW)
                    Gmax = 100;
                else
                    Gmax = apd.GainBW/minBW; 
                end
            end
        end 
        
        function [Pmean, mpam] = optimize_PAM_levels(this, Gapd, mpam, Tx, Fiber, Rx, sim)
            %% Calculate Optimal Level Spacing
            % Performed iteratively since the whitening filter depends on optical power.
            this.Gain = Gapd;
                       
            Pmean = [0 1];
            tol = 1e-6;
            n = 1;
            maxIterations = 50;
            Tx.Ptx = 1e-3;       % Initial power 
            Fiber.att = @(l) 0;  % Disregard attenuation
            Pdiff = Inf;
            mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB); % Starting level spacing
            
            while Pdiff > tol && n < maxIterations  
                [~, noise_std] = ber_apd_awgn(mpam, Tx, Fiber, this, Rx, sim);
                % Optimize levels at the transmitter based on receiver noise std
                mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, Tx.Mod.rexdB, noise_std); 
                
                % Required power at the APD input
                Pmean(n+1) = mean(mpam.a)/this.Geff;
                Tx.Ptx = Pmean(n+1);
                Pdiff = abs((Pmean(n+1) - Pmean(n))/Pmean(n));
                n = n+1; 
            end
            
            Pmean = Pmean(end);
            if n >= maxIterations
                warning('apd/optimize_PAM_levels: optimization did not converge')
            end
        end                
    end
           
    methods (Access=private)     
        function ber = calc_apd_ber(this, PtxdBm, Gapd, mpam, Tx, Fiber, Rx, sim)
            %% Iterative BER Calculation
            % Used only when M-PAM has equally spaced levels.
            Tx.Ptx = dBm2Watt(PtxdBm);
            this.Gain = Gapd; % Linear units
            ber = ber_apd_enumeration(mpam, Tx, Fiber, this, Rx, sim);
        end
    end
    
    methods
        %% Getters and Setters
        function Fa = get.Fa(this)
            % Excess noise factor (Agrawal 4.4.18, 4th edition)
            Fa = this.ka*this.Gain + (1 - this.ka)*(2 - 1/this.Gain); 
        end
               
        function GaindB = get.GaindB(this) 
            GaindB = 10*log10(this.Gain);
        end
        
        function BW = get.BW(this)
            % Operational bandwidth (limited by avalanche buildup time if GainBW/Gain < BW0)
            BW = min(this.BW0, this.GainBW/this.Gain);
        end
        
        function Geff = get.Geff(this)
            Geff = this.Gain*this.R;
        end
        
        function b = get.b(this) 
            b = 1/(1-this.ka);
        end
        
        function a = get.a(this) 
            a = 1/this.b*(1-1/this.Gain);
        end
               
        function this = set.GaindB(this, GdB)
            this.Gain = 10^(GdB/10); 
        end
        
        function self = setGain(self, Gain)
            self.Gain = Gain;
        end
    end
end

