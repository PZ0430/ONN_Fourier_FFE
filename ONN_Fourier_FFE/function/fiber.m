classdef fiber < handle
    %% Single-mode fiber class defining physical properties and propagation behavior
    properties
        L 
        att 
        D 
        PMD = false; 
        gamma = 1.4e-3 
        meanDGDps = 0.1; 
        PMD_section_length = 1e3  
        PDL = 0 
    end
    
    properties(Dependent)
        tauDGD 
    end
        
    properties(GetAccess=protected)
        JonesMatrix 
    end
    
    properties(Constant) 
        S0 = 0.092*1e3;    
        lamb0 = 1310e-9;   
    end
    
    properties(Constant, Hidden)
        c = 299792458;  
    end    
    
    methods
        function obj = fiber(L, att, D)
            %% Class constructor 
            % Inputs: 
            % - L: Fiber length (m)
            % - att (optional): Function handle of fiber attenuation (dB/km) as a function of wavelength
            % - D (optional): Function handle of fiber chromatic dispersion (s/m^2) as a function of wavelength
            if nargin >=1
                obj.L = L;
            else
                obj.L = 0;
            end
            
            if nargin >= 2 
                obj.att = att;
            else 
                obj.att = @(lamb) 0;
            end
            
            if nargin == 3
                obj.D = D;
            else          
                obj.D = @(lamb) fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); 
            end           
        end
        
        function Fibertable = summary(self, lamb)
            %% Generate a summary table of the fiber class parameters
            disp('-- Fiber class parameters summary:')
            lambnm = lamb*1e9;
            rows = {'Length'; sprintf('Attenuation at %.2f nm', lambnm);...
                sprintf('Total dispersion at %.2f nm', lambnm); 'PMD included?';...
                'Total DGD'; 'Polarization dependent loss'};
            Parameters = {'L'; 'att'; 'DL'; 'PMD'; 'tauDGD'; 'PDL'};
            Values = [self.L/1e3; self.att(lamb); self.D(lamb)*self.L*1e3;...
                self.PMD; self.tauDGD*1e12; self.PDL];
            Units = {'km'; 'dB/km'; 'ps/nm'; ''; 'ps'; 'dB'};
            Fibertable = table(Parameters, Values, Units, 'RowNames', rows);
        end
             
        function tauDGD = get.tauDGD(self)
            %% Calculate total Differential Group Delay (DGD) induced by PMD
            tauDGD = double(self.PMD)*self.meanDGDps*1e-12*sqrt(self.L/1e3); 
        end

        function [Ncd, Npmd] = Ntaps(self, Rs, ros, lambda)
            %% Estimate the number of DSP taps required to compensate for CD and PMD in a coherent link
            % Based on Ip, E., & Kahn, J. M. (2007). Digital equalization of chromatic dispersion and polarization mode dispersion. 
            % Journal of Lightwave Technology, 25(8), 2033-2043.
            Ncd = 2*pi*abs(self.beta2(lambda))*self.L*Rs^2*ros; 
            Npmd = self.tauDGD*ros*Rs; 
        end
        
        function b2 = beta2(this, lamb)
            %% Calculate the second-order dispersion coefficient (beta2)
            b2 = -this.D(lamb).*(lamb.^2)/(2*pi*this.c); 
        end    
        
        function Leff = effective_length(self, lamb)
            %% Calculate the effective fiber length accounting for attenuation
            alpha = log(10)*self.att(lamb)/1e4; 
            Leff = (1 - exp(-alpha*self.L))/alpha;
        end
        
        function [link_att, link_attdB] = link_attenuation(this, lamb)
            %% Calculate link attenuation in both linear units and decibels (dB)
            % Input:
            % - lamb: Wavelength (m)
            link_att = 10.^(-this.att(lamb)*this.L/1e4); 
            link_attdB = this.L/1e3*this.att(lamb);
        end
        
        function Eout = linear_propagation(this, Ein, f, lambda)
            %% Simulate linear propagation including chromatic dispersion and first-order PMD
            % Inputs: 
            % - Ein: Input electric field
            % - f: Frequency vector (Hz)
            % - lambda: Wavelength (m)
            % Outputs: 
            % - Eout: Output electric field
            if this.L*this.D(lambda) == 0
                Eout = Ein;
                return
            end
            
            Ein = this.enforceDimConvention(Ein);
            f = this.enforceDimConvention(f);
         
            two_pols = (size(Ein, 1) == 2);
            Einf = fftshift(fft(Ein, [], 2), 2);
            
            if this.PMD
                if isempty(this.JonesMatrix) 
                    this.generateJonesMatrix(2*pi*f); 
                end
                
                if not(two_pols)
                    Einf = [Einf; zeros(size(Einf))];
                    two_pols = true;
                end
                for k = 1:length(f)
                    Einf(:, k) = this.JonesMatrix(:,:,k)*Einf(:, k);
                end
            end

            Hele = this.Hdisp(f, lambda);
            if two_pols
                Eout = Einf;
                Eout(1, :) = ifft(ifftshift(Hele.*Einf(1, :)));
                Eout(2, :) = ifft(ifftshift(Hele.*Einf(2, :)));
            else
                Eout = ifft(ifftshift(Hele.*Einf));
            end

            if two_pols && this.PDL ~= 0
                a = 10^(-this.PDL/10);
                Eout(2, :) = a*Eout(2, :);
            end
            
            Eout = Eout*sqrt(this.link_attenuation(lambda));
        end
        
        function Hele = Hdisp(this, f, lambda)
            %% Calculate chromatic dispersion frequency response: Hele(f) = Eout(f)/Ein(f)
            % Inputs: 
            % - f: Frequency vector (Hz)
            % - lambda: Wavelength (m)
            % Outputs:
            % - Hele: Frequency response matrix
            beta2 = this.beta2(lambda);
            w = 2*pi*f;
            Dw = -1j*1/2*beta2*(w.^2);
            Hele = exp(this.L*Dw);
        end
        
        function h = hdisp(self, t, lambda)
            %% Calculate chromatic dispersion impulse response via inverse Fourier transform of Hdisp(f)
            % Inputs: 
            % - t: Time vector (s)
            % - lambda: Wavelength (m)
            b = 1/2*self.beta2(lambda)*self.L;
            h = sqrt(-pi*1j/b)/(2*pi)*exp(1j*t.^2/(4*b));
        end
               
        function Hf = Himdd(self, f, wavelength, alpha, type)
            %% Calculate fiber power frequency response for an IM-DD system with transient chirp dominance
            % Transfer function defines optical power ratio: Hfiber(f) = Pout(f)/Pin(f)
            % Inputs:
            % - f: Frequency vector (Hz)
            % - wavelength: Wavelength (m)
            % - alpha (optional): Chirp parameter (alpha > 0 for DML)
            % - type (optional): Type of frequency response ('small signal' or 'large signal')
            if not(exist('alpha', 'var')) 
                alpha = 0;
            end
            
            beta2 = self.beta2(wavelength);
            theta = -1/2*beta2*(2*pi*f).^2*self.L; 
            
            if exist('type', 'var') && strcmpi(type, 'large signal')
                mIM = 0.7; 
                Dphi = pi/2; 
                mFM = alpha/2*mIM; 
                u = 2*mFM*sin(theta);
                Hf = cos(theta).*(besselj(0, u) - besselj(2, u)*exp(1j*Dphi)) - 2*exp(1j*Dphi)/(1j*mIM)*besselj(1, u);                  
            elseif not(exist('type', 'var')) || strcmpi(type, 'small signal')
                Hf = cos(theta) - alpha*sin(theta);  
            else
                error('fiber/Hf: undefined type of fiber frequency response')
            end
        end
        
        function tau = calcDGD(self, omega)
            %% Calculate differential group delay (DGD) from the Jones Matrix
            if ~self.PMD 
                tau = zeros(size(omega));
                warning('fiber/calcDGD: PMD is disabled')
                return
            end
                
            if isempty(self.JonesMatrix)
                tau = [];
                warning('fiber/calcDGD: Jones Matrix has not been calculated yet')
                return
            end
            
            tau = zeros(1,length(omega));
            dw = abs(omega(1)-omega(2));
            for m = 1:length(omega)-1
                tau(m) = 2/dw*sqrt(det(self.JonesMatrix(:,:,m+1)-self.JonesMatrix(:,:,m)));
            end
            tau(end) = tau(end-1);
         end
        
    end  
    
    methods(Access=private)
        function V = enforceDimConvention(~, V)
            %% Ensure vector dimensions comply with standard matrix orientation requirements
            if size(V, 1) > size(V, 2)
                V = V.';
            end
        end
        
        function M = generateJonesMatrix(self, omega)
            %% Generate random Jones Matrix representing PMD sections
            Nsect = ceil(self.L/self.PMD_section_length);
                       
            dtau = self.tauDGD/sqrt(Nsect);
            M = zeros(2,2,length(omega));
            M(:, :, 1) = randomRotationMatrix(); 
            
            for k = 1:Nsect
                U = randomRotationMatrix();
                for m = 2:length(omega)
                    Dw = [exp(1j*dtau*omega(m)/2), 0; 0, exp(-1j*dtau*omega(m)/2)]; 
                    M(:,:,m) = U*Dw*U';
                end
            end
            
            self.JonesMatrix = M;
            
            function U = randomRotationMatrix()
                phi = rand(1, 3)*2*pi;
                U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
                U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
                U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];
                U = U1*U2*U3;
            end            
        end
    end    
end