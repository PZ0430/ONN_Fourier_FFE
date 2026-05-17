%% Positive-Negative Photodiode (PIN) Class
classdef pin < apd
    %% PIN photodiode class inheriting methods and properties from apd
    properties
        % Inherited from apd.m:
        % R    - Responsivity
        % Id   - Dark current  
        % BW   - Bandwidth
    end
    
    methods
        function self = pin(R, Id, BW)
            %% Constructor
            % Inputs:
            %   - R  : (Optional, default = 1 A/W) Responsivity
            %   - Id : (Optional, default = 0 A) Dark current
            %   - BW : (Optional, default = Inf) Bandwidth of 1st-order frequency response
     
            if ~exist('R', 'var')
                R = 1;
            end
            
            if ~exist('BW', 'var')
                BW = Inf;
            end
            
            if ~exist('Id', 'var')
                Id = 0;
            end
            
            % Call superclass constructor: apd(GaindB, ka, BW, R, Id)
            % PIN photodiode has 0 dB Gain and 0 impact ionization factor (ka)
            self@apd(0, 0, BW, R, Id);
        end       
    
        function PINtable = summary(self)
            %% Generate Parameter Summary Table
            disp('-- PIN class parameters summary:')
            rows = {'Responsivity'; 'Dark current'; 'Bandwidth'};
            Parameters = {'R'; 'Id'; 'BW'};
            Values = [self.R; self.Id*1e9; self.BW/1e9];
            Units = {'A/W'; 'nA'; 'GHz'};
            PINtable = table(Parameters, Values, Units, 'RowNames', rows);
        end
    end
end