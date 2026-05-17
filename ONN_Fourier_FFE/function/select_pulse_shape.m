function pulse_shape = select_pulse_shape(type, sps, varargin)
%% Generate struct pulse_shape containing pulse shaping parameters
% Inputs:
% - type: Either 'rect' (rectangular), 'rrc' (root-raised cosine), or 'rc' (raised cosine)
% - sps: Samples per symbol
% - varargin{1} (rolloff): Rolloff factor (0, 1]. Required if type = 'rrc' or 'rc'
% - varargin{2} (span): Number of symbols over which the pulse shape spans. Required if type = 'rrc' or 'rc'

    pulse_shape.type = type; 
    pulse_shape.sps = sps; 
    
    if not(strcmpi(type, 'rect'))
        assert(length(varargin) == 2, 'select_pulse_shape: invalid number of inputs. If type = rc or rrc, rolloff and span must be provided.')
    end
    
    switch lower(pulse_shape.type) 
        case 'rect' 
            pulse_shape.rolloff = 1;
            pulse_shape.h = ones(1, pulse_shape.sps); 
        case 'rrc' 
            pulse_shape.rolloff = varargin{1}; 
            pulse_shape.span = varargin{2}; 
            pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'sqrt'); 
        case 'rc' 
            pulse_shape.rolloff = varargin{1}; 
            pulse_shape.span = varargin{2}; 
            pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'normal'); 
        otherwise
            error('select_pulse_shape: Invalid pulse shape type')
    end
end