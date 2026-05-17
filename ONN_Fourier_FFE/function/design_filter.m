function filt = design_filter(type, order, fcnorm, verbose)
%% Design continuous-time analog filters and map them to discrete-time domains
% Inputs:
% - type: Filter family selection {'butter', 'cheby1', 'ellipt', 'two-pole', 'matched', 'fir1', 'gaussian', 'bessel', 'lorentzian', 'fbg'}
% - order: Filter order specification (or functional handle depending on selection)
% - fcnorm: Normalized cutoff frequency, defined as f3dB / (fs/2), where 0 < fcnorm < 1
% - verbose: Boolean flag to enable visual performance plots
% Output:
% - filt: Struct containing coefficients (num, den), group delay, impulse response h, and analytical bandwidth tracking

    if not(exist('verbose', 'var'))
        verbose = false;
    end

    maxMemoryLength = 2^10; 
    threshold = 1 - 1e-6;      
    type = lower(type); 

    switch type       
        case 'butter'
            [num, den] = butter(order, fcnorm);
         
        case 'bessel'       
            % The conversion factor wc2wo maps the -3 dB cutoff frequency to the 
            % Bessel constant delay optimization bound. Derived solving:
            % abs(polyval(b, 1j*2*pi*f3dB)/polyval(a, 1j*2*pi*f3dB))^2 == 0.5
            wc2wo = 1.621597623980423; 
            wow = wc2wo*(2*pi*fcnorm/2);
            [nums, dens] = besself(order, wow);          
            [num, den] = impinvar(nums, dens, 1);

            if verbose
                verbose = false;
                f = linspace(0, 1/2);
                Hs = freqs(nums, dens, 2*pi*f);
                Hz = freqz(num, den, f, 1);
                plot_transform(Hs, Hz, f, fcnorm/2)
            end
    
        otherwise
            error('Unknown filter type = %s!', type)
    end

    filt.type = type;
    filt.order = order;
    filt.num = num;
    filt.den = den;
    filt.grpdelay = grpdelay(num, den, 1);
    filt.fcnorm = fcnorm;
    filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);

    if exist('nbw', 'var')
        filt.noisebw = nbw;
    else
        filt.noisebw = @(fs) noisebw(num, den, 2^15, fs); 
    end

    if den == 1 
        filt.h = num;
    else
        x = zeros(1, maxMemoryLength+1);
        x(1) = 1;
        y = filter(filt.num, filt.den, x);
        E = cumsum(abs(y).^2)/sum(abs(y).^2);
        y(E > threshold) = [];
        y = y/abs(sum(y)); 
        filt.h = y;
    end

    if verbose
        plot_filter(filt)
    end
end

