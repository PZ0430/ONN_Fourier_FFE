function [xa, xqzoh, xq] = dac(x, DAC, sim, verbose)
    % DAC: Digital-to-analog conversion (quantization, ZOH, and analog filtering)

    %% 1. Parameter Initialization
    if nargin < 4 || isempty(verbose)
        verbose = false;
    end

    %% 2. Quantization & Clipping
    if isfield(sim, 'quantiz') && sim.quantiz && ~isinf(DAC.resolution)
        enob = DAC.resolution;
        
        % Get clipping ratio
        rclip = 0;
        if isfield(DAC, 'rclip')
            rclip = DAC.rclip;
        end
        
        % Get excursion limits
        if isfield(DAC, 'excursion') && ~isempty(DAC.excursion)
            xmin = DAC.excursion(1);
            xmax = DAC.excursion(2);
        else
            xmin = min(x);
            xmax = max(x);
        end
        
        % Apply clipping boundaries
        xamp = xmax - xmin;
        xmin = xmin + xamp * rclip; 
        xmax = xmax - xamp * rclip;
        xamp = xmax - xmin;
        
        % Check clipping probability
        Pc = mean(x > xmax | x < xmin);
        if Pc ~= 0
            fprintf('DAC: clipping probability = %G\n', Pc);
        end
        
        % Quantize
        dx = xamp / (2^enob - 1);
        codebook = xmin:dx:xmax;
        partition = codebook(1:end-1) + dx/2;
        [~, xq] = quantiz(x, partition, codebook); 
    else
        xq = x;
    end

    %% 3. Zero-Order Hold (ZOH)
    Nhold = sim.Mct / DAC.ros;
    assert(floor(Nhold) == ceil(Nhold), ...
        'dac: oversampling ratio of DAC (DAC.ros) must be an integer multiple of oversampling ratio of continuous time (sim.Mct)');
    
    xqzoh = upsample(xq, Nhold);
    xqzoh = filter(ones(1, Nhold), 1, xqzoh); 

    %% 4. Analog Filtering & Delay Correction
    if ~isempty(DAC.filt)
        % Filter frequency response
        Hdac = ifftshift(DAC.filt.H(sim.f / sim.fs)); 
        
        % Calculate total time shift (ZOH group delay + user offset)
        offset = (Nhold - 1) / 2;
        if isfield(DAC, 'offset')
            offset = DAC.offset + offset;
        end
        Hshift = ifftshift(exp(1j * 2 * pi * sim.f / sim.fs * offset));

        % Apply filtering in frequency domain
        xa = real(ifft(fft(xqzoh) .* Hdac .* Hshift)); 
    else
        xa = xqzoh;
    end
                                    
    %% 5. Visualization
    if verbose
        Ntraces = 100;
        Nstart = sim.Ndiscard * sim.Mct + 1;
        Nend = min(Nstart + Ntraces * 2 * sim.Mct, length(xa));
        
        figure(301);
        
        subplot(221); box on;
        eyediagram(xqzoh(Nstart:Nend), 2 * sim.Mct);
        title('Eye diagram after ZOH');
        
        subplot(222); box on;
        eyediagram(xa(Nstart:Nend), 2 * sim.Mct);
        title('DAC output eye diagram');
        
        subplot(223); box on;
        plot(sim.f / 1e9, abs(fftshift(fft(xa - mean(xa)))).^2);
        xlabel('Frequency (GHz)');
        ylabel('|X(f)|^2');
        title('DAC output spectrum');
        
        a = axis;
        axis([-sim.fs/1e9, sim.fs/1e9, a(3), a(4)]);
        drawnow;
    end
end