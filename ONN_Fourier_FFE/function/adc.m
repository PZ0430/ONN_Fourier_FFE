function [ys, varQ, xa] = adc(x, ADC, sim, Hrx)
    % ADC: Analog-to-digital conversion (filtering, downsampling, and quantization)

    %% 1. Receiver Filter Response
    if ~exist('Hrx', 'var') || isempty(Hrx)
        % Use default antialiasing filter if Hrx is not provided
        Hrx = ifftshift(ADC.filt.H(sim.f / sim.fs)); 
    else
        % Use provided filter (typically a matched filter)
        Hrx = ifftshift(Hrx); 
    end
    
    % Ensure Hrx vector length matches the input signal
    Hrx = Hrx(1:length(x));

    %% 2. Frequency-Domain Delay Shift & Filtering
    if isfield(ADC, 'offset') && ADC.offset ~= 0
        Hshift = ifftshift(exp(1j * 2 * pi * sim.f / sim.fs * ADC.offset));
    else
        Hshift = 1;
    end
    
    % Apply filtering and transform back to time domain
    xa = real(ifft(fft(x) .* Hrx .* Hshift));

    %% 3. Time Reference Alignment via Cross-Correlation
    if isfield(ADC, 'timeRefSignal')
        [c, lags] = xcorr(ADC.timeRefSignal, xa);
        [~, idx] = max(abs(c));
        if lags(idx) ~= 0
            xa = circshift(xa, [0, lags(idx)]);
            fprintf('ADC: input signal was delayed by %d samples (%.2f ps) to match time reference signal\n', ...
                lags(idx), 1e12 * lags(idx) / sim.fs);
        end
    end

    %% 4. Downsampling
    % Check if resampling is required based on rates compatibility
    if mod(sim.Mct, ADC.ros) == 0
        xs = xa(1 : sim.Mct/ADC.ros : end);
    else
        [N, D] = rat(ADC.ros / sim.Mct);
        fprintf('ADC: sim.Mct/ADC.ros is not integer, so signal was resampled by %d/%d.\n', N, D);
        xs = resample(xa, N, D);
    end

    %% 5. Quantization & Clipping
    if isfield(sim, 'quantiz') && sim.quantiz && ~isinf(ADC.ENOB)
        enob = ADC.ENOB;
        
        % Get clipping ratio
        rclip = 0;
        if isfield(ADC, 'rclip')
            rclip = ADC.rclip;
        end
        
        % Get excursion boundaries
        if isfield(ADC, 'excursion') && ~isempty(ADC.excursion)
            xmin = ADC.excursion(1);
            xmax = ADC.excursion(2);
        else
            xmin = min(xs);
            xmax = max(xs);
        end
        
        % Apply clipping
        xamp = xmax - xmin;    
        xmin = xmin + xamp * rclip; 
        xmax = xmax - xamp * rclip;
        xamp = xmax - xmin;
        
        % Check clipping probability
        Pc = mean(xs > xmax | xs < xmin);
        if Pc ~= 0
            fprintf('ADC: clipping probability = %G\n', Pc);
        end
        
        % Perform quantization
        dx = xamp / (2^enob - 1);
        codebook = xmin:dx:xmax;
        partition = codebook(1:end-1) + dx/2;
        [~, ys, varQ] = quantiz(xs, partition, codebook); 
        
        % Display quantization noise metrics
        fprintf('ADC: quantization noise variance = %G | Signal-to-quantization noise ratio = %.2f dB\n', ...
            varQ, 10 * log10(var(xs) / varQ));
    else
        ys = xs;
        varQ = 0;
    end
end