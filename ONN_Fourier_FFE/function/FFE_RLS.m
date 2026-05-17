function [weights, equalizedSig] = FFE_RLS(noisySig, trainSig, nFwdWts, refTap, forgetFactor)
    % FFE_RLS: Trains Feed-Forward Equalizer (FFE) weights using the RLS algorithm

    %% 1. Parameter Initialization
    lambda = forgetFactor;            % Forgetting factor
    delta = 1e-2;                     % Initialization factor for covariance matrix
    P = delta * eye(nFwdWts);         % Initialize inverse correlation matrix
    weights = zeros(nFwdWts, 1);      % Initialize equalizer weights
    N = length(noisySig);             % Total signal length

    % Ensure the input signal length is sufficient
    if N < nFwdWts
        error('Input signal length must be greater than or equal to the number of FFE weights.');
    end

    %% 2. RLS Adaptive Training Loop
    for n = nFwdWts:N
        % Extract the input vector (current memory buffer for FFE taps)
        x = noisySig(n:-1:n-nFwdWts+1);
        
        % Calculate current equalizer output
        y = weights' * x;
        
        % Compute the prior estimation error
        err = trainSig(n-refTap+1) - y;
        
        % Compute the Kalman gain vector
        k = (P * x) / (lambda + x' * P * x);
        
        % Update equalizer weights
        weights = weights + k * err;
        
        % Update the inverse correlation matrix (Riccati equation)
        P = (P - k * x' * P) / lambda;
    end

    %% 3. Signal Equalization
    % Filter the input signal using the final trained weights
    equalizedSig = filter(weights, 1, noisySig);
end