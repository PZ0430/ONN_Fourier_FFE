function [f, t] = freq_time(N, fs)
%% Generate frequency and time coordinate vectors for discrete signals
% Inputs:
% - N  : Total number of signal sample points
% - fs : Sampling frequency (Hz)
% Outputs:
% - f  : Frequency vector centered from -fs/2 to fs/2-df (Hz)
% - t  : Time vector starting from 0 to (N-1)*dt (s)

    dt = 1/fs;
    t = 0:dt:(N-1)*dt;
    
    df = 1/(dt*N);
    f = -fs/2:df:fs/2-df;
end