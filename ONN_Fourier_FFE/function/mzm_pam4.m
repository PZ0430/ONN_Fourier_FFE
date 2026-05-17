function Eout = mzm_pam4(Ein, Vin, Mod)
    %% Signal filtering and driver response application
    Hmod = ifftshift(Mod.H);    
    Vx = Vin(1, :);
    
    Vin_normalized = (Vx - min(Vx)) / (max(Vx) - min(Vx));
    Vx = real(ifft(fft(Vin_normalized) .* Hmod));
    
    %% Optical modulation
    Enorm = 1;  
    Vout = Enorm * sin(pi * Vx / 2);    
    
    Ein = Ein(1:length(Vout));
    Eout = Ein .* Vout;
end