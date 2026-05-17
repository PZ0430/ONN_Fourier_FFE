%% Compute analytical mapping coefficients using Fourier expansion for fiber dispersion
function [Cn, n_99] = get_fourier_coefficients(L)
    %% Physical parameters 
    c = 299792458;
    lambda = 1550e-9;
    D = 17e-6;  
    beta2 = -D*lambda^2/(2*pi*c);  
    tau = 10e-12;  

    %% Frequency settings 
    N = 25;                         
    Nw = 200001;                     
    w = linspace(-pi/tau, pi/tau, Nw);
    omega = linspace(-pi, pi, Nw).';    
    
    %% Target transfer function 
    H = exp(1j * beta2 * L / 2 .* w.^2);

    %% Compute Fourier coefficients (Cn)
    Cn = zeros(2*N+1, 1);
    n_vec = (-N:N).';
    domega = omega(2) - omega(1);
    omega = omega(:)';
    for k = 1:length(n_vec)
        n = n_vec(k);
        Cn(k) = (1/(2*pi)) * sum( H .* exp(-1j*n*omega) ) * domega;
    end

    %% Tap statistics and energy coverage analysis
    Cn_mag = abs(Cn);
    energy = Cn_mag.^2;
    energy_norm = energy / sum(energy);
    
    [~, idx_sort] = sort(abs(n_vec));
    energy_sorted = energy_norm(idx_sort);
    cum_energy = cumsum(energy_sorted);
    idx_99 = find(cum_energy >= 0.99, 1);
    n_99 = max(abs(n_vec(idx_sort(1:idx_99))));
    num_99_taps = 2*n_99 + 1;

    fprintf('-----------------------------------------\n');
    fprintf('Number of center taps covering 99%% energy: %d\n', num_99_taps);
    fprintf('Corresponding tap range: n = [-%d, %d]\n', n_99, n_99);
    fprintf('-----------------------------------------\n');
end