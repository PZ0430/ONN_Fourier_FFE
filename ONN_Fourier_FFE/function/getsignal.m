function [Erx, xt_ref, dataTX] = getsignal (sim,Tx,mpam,Fibers)
%% Transmitter signal generation
mpam = mpam.unbias;
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); 

%% Pulse shape
xd = real(pammod(dataTX,mpam.M)); 
xd_up = upsample(xd, Tx.pulse_shape.sps);                   
rrcFilter = rcosdesign(Tx.pulse_shape.rolloff, Tx.pulse_shape.span, Tx.pulse_shape.sps, 'sqrt');  
xd = conv(xd_up, rrcFilter, 'same');            

%% DAC
Tx.DAC.offset = sim.Mct/mpam.pulse_shape.sps*(length(mpam.pulse_shape.h)-1)/2; 
xt = dac(xd, Tx.DAC, sim);

%% Driver
xt = Tx.Vgain*(xt - mean(xt)) + Tx.VbiasAdj*mean(xt); 
xt_ref = xt;

%% Optical Modulator
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);

Ecw = Tx.Laser.cw(sim);
Etx = mzm_pam4(Ecw, xt, Tx.Mod); 

Padj = sqrt(Tx.Ptx/mean(abs(Etx).^2));
if abs(Padj-1) > 0.01 
    Etx = Etx*Padj;
    fprintf('Laser output power was adjusted by %.2f to meet target launched power\n', Padj);
end

%% Fiber Link Propagation
Erx = Etx;
attdB = 0;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    attdB = attdB + fiberk.att(Tx.Laser.wavelength)*fiberk.L/1e3;
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); 
end

end