# Optical Fiber Transmission Link Simulation with ONN & FFE-RLS Equalization

This MATLAB repository simulates a high-speed, direct-detection optical communication system using **4-PAM (PAM-4)** modulation. The simulation encompasses a complete transmitter (TX) processing chain, transmission over a Single-Mode Fiber (SMF), an Optical Neural Network (ONN) structure implemented via Fourier coefficients, a Photodiode (PD) receiver, an Analog-to-Digital Converter (ADC), and a Recurrent Least Squares Feed-Forward Equalizer (FFE-RLS) at the receiver DSP stage.

## System Architecture Overview

The simulation follows a classic intensity-modulation direct-detection (IM-DD) architecture enhanced with optical processing:

1. **Transmitter (TX):** Generates a 200 Gbps PAM-4 signal, applies pulse shaping (RRC), and passes it through an 8-bit DAC and an electrical Mach-Zehnder Modulator (MZM).
2. **Channel:** Simulates impairments over a Single-Mode Fiber (SMF) channel including Phase Noise, RIN, and PMD.
3. **Optical Neural Network (ONN):** Processes the optical field prior to detection using calculated channel Fourier coefficients to mitigate fiber nonlinearities/dispersion optically.
4. **Receiver (RX):** Photodiode direct-detection, Transimpedance Amplifier (TIA) thermal noise injection, and 8-bit ADC quantization.
5. **Rx DSP (Equalization):** FFE adaptive equalization powered by the RLS algorithm, followed by PAM-4 demodulation and Bit Error Rate (BER) counting.


## Key Parameters Configuration

You can alter the system performance directly within the script via these parameters:

| Category | Parameter | Default Value | Description |
| :--- | :--- | :--- | :--- |
| **Signal** | `PtxdBm` | `1 dBm` | Transmitter optical launch power. |
| | `fiblen` | `2000 m (2 km)` | Total Single-Mode Fiber transmission distance. |
| | `bitRate`| `200 Gbps` | Target net bit rate (100 Gbaud for PAM-4). |
| **Control**| `sim.quantiz` | `true` | Enforces hardware DAC/ADC bit limitations. |
| | `sim.RIN` / `sim.PMD` | `true` | Toggles laser Relative Intensity Noise & Polarization Mode Dispersion. |
| **ONN** | `parallel` | `1` | Decimation/parallel processing factor for ONN taps. |
| **DSP** | `ntaps` / `refTap` | `7 / 2` | FFE total taps and main cursor reference tap position. **More taps are needed for longer fiber length**|
| | `forgetFactor` | `0.999999` | RLS tracking algorithm memory factor ($\lambda$). |

---

## Quick Start

1. Place your main script in the root directory alongside your dependency folders.
2. Open MATLAB and run your main execution script.
3. The command window will output the received optical power and the final System Bit Error Rate (BER):

