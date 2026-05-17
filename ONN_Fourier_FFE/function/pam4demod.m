function dataRX = pam4demod(datain)
dataRX = datain;

dataRX(datain < -1) = 0;
dataRX(datain >= -1 & datain < 0) = 1;
dataRX(datain >= 0 & datain < 1) = 2;
dataRX(datain >= 1) = 3;
end