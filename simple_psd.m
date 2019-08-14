function [freq_axis, psd_axis] = simple_psd(signal, Fs)
% function [freq_axis, psd_axis] = simple_psd(signal)
%
% Compute a simple PSD given the signal and sample frequency.
%
% David Huberdeau, 01/10/2019

x_ = signal;

if mod(length(x_), 2) > 0
    x = x_(1:(length(x_) - 1));
else
    x = x_;
end
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

freq_axis = freq;
psd_axis = 10*log10(psdx + eps); %eps keeps log(.) from being Inf if 0

