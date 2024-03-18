function [x,p] = SolveViscodamperOscillatorFreqDom(m,k,lambd,C0,q,r,...
                                                   f,tsample)

% m,k,lambd,C0,q,r are parameters of the viscodamper oscillator
% tsample = sample time for numerical solution
% f = forcing function sampled at tsample; length(f) is to be even 

nfft = length(f);
dfreq = 1/tsample/nfft; % frequency resolution
omeg = (-nfft/2:nfft/2-1)'*dfreq*2*pi;

F = fftshift(fft(f)); % DFT of forcing function

kbar = C0*(1i*omeg).^q./(1 + lambd*(1i*omeg).^r);

X = F./(-m*omeg.^2 + k + kbar);
P = kbar.*X;

x = real(ifft(ifftshift(X))); % imag part is due to roundoff and is removed
p = real(ifft(ifftshift(P)));
