% References: 
% [1] Makris, N. and Constantinou, M. C. (1991). 
%     "Fractional-derivative maxwell model for viscous dampers",
%     ASCE Journal of Structural Engineering, 117(9): 2708-2724

% Viscodamper model from [1] (pg 2711, eq (3))
% P + \lambda D^r P = C_0 D^q X
% Oscillator model [1] (pg 2714, eq (11))
% m\ddot(X) + k X + P = F

% Definitions (pg 2715, eq (13)--(16))
% omeg0 = sqrt(k/m)
% xid = C0/2/sqrt(k*m)
% nu = lamd*omeg0^r

% Parameters for impulse loading example (pg 2719, Fig 9)
omeg0 = 2*pi*5; % 5Hz
q = 1;
r = 0.6;
lambd = 0.3; % s^r (nu = lambd*omeg0^r works out to 2.374 as in the figure)
xid = 1.52;
m = 60000; %  kg - mass of the hammer (pg 2718)
k = omeg0^2*m; % stiffness (N/m)
C0 = 2*xid*sqrt(k*m);

% Loading - half sine force so that the final velocity is V0 = 0.4m/s
% period of half sine, T = 20ms, so that duration is 10ms
T = 0.02; % s
V0 = 0.4; % m/s
F0 = pi*m*V0/T; % N

% construct load fuction
tsample = 0.001; % s
tload = 0:tsample:T/2;
f = F0*sin(2*pi/T*tload);

% Pad zeros on both sides
f = [zeros(2048-length(f),1); f'; zeros(2048,1)];
t = (0:length(f)-1)'*tsample; % time for plotting

[x,p] = SolveViscodamperOscillatorFreqDom(m,k,lambd,C0,q,r,f,tsample);
figure(101),
    plot(t, x)