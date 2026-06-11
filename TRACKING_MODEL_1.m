%% Code to implement first tracking model.  
clear all
clc 
close all
rng(1) %fix rng 
tic
%% Neuronal Parameters
N = 400; %total number of neurons 
T = 900; %total simulation time (s)   
dt = 0.00005; %time step (s) 
nt = round(T/dt); 
tref = 0.002; %Refractory time constant in seconds 
tm = 0.01; %Membrane time constant 
vreset = -65; %Voltage reset 
vpeak = -40; %Voltage peak. 
td = 0.02; %decay time; 
tr = 0.002; %rise  time;  
nq = 20; %only store every nq time steps 
tdo = 0.8; %decay time of synthetic calcium 
tro = 0.02; %rise time of synthetic calcium  
%% Coupling Weight Matrix between CRH neurons 
OMEGA = zeros(N,N);
%% initialization parameters for the network 
IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
ro = r; %used for integrating calcium traces 
hro = hr; %used for integrating calcium traces,
JD = 0*IPSC; %storage variable required for integrating connections 
tspike = zeros(10*nt,2); %store spike times in tspike 
ns = 0; %Number of spikes, counts during simulation  
v = vreset + rand(N,1)*(10); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  
tlast = zeros(N,1); %variable sed to implement refractory period 
kd = 0; %used to initialize storage every nq time steps, kd = kd + 1 if mod(i,nq)==1; 
nq = 20; %store every nq time steps 
%REC2 = zeros(nt/nq,N); %store calcium traces in REC2
REC = zeros(nt,20); %store voltage traces in REC 
BIAS = vpeak-1; %background current
sig = 0.5; %white noise standard deviation  
idstore = randperm(N,20); %store a random subset of neurons 

%% set input weights for the burst/home-cage/novel environment signal  
WNV = randn(N,1);
WNV = sort(WNV);
NOV = 1; %amount of increase in current due to novel environment 
%% generate bursts during home-cage sessions 
tcrit1 = 300; %transfer from home-cage to novel environment 
tcrit2 = 600; %transfer from novel-environment back to home-cage 
dfs = 0.05; %approximate burst duration 
time = ((1:nt)*dt)';
dtimes = [tcrit1*rand(10,1);tcrit2 + (T-tcrit2)*rand(10,1)]; %generate random burst times 
dtimes = sort(dtimes); 
I_INT = zeros(nt,1); %burst storage variable 
for i = 1:length(dtimes) 
I_INT = I_INT + exp(-((time - dtimes(i)).^2)/dfs);     
end

%%  START INTEGRATION 
for i = 1:1:nt 

I = IPSC + BIAS + WNV*NOV*(dt*i>tcrit1)*(dt*i<tcrit2) + WNV*I_INT(i); %implement + home cage changes 
dv = (dt*i>tlast + tref).*( (-v+I)/tm)    ; %Voltage equation with refractory period 
v = v + dt*(dv) + sqrt(dt)*randn(N,1)*sig/sqrt(tm); 
index = find(v>=vpeak);  %Find the neurons that have spiked 


%Store spike times, and get the weight matrix column sum of spikers 
if length(index)>0
JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking  
tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i]; %store spike times 
ns = ns + length(index);  % total number of spikes so far
end

tlast = tlast + (dt*i -tlast).*(v>=vpeak);  %Used to set the refractory period of LIF neurons 

% Code if the rise time is 0, and if the rise time is positive 
if tr == 0  
    IPSC = IPSC*exp(-dt/td)+   JD*(length(index)>0)/(td);
    r = r *exp(-dt/td) + (v>=vpeak)/td;
else
    IPSC = IPSC*exp(-dt/tr) + h*dt;
h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);  %Integrate the current
r = r*exp(-dt/tr) + hr*dt; 
hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
%%integrate calcium resposnes
ro = ro*exp(-dt/tro) + hro*dt ; 
hro = hro*exp(-dt/tdo) + (v>=vpeak)/(tro*tdo);
end

v = v + (30 - v).*(v>=vpeak); %rest the voltage and apply a cosmetic spike.  
REC(i,:) = v(idstore); %store voltages 
v = v + (vreset - v).*(v>=vpeak); %reset spike time 




%if mod(i,nq)==1
%kd = kd + 1;
%REC2(kd,:) = ro; 
%end



%% plotting results 
    if mod(i,round(10/dt))==1
        prog = dt*i/T;
    end
    
end
%% 
% save tracking_model_1.mat -v7.3 

toc
