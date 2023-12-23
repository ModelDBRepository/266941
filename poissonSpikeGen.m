function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, dt)
% dt = 1/1000; % 1 ms
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;