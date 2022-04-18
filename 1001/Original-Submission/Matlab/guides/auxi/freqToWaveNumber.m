function k = freqToWaveNumber(frequencyHz)
cLight = 2.99792458E8; % light velocity in vacuum
k = 2*pi*frequencyHz/cLight; % wave number
end
