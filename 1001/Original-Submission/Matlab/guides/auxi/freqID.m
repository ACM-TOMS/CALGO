function [frequencyId,frequencyHz] = freqID(frequencies,frequencyHz)

frequencyId = find(frequencies == frequencyHz);
if numel(frequencyId) == 0
    frequencyId = numel(frequencies);
    frequencyHz = frequencies(frequencyId);
end

end
