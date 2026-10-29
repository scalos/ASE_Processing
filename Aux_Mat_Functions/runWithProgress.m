function runWithProgress(fcnHandle, subjName, subjIdx, nSubj)
    persistent totalTime compCount

    if isempty(totalTime) || subjIdx == 1
        totalTime = 0;
        compCount = 0;
    end

    t0 = tic;
    fcnHandle();
    dt = toc(t0);

    totalTime = totalTime + dt;
    compCount = compCount + 1;
    avgTime = totalTime / compCount;

    subjLeft = nSubj - subjIdx;
    estSec = avgTime * subjLeft;

    pctDone = round(100 * subjIdx / nSubj);

    fprintf('Current Subject: %s (%d/%d)\n', subjName, subjIdx, nSubj);
    fprintf('Progress: %d%% (%d/%d)\n', pctDone, subjIdx, nSubj);
    fprintf('Est. time remaining: %s\n\n', string(seconds(estSec),'hh:mm:ss'));
end