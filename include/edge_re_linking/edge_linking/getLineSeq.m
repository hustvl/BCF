% written by ChengEn Lu
% This function gets sequence from the start point to the edpoint in a
% staight line
function [lineseq] = getLineSeq(startpoint, endpoint)
    startpt = round(startpoint);
    endpt = round(endpoint);
    subdst = endpt - startpt;
    [maxval, idx] = max(abs(subdst));
    if(maxval == 0)
        lineseq = startpoint;
        return;
    else
        dist1 = subdst(idx);
        longidx = idx;
        if(idx == 1)
            shortidx = 2;
            dist2 = subdst(2);
        else
            shortidx = 1;
            dist2 = subdst(1);
        end
        ra = dist2/dist1;
        retlist = zeros(maxval, 2);
        if(subdst(idx) < 0)
            step = -1;
        else
            step = 1;
        end
        i = startpoint(longidx) : step : endpoint(longidx);
        seq = round((i-startpoint(longidx)).*ra) + startpoint(shortidx);
        lineseq(:,longidx) = i;
        lineseq(:,shortidx) = seq;
    end
end