function seIDX = sz_FindEvents(sig,time,detThresh,shoulderThresh,ttc)

shoulderLog = sig > shoulderThresh;
pStarts(:,1) = find(diff(shoulderLog)>0)+1; % potential epoch starts
pEnds(:,1) = find(diff(shoulderLog)<0)+1; % potential epoch ends
if pEnds(1) < pStarts(1)
    pEnds(1) = [];
end
if pStarts(end) > pEnds(end)
    pStarts(end) = [];
end
%%
pseIDX = [pStarts,pEnds]; % potential  epochs

detLog = sig>detThresh;
dtIDX = find(detLog); % find indices of points beyond detection threshold
pseLab = false(size(pseIDX,1),1); % potential  epoch label
for psii = 1:size(pseIDX,1)
    currIDXs = pseIDX(psii,1):pseIDX(psii,2);
    if any(ismember(currIDXs,dtIDX))
        pseLab(psii) = true;
    end
end
%
nonnIDX = pseIDX(~pseLab,:);
nonnTime = time(nonnIDX);
pseTime = time(pseIDX(pseLab,:));
seIDX = pseIDX(pseLab,:);

while ~isempty(nonnTime)
    tfs = pseTime(:,1)-nonnTime(1,2); % is the epoch end close to an actual sleep epoch start?
    tkIDX = find(tfs>0); % check for true epochs that start after current epoch end
    if ~isempty(tkIDX) % if there are true epochs after the start of this one, continue
        subVec = tfs(tkIDX);
        [minVAL,minIDX] = min(subVec); % look for closest event
        if minVAL < ttc % if within the "time to combine" interval, then combine
            cIDX = tkIDX(minIDX);
            seIDX(cIDX,1) = nonnIDX(1,2);
        end
    end

    tfe = nonnTime(1,1) - pseTime(:,2); % is the epoch start close to an actual sleep epoch end?
    tkIDX = find(tfe>0);
    if ~isempty(tkIDX) % if there are true epochs after the start of this one, continue
        subVec = tfe(tkIDX);
        [minVAL,minIDX] = min(subVec); % look for closest event
        if minVAL < ttc % if within the "time to combine" interval, then combine
            cIDX = tkIDX(minIDX);
            seIDX(cIDX,2) = nonnIDX(1,1);
        end
    end


    nonnTime(1,:) = []; % remove this epoch from the list
    nonnIDX(1,:) = [];
end

%% === Merge events that happen close in time === %%
seTime = time(seIDX);
cfm = true; % check for merge flag
while cfm
    eints = seTime(2:end,1) - seTime(1:end-1,2); % intervals between epochs
    cIDX = find(eints<ttc,1,'first');
    if ~isempty(cIDX)
        seTime(cIDX,2) = seTime(cIDX+1,2); % merge
        seTime(cIDX+1,:) = []; % remove 2nd epochs from list
        seIDX(cIDX,2) = seIDX(cIDX+1,2); % merge
        seIDX(cIDX+1,:) = []; % remove 2nd epochs from list
    else
        cfm = false;
    end
end

end % function end