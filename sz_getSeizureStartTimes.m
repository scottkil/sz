function seizureStartTimes = sz_getSeizureStartTimes(seizures)

% --- Keep only Type 1 events --- %
keepLog = strcmp({seizures.type},'1');
seizures = seizures(keepLog);


for szii = 1:numel(seizures)
    startt = seizures(szii).time(seizures(szii).trTimeInds(1));
    
    seizureStartTimes(szii) = startt;
end


end % function end