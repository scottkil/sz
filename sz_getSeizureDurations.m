function seizureDurations = sz_getSeizureDurations(seizures)

% --- Keep only Type 1 events --- %
keepLog = strcmp({seizures.type},'1');
seizures = seizures(keepLog);

for szii = 1:numel(seizures)
    startt = seizures(szii).time(seizures(szii).trTimeInds(1));
    endd = seizures(szii).time(seizures(szii).trTimeInds(end));
    seizureDurations(szii) = endd-startt;
end


end % function end