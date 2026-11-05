function [roi, code, name] = getRegion(segVol, LUT, key)
    % key: numeric -> look up by code; string/char -> look up by name
    if isnumeric(key)
        row = find(LUT.code == key, 1);
    else
        row = find(string(LUT.name) == string(key), 1);
    end
    if isempty(row)
        error('getRegion:notFound', 'Region "%s" not found in LUT.', string(key));
    end
    code = LUT.code(row);
    name = char(LUT.name(row));

    roi = any(segVol == code, 4);     % works for 3D (1 layer) and layered volumes
    if ~any(roi, 'all')
        warning('getRegion:empty', 'ROI %s is empty!', name);
    end
end