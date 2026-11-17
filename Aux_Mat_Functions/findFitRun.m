function [runDir, info] = findFitRun(fitsDir, tag)
% Picks the fit run folder to use for one subject.
%
%   fitsDir : the subject's fits folder (e.g. <subj>/session_ase/fits)
%   tag     : run label given to aseFit.sh --tag; "" = any run
%
% Run folders are named <yyyymmdd_HHMMSS> or <yyyymmdd_HHMMSS>_<tag> (aseFit.sh). With a tag,
% only folders ending in exactly "_<tag>" match; without one, every run folder matches. The
% latest match (by its timestamp) is returned. There is no fallback to an older run.
%
%   runDir : full path of the chosen run, "" if none matched
%   info   : struct with fields reason ("" when found), nMatches, candidates (all run names)
    arguments
        fitsDir {mustBeTextScalar}
        tag {mustBeTextScalar} = ""
    end
    fitsDir = string(fitsDir);  tag = string(tag);
    runDir = "";
    info = struct('reason',"", 'nMatches',0, 'candidates',strings(0,1));

    if ~isfolder(fitsDir)
        info.reason = "no fits folder (" + fitsDir + ")";
        return
    end
    d = dir(fitsDir);
    d = d([d.isdir] & ~startsWith({d.name}, '.'));
    names = string({d.name}).';
    tok = regexp(cellstr(names), '^(\d{8}_\d{6})((?:_.+)?)$', 'tokens', 'once');   % cell of cellstr
    isRun = ~cellfun(@isempty, tok);
    names = names(isRun);  tok = tok(isRun);
    info.candidates = names;
    if isempty(names)
        info.reason = "no run folders in " + fitsDir;
        return
    end
    stamps = strings(numel(tok),1);  tags = strings(numel(tok),1);
    for k = 1:numel(tok)
        stamps(k) = tok{k}{1};
        if numel(tok{k}) > 1 && ~isempty(tok{k}{2}), tags(k) = extractAfter(string(tok{k}{2}), 1); end   % drop the '_'
    end
    if strlength(tag) > 0
        match = tags == tag;
    else
        match = true(size(names));
    end
    info.nMatches = nnz(match);
    if ~any(match)
        info.reason = "no run tagged '" + tag + "' (runs: " + strjoin(names, ', ') + ")";
        return
    end
    cand = find(match);
    [~, order] = sort(stamps(cand) + "_" + names(cand));   % timestamp first; name breaks ties
    runDir = fullfile(fitsDir, names(cand(order(end))));
end
