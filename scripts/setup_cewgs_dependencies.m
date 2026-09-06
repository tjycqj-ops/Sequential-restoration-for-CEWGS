function dependencies = setup_cewgs_dependencies(package_root)
%SETUP_CEWGS_DEPENDENCIES Locate YALMIP and the Gurobi MATLAB interface.
% Paths added here are temporary because the caller restores the MATLAB path.

validateattributes(package_root,{'char','string'},{'scalartext'});
package_root = char(package_root);

if exist('sdpvar','file')~=2
    yalmip_candidates = compact_candidates({ ...
        getenv('CEWGS_YALMIP_ROOT'), ...
        fullfile(package_root,'dependencies','YALMIP')});
    for candidate_index = 1:numel(yalmip_candidates)
        candidate = yalmip_candidates{candidate_index};
        if exist(fullfile(candidate,'@sdpvar','sdpvar.m'),'file')==2
            addpath(genpath(candidate));
            break
        end
    end
end

if isempty(which('gurobi'))
    gurobi_home = getenv('GUROBI_HOME');
    gurobi_candidates = compact_candidates({ ...
        getenv('CEWGS_GUROBI_MATLAB_DIR'), ...
        fullfile(gurobi_home,'matlab'), ...
        fullfile(gurobi_home,'win64','matlab'), ...
        fullfile(package_root,'dependencies','gurobi','win64','matlab')});
    for candidate_index = 1:numel(gurobi_candidates)
        candidate = gurobi_candidates{candidate_index};
        if contains_gurobi_interface(candidate)
            addpath(candidate);
            break
        end
    end
end

assert(exist('sdpvar','file')==2, ...
    ['YALMIP was not found. Set CEWGS_YALMIP_ROOT to the folder that ', ...
     'contains @sdpvar.']);
assert(~isempty(which('gurobi')), ...
    ['The Gurobi MATLAB interface was not found. Set ', ...
     'CEWGS_GUROBI_MATLAB_DIR to Gurobi''s win64\matlab folder.']);

dependencies = struct('yalmip',which('sdpvar'),'gurobi',which('gurobi'));
fprintf('[DEPENDENCIES] YALMIP: %s\n',dependencies.yalmip);
fprintf('[DEPENDENCIES] Gurobi: %s\n',dependencies.gurobi);
end

function candidates = compact_candidates(raw_candidates)
candidates = raw_candidates(~cellfun(@isempty,raw_candidates));
candidates = unique(candidates,'stable');
end

function available = contains_gurobi_interface(candidate)
available = exist(candidate,'dir')==7 && (...
    exist(fullfile(candidate,'gurobi.m'),'file')==2 || ...
    exist(fullfile(candidate,'gurobi.mexw64'),'file')==2);
end
