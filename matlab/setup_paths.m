function root = setup_paths()
    % SETUP_PATHS Add the example functions/tests without shadowing Coder helpers.
    % Run from any working directory: root = setup_paths(). Returns repository root.
    % Author: Reza Sameni | Emory University
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(here, 'functions'));
    addpath(fullfile(here, 'tests'));
    addpath(fullfile(here, 'examples'));
end
