function demo
    libs = { ...
        'ImageMelding', ...
        fullfile('EdgeMap', 'models', 'forest'), ...
        fullfile('PatchMatch', 'dist'), ...
        fullfile('PiotrToolbox', 'channels'), ...
        fullfile('SelfTuningTexOpt', 'Utils'), ...
        fullfile('vlfeat', 'toolbox', 'mex', mexext), ...
    };
    % check that the lib folders exist (this does not mean the mex files do)
    for i = 1:length(libs)
        if ~exist(fullfile('..', libs{i}), 'dir')
            error('Missing a dependency: %s', libs{i});
        end
    end
    
    d = str2num(input('please input id range :', 's'));
    s = d(1);
    e = d(2);
    
    for i = s:e
        % remove out_directory and other variables
        str = num2str(i);
        input_file = ['exemplars' filesep str '.jpg'];
        output_dir = 'Results';
        if ~exist(input_file, 'file')
            continue;
        end   
        
        % set parameters  
        params.rand_seed            = 2; % random seed, set to number of 'shuffle'
        params.extra_channels       = 2; % number of extra channels to generate (0 = none, 1 = F, 2 = F+E) 
        params.hist_params          = 'default'; % histogram voting parameters
        params.nnf_weights          = 'auto'; % weights of patch channels in the nnf udpate
        params.alignment_channels   = [1, 7]; % L + F  //modified by huajie 2015-9-6
        params.nnf_channels         = [1, 2, 3, 8, 9, 10]; % L + a + b + F
        params.vote_channels        = [1, 2, 3]; % L + a + b modified by huajie
        params.vote_method          = 'histogram'; % early voting strategy 
        params.vote_method_until    = 5; % threshold between early and late
        params.vote_method_then     = 'default'; % late voting strategyb b  
        params.smart_init           = 1; % whether to use smart initialization
        params.smart_init_margin    = 25; % margin of alignment
        params.smart_init_trys      = 5; % number of tentatives for block alignment
        params.high_weight          = 3; % high weight value
        params.rand_search          = 10; % number of random search
        params.incomp_search        = 6; % number of low-completeness search
        params.comp_penalty         = 10; % completeness weight
        % call synthesis
        synth_func(input_file, output_dir, params);
        clear all;
    end
end