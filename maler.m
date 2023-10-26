%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  FUNCTION MALER: GET PROBLEM SPECIFICATION, SETUP GUI  %%%%%%%%%%%
%%%%%%%%                  EVALUATE PROBLEM NUMERICALLY AND PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*NASGU>  
%#ok<*CTCH>  
%#ok<*ASGLU> 
%#ok<*DEFNU>
%#ok<*AGROW>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maler(varargin)
%%% argumnets (optional) are:
%%% INPUT_FNAME - name of the input text-file,        {<name>,''},  string
%%% OUTPUT_DIR  - name of the folder to save figures, {<name>,''},  string
%%% MAKESYMFLAG - make symbolic calculations (if possible) {1, 0},  logical

    %%% get input argumnets
    INPUT_FNAME = '';  if ( nargin > 0 ),  INPUT_FNAME = varargin{1};  end,
    
    OUTPUT_DIR = 'out';     % default value   
    examples = {'example_txt','example_1d','example_2d'};
    for j = 1 : 1 : length(examples)
        if (~isempty(regexp(INPUT_FNAME,examples{j},'once')))
            OUTPUT_DIR = fullfile('ex',examples{j});    break,
        end
    end
    if ( ( nargin > 1 ) && (~isempty(varargin{2})) )
        OUTPUT_DIR = varargin{2};   
    end
    
    MAKESYMFLAG = f_check_toolbox('symbolic_toolbox');  % set default value
    if (MAKESYMFLAG)&&(nargin > 2), MAKESYMFLAG = logical(varargin{3}); end,
     
    % developer utility settings
    fUpdateFlag = false;             % default-settings file update control 
    inputDelims = {',',';','\s'};  % allowed keyboard input delimiters
    maxDimLimit = 2;             % limiter of maximal problem dimension
    maxTbarUBut = 7;           % maximal number of users' toolbar buttons
    gridParForm = int8(0);   % 0 ~ [xmin, xmax, length]; 1 ~ [xmin,dx,xmax]
    progVersion = 'Maler v6.04';
    makeIndlg   = false;     % create file 'maler_indlg.m'
    inFileType = '';         % {'txt','m','mat'}, undefined by now
 
    global MALER 
    MALER.DIM = [];                  % dimension of the problem
    MALER.WSP = struct();            % workspace variables
    MALER.GUI = struct();            % gui settings
    MALER.HDL = struct();            % gui handles
    MALER.PPR = struct();            % plot properties
    MALER.uset = struct();           % maler utility settings 
        
    MALER.WSP.Grid.name = {};        % names of the grid variables
    MALER.WSP.Grid.par  = {};        % [xmin, xmax, lth] or [xmin,dx,xmax]       
    MALER.WSP.Grid.val  = {};        % vectors (1d) and arrays (2d)
    
    MALER.WSP.Param.name = {};       % names of parameters       
    MALER.WSP.Param.val  = {};       % parameter values: [value,vmin,vmax]
    MALER.WSP.Param.dep  = {};       % list of depending functions
    
    MALER.WSP.Const.name = {};       % names of constants
    MALER.WSP.Const.val  = {};       % values of constants  
    
    MALER.WSP.Func.name = {};        % names of functions
    MALER.WSP.Func.expr = {};        % function formulas (string)
    MALER.WSP.Func.val  = {};        % numerical values at the grid
    MALER.WSP.Func.dep  = {};        % dependences
    MALER.WSP.Func.slc  = {};        % numerical values at the slice
    
    MALER.PPR.type = '';             % {'1d','2d','slc'}
    MALER.PPR.func = '';             % {'plot','contour',...}
    MALER.PPR.List.show = {};        % list of functions shown in listbox
    MALER.PPR.List.hide = {};        % list of hidden functions 
    MALER.PPR.Slc.lhside = {};       % slice expression: left-hand side
    MALER.PPR.Slc.rhside = {};       % slice expression: right-hand side
    MALER.PPR.linecolor  = {};       % color scheme for multy-line plots
    
    %%% create default settings folder if not exist    
    MALER.GUI.Defaults.dir  = fullfile(pwd,'auxiliary');   
    if (~exist(MALER.GUI.Defaults.dir,'dir'))
        [status,msg] = mkdir(MALER.GUI.Defaults.dir);
        if (~status),   disp(msg);
        else
            addpath(MALER.GUI.Defaults.dir);
            fUpdateFlag = true;
        end
    else
        addpath(MALER.GUI.Defaults.dir);
    end
    
    %%% setup gui, create default settings file if not exist
    MALER.GUI.Defaults.file = 'maler_defaults.in';
    if (~exist(fullfile(MALER.GUI.Defaults.dir,MALER.GUI.Defaults.file),'file'))
        f_gui_defaults_setup(true,maxTbarUBut);
    else
        f_gui_defaults_setup(fUpdateFlag,maxTbarUBut);
    end
    if ( isempty(MALER.GUI.Settings) )
        uiwait(errordlg({'Run terminated: gui setup failed.'}));
        disp([mfilename,': FATAL ERROR. Gui setup failed.']);
        return        
    end

    %%% create 'maler_indlg.m' if not exist
    f_create_indlg(MALER.GUI.Defaults.dir,makeIndlg);

    %%% create file 'maler_userbutt_callbacks.m' if not exist
    f_create_user_buttons_callbacks(fUpdateFlag,MALER.GUI.Defaults.dir,MALER.GUI.Settings.userbutts);
    
    %%% create output directory
    if (~exist(OUTPUT_DIR,'dir')) 
        [status,msg] = mkdir(OUTPUT_DIR);
        if (~status),         disp(msg); 
           uiwait(errordlg({['Error! Cannot create directory ',OUTPUT_DIR];...
                             'Use current folder for output'}));
           OUTPUT_DIR = '';
        else      disp(['Output directory is <',OUTPUT_DIR,'>']);
        end
    end
    
    %%% load input file if exist
    if (~isempty(INPUT_FNAME))  
        inFileType = f_load_input_file(INPUT_FNAME,MALER.GUI.Settings.fontsize,inFileType);    
    end
    
    %%% get <MALER> from workspace if exist, creat 'WORKSPACE' if not exist
    f_get_ws_maler();
    
    % store utility settings
    if (~isfield(MALER,'uset')),    MALER.uset = struct();      end,
    MALER.uset.progVersion = progVersion;
    MALER.uset.fUpdateFlag = fUpdateFlag;
    MALER.uset.inputDelims = inputDelims;
    MALER.uset.gridParForm = gridParForm;
    MALER.uset.maxTbarUBut = maxTbarUBut;
    MALER.uset.maxDimLimit = maxDimLimit;
    MALER.uset.numSpec = ['%.',num2str(MALER.GUI.Settings.dignumber),'g'];
    MALER.uset.saveFname = struct();    % save file name prefix & suffix
    MALER.uset.makeIndlg = makeIndlg;
    try
        MALER.uset.defGraphFsize = get(groot,'factoryUicontrolFontSize');
    catch
        try
            MALER.uset.defGraphFsize = get(0,'factoryUicontrolFontSize');
        catch
            MALER.uset.defGraphFsize = 8;
        end
    end
    MALER.uset.inFileType = inFileType;
    MALER.uset.matlabRelease = version('-release');
    MALER.uset.specifReleases = {'2013a','2015a'};
    MALER.uset.releaseCond = ( ismember(MALER.uset.matlabRelease,MALER.uset.specifReleases) );
    
    if ( nargin > 1 ),     MALER.uset.outputDir = OUTPUT_DIR;
    else
        if (isfield(MALER.uset,'outputDir')),  OUTPUT_DIR = MALER.uset.outputDir;  
        else                                   MALER.uset.outputDir = OUTPUT_DIR;
        end
    end
    if ( nargin > 2 ),     MALER.uset.makeSym = MAKESYMFLAG;
    else
        if (isfield(MALER.uset,'makeSym')),   MAKESYMFLAG = MALER.uset.makeSym;
        else                                  MALER.uset.makeSym = MAKESYMFLAG;
        end
    end
    
    %%% detect existence of independent variables and assign values
    if      (~isfield(MALER,'WSP'))     
        MALER.WSP.Grid.name = f_assign_vars('WORKSPACE.VARIABLES',MAKESYMFLAG);
    elseif  (~isfield(MALER.WSP,'Grid'))
        MALER.WSP.Grid.name = f_assign_vars('WORKSPACE.VARIABLES',MAKESYMFLAG);
    elseif  (~isfield(MALER.WSP.Grid,'name'))
        MALER.WSP.Grid.name = f_assign_vars('WORKSPACE.VARIABLES',MAKESYMFLAG);
    elseif ( isempty(MALER.WSP.Grid.name) )
        MALER.WSP.Grid.name = f_assign_vars('WORKSPACE.VARIABLES',MAKESYMFLAG);
    end
        
    %%% keyboar input of variables, create symbolic variables
    if (isempty(MALER.WSP.Grid.name)) 
        MALER.WSP.Grid.name = f_keyboard_input_vars(MALER.uset.defGraphFsize,MALER.GUI.Settings.fontsize,inputDelims,MAKESYMFLAG);
    end
    
    %%% remove non-symbolic variables from the list
    if (MAKESYMFLAG)
        MALER.WSP.Grid.name = f_remove_nontype(MALER.WSP.Grid.name,'base','sym','VARIABLES',MALER.GUI.Settings.fontsize);
    end
    
    % minimum problem dimension control
    MALER.DIM = numel(MALER.WSP.Grid.name);
    if ( MALER.DIM < 1 )
        h = errordlg({'Run terminated: empty list of variables'});
        h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize);  uiwait(h);
        disp([mfilename,': FATAL ERROR. Empty list of variables']);
        return
    end
    %%%--------------------------------------------------------------------
    
    %%% detect existence of parameters and assign names
    if      (~isfield(MALER,'WSP'))     
        MALER.WSP.Param.name = f_assign_vars('WORKSPACE.PARAMETERS',MAKESYMFLAG);
    elseif  (~isfield(MALER.WSP,'Param'))
        MALER.WSP.Param.name = f_assign_vars('WORKSPACE.PARAMETERS',MAKESYMFLAG);
    elseif  (~isfield(MALER.WSP.Param,'name'))
        MALER.WSP.Param.name = f_assign_vars('WORKSPACE.PARAMETERS',MAKESYMFLAG);
    elseif ( isempty(MALER.WSP.Param.name) )
        MALER.WSP.Param.name = f_assign_vars('WORKSPACE.PARAMETERS',MAKESYMFLAG);
    end

    %%% keyboar input of parameters, create symbolic parameters
    if (isempty(INPUT_FNAME)) 
        MALER.WSP.Param.name = f_keyboard_input_params(MALER.uset.defGraphFsize,MALER.GUI.Settings.fontsize,inputDelims,MAKESYMFLAG);
    end
    
    %%% remove non-symbolic parameters from the list
    if (MAKESYMFLAG)
        MALER.WSP.Param.name = f_remove_nontype(MALER.WSP.Param.name,'base','sym','PARAMETERS',MALER.GUI.Settings.fontsize);
    end
    
    %%% remove double-defined quantities from the parameter's list
    MALER.WSP.Param.name = setdiff(MALER.WSP.Param.name,MALER.WSP.Grid.name,'stable');
    
    % maximal problem dimension control
    if ( MALER.DIM > maxDimLimit )
        h = warndlg({'Too long list of variables';'treat extra variables as parameters'});
        h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize);  uiwait(h);
        for j = MALER.DIM : -1 : maxDimLimit+1
            MALER.WSP.Param.name = [MALER.WSP.Grid.name{j} , MALER.WSP.Param.name];
            MALER.WSP.Grid.name(j) = [];
        end
        MALER.DIM = numel(MALER.WSP.Grid.name);
    end 
    
    % transfer params and vars in base ws to the <WORKSPACE> structure
    wsp = evalin('base','WORKSPACE');   
    field = 'VARIABLES';    wsp.(field) = MALER.WSP.Grid.name; 
    field = 'PARAMETERS';   wsp.(field) = MALER.WSP.Param.name;  
    assignin('base','WORKSPACE',wsp);         
    %%%--------------------------------------------------------------------
     
    %%% keyboard input of functions 
    if (isempty(INPUT_FNAME))
        f_keyboard_input_funcs(MAKESYMFLAG);
    else
        f_eval_char_func('base',MAKESYMFLAG);
    end    
        
    %%% create cell arrays of structures, containing workspace funcs & const
    if      (~isfield(MALER,'WSP'))     
        [MALER.WSP.Func,MALER.WSP.Const] = f_get_ws_functions('base',union(MALER.WSP.Grid.name,MALER.WSP.Param.name));
    elseif  (~isfield(MALER.WSP,'Func'))
        [MALER.WSP.Func,MALER.WSP.Const] = f_get_ws_functions('base',union(MALER.WSP.Grid.name,MALER.WSP.Param.name));
    elseif  (~isfield(MALER.WSP.Func,'name'))
        [MALER.WSP.Func,MALER.WSP.Const] = f_get_ws_functions('base',union(MALER.WSP.Grid.name,MALER.WSP.Param.name));
    elseif ( isempty(MALER.WSP.Func.name) )
        [MALER.WSP.Func,MALER.WSP.Const] = f_get_ws_functions('base',union(MALER.WSP.Grid.name,MALER.WSP.Param.name));
    end  
    
    if (isempty(MALER.WSP.Func.name))
        h = errordlg({'Run terminated: empty list of functions'});
        h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize);  uiwait(h);
        disp([mfilename,': FATAL ERROR. Empty list of functions']);
        return
    end
    
    %%% synhronize base workspace
    f_sync_base_ws();
    
    %%% correct function expressions (atan-->atan2, operands)
    MALER.WSP.Func = f_correct_functions(MALER.WSP.Func);
    %%%--------------------------------------------------------------------    
      
    %%% get variable limits from the workspace or from keyboard input
    if     (~isfield(MALER,'WSP'))
        MALER.WSP.Grid.par = f_get_field('WORKSPACE','GRIDPARAM','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (~isfield(MALER.WSP,'Grid'))    
        MALER.WSP.Grid.par = f_get_field('WORKSPACE','GRIDPARAM','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (~isfield(MALER.WSP.Grid,'par'))
        MALER.WSP.Grid.par = f_get_field('WORKSPACE','GRIDPARAM','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (isempty(MALER.WSP.Grid.par))||(isempty(MALER.WSP.Grid.par{1}))
         MALER.WSP.Grid.par = f_get_field('WORKSPACE','GRIDPARAM','base',gridParForm,MALER.GUI.Settings.fontsize);
    end
    
    if (numel(MALER.WSP.Grid.par) ~= numel(MALER.WSP.Grid.name) )
        MALER.WSP.Grid.par = {};
        MALER.WSP.Grid.par = f_keyboard_input_gridparam(MALER.WSP.Grid.name,MALER.uset.defGraphFsize,MALER.GUI.Settings.fontsize,inputDelims,gridParForm);
    end
    wsp = evalin('base','WORKSPACE');   
    field = 'GRIDPARAM';    wsp.(field) = MALER.WSP.Grid.par;
    assignin('base','WORKSPACE',wsp);
    %%%--------------------------------------------------------------------
    
    %%% assign sampling numerical values to the parameters
    if     (~isfield(MALER,'WSP'))
        MALER.WSP.Param.val = f_get_field('WORKSPACE','PARSAMPVAL','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (~isfield(MALER.WSP,'Param'))    
        MALER.WSP.Param.val = f_get_field('WORKSPACE','PARSAMPVAL','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (~isfield(MALER.WSP.Param,'val'))
        MALER.WSP.Param.val = f_get_field('WORKSPACE','PARSAMPVAL','base',gridParForm,MALER.GUI.Settings.fontsize);
    elseif (numel(MALER.WSP.Param.val) ~= numel(MALER.WSP.Param.name) )
        MALER.WSP.Param.val = f_get_field('WORKSPACE','PARSAMPVAL','base',gridParForm,MALER.GUI.Settings.fontsize);    
    end
    
    if (numel(MALER.WSP.Param.val) ~= numel(MALER.WSP.Param.name) )
        MALER.WSP.Param.val = {};
        MALER.WSP.Param.val = f_keyboard_input_param_vals(MALER.WSP.Param.name,MALER.uset.defGraphFsize,MALER.GUI.Settings.fontsize);
    end
    wsp = evalin('base','WORKSPACE');   
    field = 'PARSAMPVAL';    wsp.(field) = MALER.WSP.Param.val;
    assignin('base','WORKSPACE',wsp);
    %%%--------------------------------------------------------------------
    
    %%% create computational grid
    if      (~isfield(MALER,'WSP'))
        MALER.WSP.Grid.val = f_create_grid(MALER.WSP.Grid.name,MALER.WSP.Grid.par,gridParForm,MALER.GUI.Settings.fontsize);
    elseif  (~isfield(MALER.WSP,'Grid'))
        MALER.WSP.Grid.val = f_create_grid(MALER.WSP.Grid.name,MALER.WSP.Grid.par,gridParForm,MALER.GUI.Settings.fontsize);
    elseif  (~isfield(MALER.WSP.Grid,'val'))
        MALER.WSP.Grid.val = f_create_grid(MALER.WSP.Grid.name,MALER.WSP.Grid.par,gridParForm,MALER.GUI.Settings.fontsize);
    elseif  ( isempty(MALER.WSP.Grid.val) )
        MALER.WSP.Grid.val = f_create_grid(MALER.WSP.Grid.name,MALER.WSP.Grid.par,gridParForm,MALER.GUI.Settings.fontsize);
    end        
        
    if ( isempty(MALER.WSP.Grid.val) )
        h = errordlg({'Run terminated: calculational grid is unset.'});
        h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize);  uiwait(h);
        disp([mfilename,': FATAL ERROR. Calculational grid is unset.']);
        return
    end   
    %%%--------------------------------------------------------------------
 
    %%% evaluate functions
    if      (~isfield(MALER,'WSP')),          status = f_eval_all_funcs();
    elseif  (~isfield(MALER.WSP,'Func')),     status = f_eval_all_funcs();
    elseif  (~isfield(MALER.WSP.Func,'val')), status = f_eval_all_funcs();
    elseif  ( isempty(MALER.WSP.Func.val) ),  status = f_eval_all_funcs();
    else                                      status = true;        
    end 
    if (~status)
        disp('Error! Cannot evaluate functions. Run terminated.');
        return,
    end
    if (isempty(MALER.WSP.Func.name))||(isempty(MALER.WSP.Func.name{1}))
        disp([mfilename,': Run terminated! Unable evaluate any function.']);
        return,
    end
    wsp = evalin('base','WORKSPACE');   
    field = 'FUNCTIONS';
    if (~isfield(wsp,field))||(~isfield(wsp.(field),'name'))||(numel(wsp.(field).name) ~= numel(MALER.WSP.Func.name))
        wsp.(field).name = MALER.WSP.Func.name;
        wsp.(field).expr = MALER.WSP.Func.expr;
        assignin('base','WORKSPACE',wsp);
    end
    %%%--------------------------------------------------------------------
    
    %%% find dependences
    f_find_func_dependences();      % find functional dependences
    f_find_dependlist_param();      % find denedence lists for parameters
    %%%--------------------------------------------------------------------
    
    %%% setup plotter
    if      (~isfield(MALER,'PPR')),                f_plotter_setup();
    elseif  (~isfield(MALER.PPR,'type')),           f_plotter_setup();
    elseif  ( isempty(MALER.PPR.type) ),            f_plotter_setup();
    end 
    %%%--------------------------------------------------------------------
    
    %%% setup line colors
    if  (~isfield(MALER.PPR,'linecolor'))
        MALER.PPR.linecolor = f_color_setup(length(MALER.WSP.Func.name)+1,MALER.GUI.Settings.dir);
    elseif (isempty(MALER.PPR.linecolor))
        MALER.PPR.linecolor = f_color_setup(length(MALER.WSP.Func.name)+1,MALER.GUI.Settings.dir);
    end
    %%%--------------------------------------------------------------------   
    
    %%% setup gui
    if      (~isfield(MALER,'GUI')),                f_gui_setup();
    elseif  (~isfield(MALER.GUI,'AxeI')),           f_gui_setup();
    elseif  ( isempty(MALER.GUI.AxeI) ),            f_gui_setup();
    end     
    %%%--------------------------------------------------------------------    
    
    %%% plot figure
    f_gui_draw();
    f_plot_selected();     
end    % end of the 'maler' function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   E N D   O F   T H E   < M A L E R >   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load or keyboard input the problem 

function installed = f_check_toolbox(tboxname)
%%% check if toolbox <tboxname> is installed
    installed = false;
    
    try     installed = logical(license('checkout',tboxname));
    catch
        if (strcmp(tboxname,'symbolic_toolbox'))
            try 
                evalin('base',      'sym q_w_e_r_t_y_g_b_v_c_x_z');
                evalin('base','clear sym q_w_e_r_t_y_g_b_v_c_x_z');
                installed = true;
            catch
            end
        end
    end
end     % end of the function <f_check_toolbox>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gui_defaults_setup(fupdateflag,maxbuttnum)
%%% setup default gui parameters
    global MALER
    defaultdir = MALER.GUI.Defaults.dir;    defaultfile = MALER.GUI.Defaults.file;

    MALER.GUI.Settings = struct();

    %%% create directory for default values if not exist
    if (exist(defaultdir,'dir')),     mkdir_stat = true; 
    else    [mkdir_stat, mkdir_msg] = mkdir(defaultdir);       end,
    if ( ~mkdir_stat )
        disp([mfilename,' FATAL ERROR: cannot create directory <',defaultdir,'>!']);
        disp(mkdir_msg);
        return
    end
    
    %%% get screen parameters
    try
        MALER.GUI.SCR.Units = 'pixels';
        set(0,'Units',MALER.GUI.SCR.Units);
        MALER.GUI.SCR.Size = get(0,'ScreenSize');       % monitor rezolution
        MALER.GUI.SCR.Width = MALER.GUI.SCR.Size(3);          % width
        MALER.GUI.SCR.Height = MALER.GUI.SCR.Size(4);         % height    
            MALER.GUI.SCR.Prop = int8(1);  % flag of sucsessful getting of properties
    catch  
        MALER.GUI.SCR.Prop = int8(0);  % flag of failure
        MALER.GUI.SCR.Width = 1600;
        MALER.GUI.SCR.Height = 1200;
    end

    MALER.GUI.Defaults.position = [0, 0.047, 1, 0.87];  % figure position
    
    MALER.GUI.Defaults.fontsize = 14;     % default fontsize
    if (MALER.GUI.SCR.Prop)
        MALER.GUI.Defaults.fontsize = max(8,...
            round(MALER.GUI.Defaults.fontsize * MALER.GUI.SCR.Height / 1050) );
    end
    
    MALER.GUI.Defaults.dignumber = 3;     % default number of digits behind comma
    
    MALER.GUI.Defaults.figsavetype = {'fig','pdf'};
    
    MALER.GUI.Defaults.userbutts = 7;     % default number of user buttons
    MALER.GUI.Defaults.userbutts = min(MALER.GUI.Defaults.userbutts,maxbuttnum);
    
    for j = 1 : 1 : MALER.GUI.Defaults.userbutts
        eval(['UserButtonsTooltips{',num2str(j),'} = ',...
              '[''User Defined Action '',num2str(j,''%02d'')];']);
        eval(['MALER.GUI.Defaults.tooltip',num2str(j,'%02d'),' = ',...
              'UserButtonsTooltips{',num2str(j),'};']);
    end

    %%% setup MALER.GUI
    MALER.GUI.Settings = MALER.GUI.Defaults;
    
    %%% update MALER.GUI default settings file if required 
    if (fupdateflag)
        f_gui_defaults_update(fupdateflag);
        return
    end

    %%% read MALER.GUI default settings from file, if it does exist
    fname = fullfile(defaultdir,defaultfile);

    if (exist(fname,'file')),      [fid,errmsg] = fopen(fname,'r');
        
        if (sign(fid) < 0)
            disp([mfilename,' ERROR: ',errmsg]);
        else
            while(~feof(fid))
                tline = fgetl(fid);    spaceind = isspace(tline);
                sepind = find(spaceind,1,'first');
                if ( isempty(sepind) ),     continue,       end,

                propname = tline(1:sepind-1);
                propval  = tline(sepind+1:end); 
                
                if (~isempty(regexp((lower(propname)),'tooltip','once')))
                    MALER.GUI.Settings.(lower(propname)) = strtrim(propval);
                else
                    MALER.GUI.Settings.(lower(propname)) = eval(propval);
                end               
            end
            fclose(fid);
        end
    end    
end     % end of the function <f_gui_defaults_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gui_defaults_update(fupdateflag)
%%% update default MALER.GUI settings and containing file
    global MALER
    MALER.GUI.Defaults = MALER.GUI.Settings;
    
    defaultdir = MALER.GUI.Defaults.dir;    defaultfile = MALER.GUI.Defaults.file;
    fname = fullfile(defaultdir,defaultfile);
    
    if (fupdateflag)
            
        [fid,errmsg] = fopen(fname,'w');       
        if ( sign(fid) < 0 )
            disp([mfilename,' ERROR: ',errmsg]);
            return
        end
        
        propnames = fieldnames(MALER.GUI.Defaults);
        for j = 1 : 1 : numel(propnames)
            
            if (strcmp(propnames{j},'dir' )),    continue,       end,
            
            if (strcmp(propnames{j},'file')),    continue,       end,
            
            if (strcmp(propnames{j},'position'))
                fprintf(fid,'%s\t[%g, %g, %g, %g]\r\n',upper(propnames{j}),...
                    MALER.GUI.Defaults.(propnames{j})(1),...
                    MALER.GUI.Defaults.(propnames{j})(2),...
                    MALER.GUI.Defaults.(propnames{j})(3),...
                    MALER.GUI.Defaults.(propnames{j})(4));
                
            elseif (strcmp(propnames{j},'figsavetype'))
                fprintf(fid,'%s\t{ ',upper(propnames{j}));
                for k = 1:1:numel(MALER.GUI.Defaults.(propnames{j}))
                    fprintf(fid,' ''%s'',',MALER.GUI.Defaults.(propnames{j}){k});
                end
                fseek(fid,-1,0);
                fprintf(fid,'  }\r\n');
                
            elseif (isnumeric(MALER.GUI.Defaults.(propnames{j})))
                fprintf(fid,'%s\t%g\r\n',upper(propnames{j}),MALER.GUI.Defaults.(propnames{j}));
            else
                fprintf(fid,'%s\t%s\r\n',upper(propnames{j}),MALER.GUI.Defaults.(propnames{j}));
            end
        end
        fclose(fid);
    end
end     % end of the function <f_gui_defaults_update>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_indlg(targetdir,flag) 
%%% create 'maler_indlg.m' with user defined font size if not exist

if (~flag),     return,     end,

if (   exist(  fullfile( targetdir , 'maler_indlg.m' ),'file'  )   )
    return
end

try
    %%% create file 'maler_indlg.m'
    rootdirectory  = fullfile(matlabroot,'toolbox','matlab','uitools');
    fullfilename   = fullfile(rootdirectory ,'inputdlg.m');
    newfilename    = fullfile(targetdir ,'maler_indlg.m' );
    copyfile(fullfilename,newfilename,'f');
    
    %%% change file 'maler_indlg.m' permission 
    [status, msg] = f_change_file_permission(newfilename);
    if ( status ~= 0)
        disp([mfilename,' WARNING: ',msg]);
    end

    %%% read file 'maler_indlg.m'  
    CharData = f_read_text_file(newfilename);
        
    %%% change the first and second strings
    CharData{1} = ['function Answer = maler_indlg(Prompt, Title, ',...
                   'NumLines, DefAns, Resize, UserFontSize)'];
    CharData{2} = ['if (nargin < 6), UserFontSize = ',...
                   'get(0,''FactoryUicontrolFontSize''); end,'];
    CharData{3} = 'UserFontSize = max(6,min(UserFontSize,20));  % limiter';

    %%% add function 'getnicedialoglocation'
    CharData{length(CharData)+1} = ' ';
    CharData{length(CharData)+1} = ['function figure_size = ',...
                   'getnicedialoglocation(figure_size, figure_units)'];
    CharData{length(CharData)+1} = 'parentHandle = gcbf;';
    CharData{length(CharData)+1} = 'propName = ''Position'';';
    CharData{length(CharData)+1} = ['if isempty(parentHandle),   ',...
             'parentHandle = 0;   propName = ''ScreenSize'';    end,'];
    CharData{length(CharData)+1} = 'old_u = get(parentHandle,''Units'');';
    CharData{length(CharData)+1} = 'set(parentHandle,''Units'',figure_units);';
    CharData{length(CharData)+1} = 'container_size=get(parentHandle,propName);';
    CharData{length(CharData)+1} = 'set(parentHandle,''Units'',old_u);';
    CharData{length(CharData)+1} = ['figure_size(1) = container_size(1)',...
                    '  + 1/2*(container_size(3) - figure_size(3));'];
    CharData{length(CharData)+1} = ['figure_size(2) = container_size(2)',...
                    '  + 2/3*(container_size(4) - figure_size(4));'];

    %%% add function 'setdefaultbutton'
    CharData{length(CharData)+1} = ' ';
    CharData{length(CharData)+1} = ['function ',...
                            'setdefaultbutton(figHandle, btnHandle)'];
    CharData{length(CharData)+1} = ['if nargin<1, errordlg(''MATLAB:',...
        'setdefaultbutton:InvalidNumberOfArguments'',',...
        '''Too few arguments for setdefaultbutton'');  end,'];
    CharData{length(CharData)+1} = ['if nargin>2, errordlg(''MATLAB:',...
        'setdefaultbutton:InvalidNumberOfArguments'',',...
        '''Too many arguments for setdefaultbutton''); end,'];
    CharData{length(CharData)+1} = 'if (usejava(''awt'') == 1)';
    CharData{length(CharData)+1} = ['    useJavaDefaultButton',...
                                        '(figHandle, btnHandle);'];
    CharData{length(CharData)+1} = 'else';
    CharData{length(CharData)+1} = ['    useHGDefaultButton',...
                                        '(figHandle, btnHandle);'];
    CharData{length(CharData)+1} = 'end';
    CharData{length(CharData)+1} = ' ';
    CharData{length(CharData)+1} = ['function useJavaDefaultButton',...
                                        '(figH, btnH)'];
    CharData{length(CharData)+1} = 'fh = handle(figH);';
    CharData{length(CharData)+1} = 'fh.setDefaultButton(btnH);';
    CharData{length(CharData)+1} = ' ';
    CharData{length(CharData)+1} = ['function useHGDefaultButton',...
                                        '(figHandle, btnHandle)'];
    CharData{length(CharData)+1} = 'btnPos = getpixelposition(btnHandle);';
    CharData{length(CharData)+1} = 'leftOffset   = btnPos(1) - 1;';
    CharData{length(CharData)+1} = 'bottomOffset = btnPos(2) - 2;';
    CharData{length(CharData)+1} = 'widthOffset  = btnPos(3) + 3;';
    CharData{length(CharData)+1} = 'heightOffset = btnPos(4) + 3;';
    CharData{length(CharData)+1} = ['h1 = uipanel(get(btnHandle, ',...
                '''Parent''), ''HighlightColor'', ''black'', ...'];
    CharData{length(CharData)+1} = ['''BorderType'', ''etchedout'',',...
                                    '''units'', ''pixels'', ...'];
    CharData{length(CharData)+1} = ['''Position'', [leftOffset ',...
                        'bottomOffset widthOffset heightOffset]);'];
    CharData{length(CharData)+1} = 'uistack(h1, ''bottom'');';

    %%% edit file strings
    for j = 1 : 1 : length(CharData)
        %%% change number of input arguments control
        if ( ~isempty( strfind(CharData{j},'0,5,nargin') ) )
            CharData{j} = 'errordlg(nargchk(0,6,nargin));';
        end
        %%% change resize definition string
        if ( ~isempty( strfind(CharData{j},'nargin==5') ) )
            CharData{j} = 'if nargin>4 && isstruct(Resize)';
        end
        %%% change font size entry in the 'CharData'
        if ( strncmp( 'TextInfo.FontSize',CharData{j},17 ) )
            CharData{j} = ['TextInfo.FontSize           = ',...
            'UserFontSize; % get(0,''FactoryUicontrolFontSize'');'];
            break
        end
   end
        
   %%% write corrected data in the file 'maler_indlg.m'
   [fileID, message] = fopen(newfilename,'w'); 
   if ( sign(fileID) < 0 )
       disp([mfilename,' ERROR: ',message]);
   else
       for j = 1 : 1 : length(CharData)
           fprintf(fileID,'%s\n', CharData{j});
       end
       fclose(fileID);
   end
                          
catch
    warndlg('WARNING: CANNOT CREATE ''INDLG.M''');
end

end   % end of the function <f_create_indlg> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, msg] = f_change_file_permission(fullfilename)
%%% change file permisiion to the full access

status = -1;        msg = '';
if ( ispc )
    cmnd = ['ICACLS "',fullfilename,'" /grant "',getenv('USERNAME'),'":F'];
elseif ( isunix )
    cmnd = ['chmod 777 ',fullfilename];
elseif ( ismac )
    cmnd = ['chmod ugo=rwx ',fullfilename];
end

try
    [status, msg] = system(cmnd);
catch
    warndlg([mfilename, ' WARNING: CANNOT CHANGE FILE <',fullfilename,'> PERMISSISONS']);
end 

end    % end of the function <f_change_file_permission> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_user_buttons_callbacks(fupdateflag,targetdir,N)
%%% create file 'maler_userbutt_callbacks.m' if not exist

    fname = fullfile(targetdir,'maler_userbutt_callbacks.m');
    if ( (fupdateflag)||(~exist(fname,'file')) )
    
        [fid,msg] = fopen(fname,'w');
        if (sign(fid) < 0)
            disp([mfilename,' ERROR: ',msg]);
            return
        end
        
        CharData{1} = 'function maler_userbutt_callbacks(varargin)';
        CharData = [CharData; '%%% callbacks for the user-defined buttons'];
        CharData = [CharData; ' '];
        CharData = [CharData; 'if ( numel(varargin) < 1 ),'];
        CharData = [CharData; '    h = errordlg(''ERROR! CANNOT GET BUTTON HANDLE!'');'];
        CharData = [CharData; '    h = f_resize_dlgbox(h,14);   uiwait(h),  return,'];
        CharData = [CharData; 'end'];
        CharData = [CharData; 'if ( numel(varargin) > 0 ),    hObject   = varargin{1};    end,'];
        CharData = [CharData; 'if ( numel(varargin) > 1 ),    eventdata = varargin{2};    end,'];
        CharData = [CharData; ' '];
        CharData = [CharData; 'set(hObject,''State'',''off'');'];
        CharData = [CharData; 'ButtonTag = get(hObject,''tag''); '];       
        CharData = [CharData; ' '];
        CharData = [CharData; 'switch ButtonTag   %%#ok<*CTCH> '];
        for j = 1:1:N
            CharData = [CharData; ' ']; 
            CharData = [CharData; ['    case ','''action',num2str(j,'%02d'),''''] ];
            CharData = [CharData; ['        % specify callback action for button ',num2str(j),' here'] ];
            CharData = [CharData; '        h = msgbox({'' '';''specify callback action in'';'' '';[''<auxiliary/'',mfilename,''.m>'']});'];
            CharData = [CharData; '        h = f_resize_dlgbox(h,14);   uiwait(h),  return,'];
        end
        CharData = [CharData; ' '];
        CharData = [CharData; '    otherwise'];
        CharData = [CharData; '        h = warndlg(''CALLBACK ACTION IS UNDEFINED'');'];
        CharData = [CharData; '        h = f_resize_dlgbox(h,14);   uiwait(h),  return,'];
        CharData = [CharData; 'end'];
        CharData = [CharData; 'end     % end of the ''maler_userbutt_callbacks'' function'];
        CharData = [CharData; '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'];
        CharData = [CharData; ' ']; 
        CharData = [CharData; 'function h = f_resize_dlgbox(h,fsize)'];
        CharData = [CharData; '%%% resize dialog boxes'];
        CharData = [CharData; '    fsize = max(6,min(fsize,20));         % limiter'];
        CharData = [CharData; '    htext = findobj(h, ''Type'', ''Text'');   % find text control in dialog box']; 
        CharData = [CharData; '    set(htext,''FontSize'',fsize);']; 
        CharData = [CharData; '    set(h,''Resize'',''on'');           pos = get(h,''Position'');'];
        CharData = [CharData; '    try'];
        CharData = [CharData; '        deffsize = get(0,''factoryUicontrolFontSize'');'];
        CharData = [CharData; '    catch   %#ok'];
        CharData = [CharData; '        deffsize = 8;'];
        CharData = [CharData; '    end'];
        CharData = [CharData; '    set(h,''Position'',[pos(1)-(pos(3)*(fsize/deffsize-1)/2), pos(2), pos(3)*fsize/deffsize, pos(4)]);'];
        CharData = [CharData; 'end    % end of the function <f_resize_dlgbox>'];        
        
        for j = 1:1:length(CharData)
            fprintf(fid,'%s\n', CharData{j});
        end
        fclose(fid);
    end
end     % end of the function <f_create_user_buttons_callbacks>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ftype = f_load_input_file(filename,fsize,ftype)
%%% load/evaluate file in a 'base' work space

    if (isempty(filename)),     return,     end,
    
    %%% try treat input file as mat-file
    try   
        evalin( 'base' , sprintf('load(''%s'')',filename) ); 
        ftype = 'mat';
        return
    catch
    end

    % check if extention = '.m', if not make copy with required extention
    [filepath,name,ext] = fileparts(filename);
    if ( (isempty(ext)) || (~strcmp(ext,'.m')) )       
        newname = fullfile(filepath,[name,'.m']); 
        ftype = 'txt';
        [status,msg,msgID] = copyfile(filename,newname);
    else
        newname = filename;
        status = true;
        ftype = 'm';
    end
    
    %%% try to evaluate input file (or copy) as m-file.
    if (status)     
        try 
            evalin( 'base' , sprintf('run(''%s'')',newname) );     
            if (~strcmp(newname,filename)),     delete(newname),    end,
            return,
        catch
            if (~strcmp(newname,filename)),     delete(newname),    end,
        end
    else
        disp([mfilename,': ERROR! cannot copy file ',filename]);
        disp(msg);       
    end
                    
    %%% if 'run' or 'copy' command failed, try to execute file string-by-string
    CharData = f_read_text_file(filename);
    if (isempty(CharData)),         return,         end,
        
    h = [];     delay = 0.1;
    for ii = 1:1:length(CharData)
        try
            evalin('base',CharData{ii});           
        catch
            errmsg = {' ';'ERROR: cannot evaluate string'; CharData{ii}};
            disp([mfilename,': ERROR! cannot evaluate string ',CharData{ii}]);               
            if (isempty(h))
                h = errordlg(errmsg);  h = f_resize_dlgbox(h,fsize);
                pos = get(h,'Position');    height = pos(4);
                pause(delay);
            else
                htext = findobj(h, 'Type', 'Text');
                str = get(htext,'String');
                str = [str ; errmsg];       set(htext,'String',str);
                pos = get(h,'Position');    
                pos(2) = pos(2)-height/2;   pos(4) = pos(4) + height;
                set(h,'Position',pos);      pause(delay);
            end               
        end
    end
    if (~isempty(h)),       uiwait(h);          end, 
    
end     % end of the function <f_load_input_file>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CharData = f_read_text_file(fname)
%%% read text file in a cell array [N,1]

    CharData = {};
    
    [fileID, message] = fopen(fname,'r');  
    if ( sign(fileID) < 0 )
        disp([mfilename,' FATAL ERROR: cannot open file ',fname]);
        disp(message);
        return
    end
    
    tline = 0;      j = 0;
    while ( tline ~= -1 )
        tline = fgetl(fileID);      % read one string
        if ( tline == -1 ),     break,        end,        % end of file
        if ( isempty(tline) ),  tline = ' ';  end,
        j = j + 1;             
        CharData{j,1} = tline;      % store the line 
    end
    fclose(fileID); 
end     % end of the function <f_read_text_file>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_eval_char_func(ws,makesym)
%%% evaluate char functions in 'ws' workspace
    if (makesym)
        allqnt = evalin(ws,'who');
        for j = 1:1:numel(allqnt)
            qnt = evalin(ws,allqnt{j});
            if (~ischar(qnt)),      continue,       end,
            try
                evalin(ws,[allqnt{j},' = ',qnt,';']);
            catch
                disp(['Error! Cannot evaluate function ',allqnt{j}]);
            end
        end
    end
end     % end of the function <f_get_ws_maler>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_get_ws_maler()
%%% get MALER from base workspace if exist, creat 'WORKSPACE' if not exist 
    global MALER
    if (  ismember( 'MALER',evalin('base','who') )  )
        maler = evalin('base','MALER');     
        mfiname = fieldnames(maler);
        Mfiname = fieldnames(MALER);
        finame = intersect(mfiname,Mfiname,'stable');
        for j = 1:1:length(finame)            
            MALER.(finame{j}) = maler.(finame{j});
        end
    end
    
    if ( ~ismember( 'WORKSPACE',evalin('base','who') )  )
        assignin('base','WORKSPACE',struct());
    end
end     % end of the function <f_get_ws_maler>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varnames = f_assign_vars(VarList,makesym)
%%% detect if variables of 'VarList' exist in a 'base' workspace and
%%% transfer them in the 'maler' workspace,
%%% otherwise create empty cell arrays in 'maler' and 'base' workspace
    varnames = {};     
    try
        VarList = evalin('base',VarList);

        for j = 1:1:length(VarList)
            var_name = VarList{j};     
            varnames = [varnames,var_name];
        end    
        %%% evaluate variables in the 'base' ws
        if (makesym)
            f_create_symvars(varnames,'base');
        end
    catch
        return
    end
end     % end of the function <f_assign_vars>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_symvars(qnt,ws)
%%% create sym vars in the workspace 
    allvars = evalin(ws,'who');
    for j = 1:1:length(qnt)
        if (~ismember(qnt{j},allvars))
            evalin(ws,[qnt{j},' = sym(''',qnt{j},''');'])
        end
    end
end     % end of the function <f_create_symvars_vars>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qnt = f_keyboard_input_vars(deffsize,fsize,delimiters,makesym)
%%% keyboard input of independent variables
    dlg_name   = 'ENTER SYMBOLIC VARIABLES';
    dlg_prompt = {'ENTER INDEPENDENT VARIABLES'};
    dlg_lines  = [1,round(length(dlg_prompt{1})*fsize/deffsize)];   
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {''};
    qnt = {};
    tline = [];
    while ( isempty(tline) )
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if (~isempty(answer))
            tline = answer{1};
            if (isempty(tline))
                h = errordlg({'ERROR: the list of independent variables cannot be empty'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);                              
                for j = 1:1:length(dlg_defans),   dlg_defans{j} = answer{j};   end,
            end
        else
            return,  
        end
    end
    
    % split input string to cell array of char such as {'1', '0.5', '3'}
    qnt = f_extract_keyboard_input(tline,delimiters);
    
    %%% assign values to <WORKSPACE.VARIABLES>
    if (~ismember('WORKSPACE',evalin('base','who')))
        assignin('base','WORKSPACE',struct());
    end
    
    wsp = evalin('base','WORKSPACE');   
    field = 'VARIABLES';    vname = ['WORKSPACE.',field];
    
    if (~isfield(wsp,field)),   evalin('base',[vname,' = {};']);     end,
    
    qnt = union(qnt,evalin('base',vname),'stable');      wsp.(field) = qnt;
    assignin('base','WORKSPACE',wsp);
    
    %%% evaluate variables in the 'base' ws
    if (makesym)
        f_create_symvars(qnt,'base');
    end    
end     % end of the function <f_keyboard_input_vars>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qnt = f_extract_keyboard_input(tline,delimiters)
%%% extract values from the input string
    qnt = {};
    for ii = 1 : 1 : length(delimiters)
        delim = delimiters{ii}; 
        if ( ~isempty(regexp(tline,delim,'once')) )
            inqnt = regexp(tline,delim,'split');
            for k = 1:1:length(inqnt)
                qnt = [ qnt ; strtrim(inqnt{k}) ];
            end
            break
        end
    end
    if ( isempty(qnt) ),    qnt = [ qnt ; strtrim(tline) ];      end,
end     % end of the function <f_extract_keyboard_input>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qntlist = f_remove_nontype(qntlist,ws,type,field,fsize)
%%% remove non-<type> expressions from the qntlist
%%% <ws> is the workspace
%%% <qntlist> is the cell array of the quantity names
%%% field is the proper field name of the base WORKSPACE structure
    S = evalin(ws,'who');
    nontype = {};
    for j = 1:1:length(qntlist)
        
        varname = qntlist{j};        
        if (  (~ismember(varname,S)) && (strcmp(type,'sym'))  )     
            f_create_symvars({varname},ws);     continue,        
        end
        
        var = evalin(ws,varname);                  
        varclass = getfield(whos('var'),'class');                %#ok<GFLD>
        if (  strcmp(varclass,type)  ),         continue,       end,
        
        nontype = [nontype,varname];
        h = warndlg({['Remove quantity ''',varname,''' form the list'];...
                        ['quantity type is non-',type]});
        h = f_resize_dlgbox(h,fsize);   uiwait(h);                                       
    end
    
    % transfer list of variables in the base workspace
    if (~isempty(nontype)),  
        qntlist = setdiff(qntlist,nontype);
        
        wsp = evalin('base','WORKSPACE');   
        vname = ['WORKSPACE.',field];
        if (~isfield(wsp,field)),   evalin('base',[vname,' = {};']);   end,
        
        wsp.(field) = qntlist;
        assignin('base','WORKSPACE',wsp);
    end
end    % end of the function <f_remove_nontype>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qnt = f_keyboard_input_params(deffsize,fsize,delimiters,makesym,varargin)
%%% keyboard input of parameters
    if (~isempty(varargin)),        stoper = varargin{1}; 
    else                            stoper = false;
    end
    dlg_name   = 'ENTER SYMBOLIC PARAMETERS';
    dlg_prompt = {'ENTER PARAMETERS'};
    dlg_lines  = [1,round(length(dlg_name)*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {''};
    qnt = {};
    tline = [];                                                 %%#ok<NASGU>
    try
        answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
    catch
        answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
    end
    if (~isempty(answer))
        tline = answer{1};
    else    return,  
    end
    if (isempty(tline)),     return,     end,
    
    % split input string
    for ii = 1 : 1 : length(delimiters)
        delim = delimiters{ii}; 
        if ( ~isempty(regexp(tline,delim,'once')) )
            inqnt = regexp(tline,delim,'split');
            for j = 1:1:length(inqnt)
                qnt = [ qnt ; strtrim(inqnt{j}) ];
            end
            break
        end
    end
    if ( isempty(qnt) ),    qnt = [ qnt ; tline ];      end,
    
    %%% interraption
    if (stoper),    return,     end,
    
    %%% assign values to <WORKSPACE.PARAMETERS>
    if (~ismember('WORKSPACE',evalin('base','who')))
        assignin('base','WORKSPACE',struct());
    end
    
    wsp = evalin('base','WORKSPACE');   
    field = 'PARAMETERS';    vname = ['WORKSPACE.',field];
    
    if (~isfield(wsp,field)),   evalin('base',[vname,' = {};']);     end,
    
    qnt = union(qnt,evalin('base',vname),'stable');      wsp.(field) = qnt;
    assignin('base','WORKSPACE',wsp);
    
    %%% evaluate variables in the 'base' ws
    if (makesym)
        f_create_symvars(qnt,'base');
    end
end     % end of the function <f_keyboard_input_params>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = f_keyboard_input_funcs(makesym)
%%% keyboard input of functions, one-by-one
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    inputlist = {};
    
    dlg_name = 'TEMPLATE: F = ...';
    dlg_prompt = {'ENTER NEW FUNCTION'};
    dlg_lines  = [1,round(length(dlg_prompt{1})*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {''};
    
    while 1
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if (isempty(answer)),       break,       end,
        if (isempty(answer{1})),    break,       end,
        
        tline = answer{1};          inqnt = regexp(tline,'=','split');
        if ( length(inqnt) < 2 )           
            h = errordlg({'The equality sign is lost! try again'}); 
            h = f_resize_dlgbox(h,fsize);     uiwait(h);          continue,
        elseif ( length(inqnt) == 2 )   
            fname = strtrim(inqnt{1}); 
            fstr  = strtrim(inqnt{2});
        else
            h = errordlg({'Too many equality signs! try again'}); 
            h = f_resize_dlgbox(h,fsize);     uiwait(h);          continue,
        end
        
        if ( ismember(fname,evalin('base','WORKSPACE.VARIABLES')) )
            h = warndlg({'The notation is occupired by variable! try again'}); 
            h = f_resize_dlgbox(h,fsize);     uiwait(h);          continue,            
        end
        
        if ( ismember(fname,evalin('base','WORKSPACE.PARAMETERS')) )
            h = warndlg({'The notation is occupired by parameter! try again'}); 
            h = f_resize_dlgbox(h,fsize);     uiwait(h);          continue,            
        end
        
        if ( ismember(fname,evalin('base','who')) || (ismember(fname,MALER.WSP.Func.name)) )
            if (~isempty(fstr))
                choice = questdlg(['Redefine existing function ',fname,'?'],...
                'Update function definition', 'Yes','No','Yes'); 
                if (strcmp(choice,'No'))
                    if (ismember(fname,MALER.WSP.Func.name))
                        continue,
                    else
                        try
                            fstr = char(evalin('base',fname));
                        catch
                            h = errordlg({['Error! Cannot evaluate function ',fname]}); 
                            h = f_resize_dlgbox(h,fsize);     uiwait(h);      continue,                
                        end
                    end
                end
            else
                choice = questdlg(['Remove existing function ',fname,'?'],...
                'Clear function definition', 'Yes','No','No'); 
                if (strcmp(choice,'No')),       continue,       end,            
            end
        end
        
        % remove function from the list or add new function        
        if (isempty(fstr))
            if (ismember(fname,MALER.WSP.Func.name))
                [C,ia] = setdiff(MALER.WSP.Func.name,fname,'stable');            
                MALER.WSP.Func.expr = MALER.WSP.Func.expr(ia);
                if (isfield(MALER.WSP.Func,'val'))&&(~isempty(MALER.WSP.Func.val))
                    MALER.WSP.Func.val = MALER.WSP.Func.val(ia);
                end
                if (isfield(MALER.WSP.Func,'dep'))&&(~isempty(MALER.WSP.Func.dep))
                    MALER.WSP.Func.dep = MALER.WSP.Func.dep(ia);
                end               
                MALER.WSP.Func.name = C;
                evalin('base',['clear ',fname]);
            elseif (ismember(fname,evalin('base','who')))
                evalin('base',['clear ',fname]);
            else
                h = errordlg({['Error! Cannot evaluate function ',fname]}); 
                h = f_resize_dlgbox(h,fsize);     uiwait(h);      continue,                
            end
            
        else 
            
            str = fstr;
            if (makesym)
                if (~strcmp(fstr(end),';')),  str = [fstr,';'];    end,
                try                     
                    evalin('base',[fname,' = ',str]);
                    fstr = char(evalin('base',fname));
                catch                         
                    h = errordlg({['Error! Cannot evaluate function ',fname]}); 
                    h = f_resize_dlgbox(h,fsize);   uiwait(h);    continue, 
                end
            else
                try
                    evalin('base',[fname,' = ','''',str,'''',';']);
                    fstr = char(evalin('base',fname));
                catch
                    h = errordlg({['Error! Cannot evaluate function ',fname]}); 
                    h = f_resize_dlgbox(h,fsize);   uiwait(h);    continue,                    
                end
            end
            
            % replace function expression or add new function
            if (ismember(fname,MALER.WSP.Func.name))
                [C,ia,ib] = intersect(MALER.WSP.Func.name,fname,'stable'); %%#ok<ASGLU>
                MALER.WSP.Func.expr{ia} = fstr;
            else
                MALER.WSP.Func.name{end+1} = fname; 
                MALER.WSP.Func.expr{end+1} = fstr;               
            end           
            inputlist = [inputlist, fname];
        end        
    end
    
    MALER.PPR.List.hide = intersect(MALER.WSP.Func.name,MALER.PPR.List.hide,'stable');
    MALER.PPR.List.show = setdiff(MALER.WSP.Func.name,MALER.PPR.List.hide,'stable');
    varargout{1} = inputlist;
    
    %%% update 'WORKSPACE'
    wsp = evalin('base','WORKSPACE');  field = 'FUNCTIONS';       
    wsp.(field).name = MALER.WSP.Func.name; 
    wsp.(field).expr = MALER.WSP.Func.expr;
    assignin('base','WORKSPACE',wsp);
end    % end of the function <f_keyboard_input_func>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Funcs,Consts] = f_get_ws_functions(ws,skiplist)
%%% get list of functions and constants defined in the workspace <ws>
    Funcs = struct();           Consts = struct();     
    Funcs.name = {};            Consts.name = {};
    Funcs.expr = {};            Consts.val  = {};
    S = evalin(ws,'whos');
    for j = 1:1:length(S)
        var = S(j);     varname = var.name;     varclass = var.class;
        if (ismember(varname,skiplist)),    continue,       end,
        
        if (strcmp(varclass,'sym'))
            fnc.name = varname;
            fnc.expr = char(evalin(ws,varname));
            Funcs.name = [Funcs.name,fnc.name];
            Funcs.expr = [Funcs.expr,fnc.expr];
            continue,
        end
        
        if (strcmp(varclass,'char'))
            fnc.name = varname;
            fnc.expr = evalin(ws,fnc.name);
            Funcs.name = [Funcs.name,fnc.name];
            Funcs.expr = [Funcs.expr,fnc.expr];            
        end
        
        if ((strcmp(varclass,'double'))&&(length(evalin(ws,varname)) == 1))
            const.name = varname;
            const.val  = evalin(ws,varname);
            Consts.name = [Consts.name,const.name];
            Consts.val  = [Consts.val ,const.val ];
            continue,
        end
    end
 
    wsp = evalin('base','WORKSPACE');  
    field = 'FUNCTIONS'; wsp.(field) = Funcs;  assignin('base','WORKSPACE',wsp);
    field = 'CONSTANTS'; wsp.(field) = Consts; assignin('base','WORKSPACE',wsp);   
end    % end of the function <f_get_ws_functions>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_sync_base_ws()
%%% synhronize base workspace with maler workspace
    global MALER
    makesym = MALER.uset.makeSym;
    basewsqnt = evalin('base','who');
    
    if (ismember('WORKSPACE',basewsqnt))
        wsp = evalin('base','WORKSPACE');
        if (isfield(wsp,'FUNCTIONS'))         
            Funcs = wsp.FUNCTIONS;
        else           
            if ((~isfield(MALER.WSP.Func,'name'))||(isempty(MALER.WSP.Func.name)))  
                return 
            end            
            Funcs = MALER.WSP.Func;
        end
    end   
    
    for j = 1:1:length(Funcs.name)
        if (ismember(Funcs.name{j},basewsqnt))
            continue,
        else
            if (makesym)
                try  
                    evalin('base',[Funcs.name{j},' = ',Funcs.expr{j},';']);
                catch
                    disp(['Error! Cannot evaluate function ',Funcs.name{j}]);
                end
            else
                try
                    evalin('base',[Funcs.name{j},' = ','''',Funcs.expr{j},'''',';']);
                catch
                    disp(['Error! Cannot evaluate function ',Funcs.name{j}]);
                end                    
            end
        end
    end
end    % end of the function <f_sync_base_ws>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Funcs = f_correct_functions(Funcs)
%%% replace atan by atan2 and operands by array operands

    for j = 1:1:length(Funcs.name)
        fname = Funcs.name{j};      instr = Funcs.expr{j};
        [runstat,outstr] = f_replace_atan(instr);
        if (sign(runstat) < 0)
            disp([mfilename,': Warning! Unable to replace <atan> by <atan2> in ',fname]);
            outstr = instr;
        elseif (~strcmp(instr,outstr))
            disp([mfilename,' Message: <atan> is replaced by <atan2> in ',fname]);
        end
        funcstr = f_matrix_operands(outstr);
        Funcs.expr{j} = funcstr;
    end
    
end    % end of the function <f_correct_functions>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runstat,outstr] = f_replace_atan(instr, varargin)
%%% replace atan by atan2 with different angle ranges
    if (isempty(varargin)),     AngleRange = [];            % [-pi, pi]   
    else                        AngleRange = varargin{1};   % [0, 2*pi]+2*pi*AngleRange
    end
    runstat = 0; 
    
    instr  = instr(~isspace(instr));
    outstr = instr;
    
    i0  = regexp(outstr,'atan(');    % indexes of the atan function entries
    rid = [];                        % replaced pattern indexes
    
    while (~isempty(i0))
                
        ind0 = max( setdiff( i0 , rid ) );
        if (isempty(ind0)),   break,     end,
        
        rid = [rid, ind0];      % replaced entry index
        
        substr = outstr(ind0:end);      counter = 0;
        atgI0  = ind0 + 5;              % first index of atan argument
        for j = 5 : 1 : length(substr)
            if ( strcmp(substr(j),'(') ),    counter = counter + 1;    end, 
            if ( strcmp(substr(j),')') ),    counter = counter - 1;    end,
            if ( counter == 0 )
                atgI1 = ind0 + j - 2;   % last  index of atan argument
                break,      
            end
        end
        atgarg = outstr(atgI0 : atgI1);    % argument of the atan function
        
        if (isempty(regexp(atgarg,'/','once'))),    continue,       end,
        
        substr = atgarg;                counter = 0;
        for j = length(substr) : -1 : 2
            if ( substr(j) == ')' ),    counter = counter - 1;      end,
            if ( substr(j) == '(' ),    counter = counter + 1;      end,
            if ( counter == 0 ) && ( substr(j-1) == '/' ),  
                denomI0 = j;        break,  
            end
            if (j == 2)
                disp([mfilename,': FATAL ERROR!'])
                runstat = -1;
                return
            end
        end
        denom = substr(denomI0:end);          % denominator
        substr(denomI0:end) = [];
        if (isempty(substr)),    continue,       end,
        
        substr(end) = [];
        if (isempty(substr)),    continue,       end,
        
        if  ( strcmp(substr(end),'.') ),     substr(end) = [];   end,
        if (isempty(substr)),    continue,       end,
        
        atg2arg = [substr,',',denom];        
        substr = outstr(ind0:atgI1+1);      % from 'a' to ')'
        
        if (isempty(AngleRange))      % [-pi,pi]
            
            newstr = ['atan2(',atg2arg,')'];
            
        else  % [0, 2*pi] + 2*pi*AngleRange  
            
            newstr = ['mod(atan2(',atg2arg,')+2*pi,2*pi)',...
                '+',num2str(round(AngleRange(2))),'*2*pi'];
        end
               
        outstr = strrep(outstr,substr,newstr);
        i0 = regexp(outstr,'atan('); % indexes of the atan function entries
    end
end    % end of the function <f_replace_atan>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function funcstr = f_matrix_operands(funcstr)
%%% replace scalar operands by array/matrix operands
    funcstr = strrep(funcstr, '/' , './');
    funcstr = strrep(funcstr, '*' , '.*');
    funcstr = strrep(funcstr, '^' , '.^');
    
    funcstr = strrep(funcstr, '..', '.' );
    
    %%% matrix operands are supposed to be entered as double operands
    funcstr = strrep(funcstr, '././', '/' );
    funcstr = strrep(funcstr, '.*.*' , '*');
    funcstr = strrep(funcstr, '.^.^' , '^');
end   % end of the function <f_matrix_operands>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = f_get_field(strname,field,ws,gridParForm,fsize)
%%% get value 'field' of the structure 'strname' in 'ws' workspace
    vals = {};
    if (ismember(strname,evalin(ws,'who')))
        wsp = evalin(ws,strname);
        if (isfield(wsp,field))
            vals = wsp.(field);
        end
    end
    
    % convert to the form [xmin , dx, xmax];
    if (  (~isempty(vals)) && (strcmp(field,'GRIDPARAM')) && (sign(gridParForm) > 0)  )
        for j = 1:1:numel(vals)
            try
                x = linspace(vals{j}(1), vals{j}(2), vals{j}(3));
                vals{j} = [ min(x) , x(2)-x(1) , max(x) ];
            catch
                h = errordlg({['ERROR: inconsistent grid parameters of variable ',num2str(j)];'try keyboard input'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
                vals = {};      return,
            end
        end
        
    elseif (  (~isempty(vals)) && (strcmp(field,'PARSAMPVAL')) && (numel(vals{1}) ~= 3)  )
        for j = 1:1:numel(vals)
            vals{j} = [vals{j}(1), min(vals{j}(1)/4, vals{j}(1)-1), max(vals{j}(1)*4, vals{j}(1)+1)]; 
        end
    end    
end     % end of the function <f_get_field>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gridparam = f_keyboard_input_gridparam(vars,deffsize,fsize,delimiters,gridParForm)
%%% keyboard input of the grid parameters      
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {'','',''};
    
    gridparam = {};
    for j = 1:1:length(vars)
        dlg_name   = ['variable ',vars{j}];
        dlg_prompt = {'enter Min, Max, Length'};
        dlg_lines  = [1,round(length(dlg_prompt{1})*fsize/deffsize)];
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        
        if ( (isempty(answer))||(isempty(answer{1})) )
            h = warndlg({'Empty input. Use default values: 0, 1, 101'});
            h = f_resize_dlgbox(h,fsize);   uiwait(h);
            tline = '0,1,101';
        else
            tline = answer{1};
        end
        
        % split input string
        qnt = [];
        for ii = 1 : 1 : length(delimiters)
            delim = delimiters{ii}; 
            if ( ~isempty(regexp(tline,delim,'once')) )
                inqnt = regexp(tline,delim,'split');
                for k = 1:1:length(inqnt)
                    try     
                        val = eval(inqnt{k});
                    catch
                        val = 0; 
                        h = warndlg({['Error: Unable evaluate ',inqnt{k}]});
                        h = f_resize_dlgbox(h,fsize);   uiwait(h);
                    end
                    qnt = [ qnt , val ];
                end
                break
            end
        end
        if ( isempty(qnt) )    
            try    
                val = eval(tline);
            catch
                val = 0; 
                h = warndlg({['Error: Unable evaluate ',tline]});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
            end
            qnt = val;      
        end
        
        switch length(qnt)
            case 0
                h = warndlg({'Empty input. Use default values: 0, 1, 101'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
                qnt = [0,1,101];
            case 1
                h = warndlg({'Right boundary and Length are not defined. Use default values'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
                qnt = [qnt(1),qnt(1)+1,101];
            case 2
                h = warndlg({'Length is not defined. Use default value 101'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
                qnt = [qnt,101];
            case 3
                % OK, do nothing
            otherwise
                h = warndlg({'Excessive input, spit extra values'});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);
                qnt(4:end) = [];
        end
        % input correction
        qnt = [min(qnt(1),qnt(2)), max(qnt(1),qnt(2)), max(1,round(abs(qnt(3)))) ];
        if ( abs(qnt(2)-qnt(1)) < 10*eps ) 
            h = warndlg({'Left and Right boundaries and the same. Use default values'});
            h = f_resize_dlgbox(h,fsize);   uiwait(h);
            qnt(2) = qnt(2)+1;
        end
        gridparam{j} = qnt;
    end
    
    wsp = evalin('base','WORKSPACE');       field = 'GRIDPARAM';
    wsp.(field) = gridparam;    assignin('base','WORKSPACE',wsp);
    
    % convert to the form [xmin , dx, xmax];
    if (  (~isempty(gridparam)) && (sign(gridParForm) > 0)  )
        for j = 1:1:numel(gridparam)
            x = linspace(gridparam{j}(1), gridparam{j}(2), gridparam{j}(3));
            gridparam{j} = [ min(x) , x(2)-x(1) , max(x) ];          
        end
    end
end     % end of the function <f_keyboard_input_gridparam>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramvals = f_keyboard_input_param_vals(params,deffsize,fsize)
%%% assign values to numerical parameters
    paramvals = {};
    if ( (isempty(params)) || (isempty(params{1})) ),    return,       end,
    
    dlg_name   = 'ASSIGN TRIAL VALUE';
    dlg_lines  = [1,round(length(dlg_name)*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {''};
    
    for j = 1:1:length(params)
        dlg_prompt = {['parameter ',params{j}]};
        
        while 1
            try
                answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
            catch
                answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
            end            
            if (isempty(answer))        % 'cancel' button pressed 
                choice = questdlg({['Warning: Parameter <',params{j},'> is unset!'];'Assign value?'}, ...
                'Unassigned parameter', 'Yes','No', 'Yes');                
                if ( strcmp(choice,'Yes') ),      continue,
                else                              qnt = [];      break,
                end                            
            end
            
            if (isempty(answer{1}))      % 'ok' button pressed                
                choice = questdlg({['Warning: Parameter <',params{j},'> is unset!'];'Assign value?'}, ...
                'Unassigned parameter', 'Zero','NaN','Other', 'Zero');                
                if     ( strcmp(choice,'Zero') ),     qnt = 0;       break,
                elseif ( strcmp(choice,'NaN' ) ),     qnt = NaN;     break,
                else    continue, 
                end                            
            end
            
            try     qnt = eval(answer{1});      break,
            catch
                h = warndlg({'Error! Cannot evaluate answer. Try again'});
                h = f_resize_dlgbox(h,fsize);     uiwait(h);     continue,
            end
        end
        paramvals{j} = [qnt , min(qnt/4,qnt-1) , max(qnt*4,qnt+1)];
    end
    
    wsp = evalin('base','WORKSPACE');       field = 'PARSAMPVAL';
    wsp.(field) = paramvals;    assignin('base','WORKSPACE',wsp);
end     % end of the function <f_keyboard_input_param_vals>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gridvars = f_create_grid(vars,gridparams,gridParForm,fsize)
%%% create computational grid
    gridvars = {};                  % cell array of grid points     
    DIM = length(vars);             % number of independent variables
    
    switch int8(DIM)
        
        case 0
            h = warndlg({'Error! execution terminated';'list of independent variables is empty'});
            h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
            
        case 1
            [a,b,c] = f_matsplit(gridparams{1}); 
            if (sign(gridParForm) < 1),     x = linspace(a,b,c);
            else                            x = a : b : c ;
            end
            
        case 2
            [a,b,c] = f_matsplit(gridparams{1});    
            if (sign(gridParForm) < 1),     x = linspace(a,b,c);
            else                            x = a : b : c ;
            end
            
            [a,b,c] = f_matsplit(gridparams{2});    
            if (sign(gridParForm) < 1),     y = linspace(a,b,c);
            else                            y = a : b : c ;
            end
            
            [X , Y] = meshgrid( x , y );            % [dimy,dimx]
            
        otherwise
            h = warndlg({'Error! execution terminated';'too long list of independent variables'});
            h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
    end
    
    if (exist('x','var')),      gridvars{1,1} = x;      end,
    if (exist('y','var')),      gridvars{2,1} = y;      end,
    
    if (exist('X','var')),      gridvars{1,2} = X;      end,
    if (exist('Y','var')),      gridvars{2,2} = Y;      end,
end     % end of the function <f_create_grid>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = f_matsplit(x,axis)
% return matrix elements as separate output arguments
% optionally can specify an axis to split along.
% example: [a1,a2,a3,a4] = split(1:4)
% example: [x,y,z] = split(zeros(10,3),2)
    if ( nargin < 2 ),       axis = 2;          end,    % default: columns
    dims = num2cell(size(x));       
    dims{axis}= ones([1 dims{axis}]);
    varargout = mat2cell(x,dims{:});
end     % end of the function <f_matsplit>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = f_eval_all_funcs(varargin)
%%% evaluate all functions
    global MALER
    varargout{1} = 0;       % execution status
    
    if ((isempty(MALER.WSP.Func.name)) || (isempty(MALER.WSP.Func.name{1})))  
        return,   
    end
    
    if ((isempty(MALER.WSP.Grid.name) ) || (isempty(MALER.WSP.Grid.name{1} )))            
        return,   
    end
    
    if (~isempty(varargin))    
        f_u_u_n_k_I_n_D = varargin{1};
        if (  size(f_u_u_n_k_I_n_D,1) > size(f_u_u_n_k_I_n_D,2)  )
            f_u_u_n_k_I_n_D = f_u_u_n_k_I_n_D.';
        end
    else
        f_u_u_n_k_I_n_D = 1:1:numel(MALER.WSP.Func.name);
    end      
    
    % create max-dimensional numerical grid 
    for y_n_d_e_k_s = 1:1:length(MALER.WSP.Grid.name)
        eval( [MALER.WSP.Grid.name{y_n_d_e_k_s} , ' = MALER.WSP.Grid.val{y_n_d_e_k_s , end};'] );    
    end 
    
    %%% create numerical parameters
    if ( ~isempty(MALER.WSP.Param.name) )&&( ~isempty(MALER.WSP.Param.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Param.name)
            try
                eval( [MALER.WSP.Param.name{y_n_d_e_k_s},' = MALER.WSP.Param.val{y_n_d_e_k_s}(1);'] );
            catch
                disp(['Error! Cannot evaluate parameter ',MALER.WSP.Param.name{y_n_d_e_k_s}]);
                h = (errordlg(['Error! Cannot evaluate parameter ',MALER.WSP.Param.name{y_n_d_e_k_s}]));
                h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize); uiwait(h);                
            end
        end
    end
    
    %%% create constants
    if ( ~isempty(MALER.WSP.Const.name) )&&( ~isempty(MALER.WSP.Const.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Const.name)
            eval( [MALER.WSP.Const.name{y_n_d_e_k_s},' = MALER.WSP.Const.val{y_n_d_e_k_s};'] );
        end
    end    
    
    %%% evaluate functions
    f_u_u_n_k_s = MALER.WSP.Func;             % struct
    f_e_i_l_e_t = []; 
    for y_n_d_e_k_s = f_u_u_n_k_I_n_D 
        f_u_u_n_k_c.name = f_u_u_n_k_s.name{y_n_d_e_k_s};
        f_u_u_n_k_c.expr = f_u_u_n_k_s.expr{y_n_d_e_k_s};
        if (~strcmpi(f_u_u_n_k_c.expr(end),';')),  f_u_u_n_k_c.expr = [f_u_u_n_k_c.expr,';'];  end,
        try     
             f_u_u_n_k_v_a_l_u = eval(f_u_u_n_k_c.expr);
             if (~isnumeric(f_u_u_n_k_v_a_l_u(1))),   f_u_u_n_k_v_a_l_u = eval(f_u_u_n_k_v_a_l_u);    end,
             if (numel(f_u_u_n_k_v_a_l_u) < 2),       f_u_u_n_k_v_a_l_u = f_u_u_n_k_v_a_l_u*ones(size(MALER.WSP.Grid.val{1,end}));    end,            
        catch
            e_r_o_r_h_e_n_d_e_l = errordlg({['ERROR: Unable evaluate function ',f_u_u_n_k_c.name]});
            e_r_o_r_h_e_n_d_e_l = f_resize_dlgbox(e_r_o_r_h_e_n_d_e_l,MALER.GUI.Settings.fontsize);  
            uiwait(e_r_o_r_h_e_n_d_e_l);
            f_u_u_n_k_v_a_l_u = NaN;
            f_e_i_l_e_t = [f_e_i_l_e_t , y_n_d_e_k_s];
        end
        
%         if (strcmpi(f_u_u_n_k_c.name,'theta'))
%             f_u_u_n_k_v_a_l_u = f_u_u_n_k_v_a_l_u + (f_u_u_n_k_v_a_l_u < 0)*2*pi;
%         end        
        
        MALER.WSP.Func.val{y_n_d_e_k_s} = f_u_u_n_k_v_a_l_u;
        if (~MALER.uset.makeSym)            
            eval([MALER.WSP.Func.name{y_n_d_e_k_s},' = MALER.WSP.Func.val{y_n_d_e_k_s};']);            
        end
    end
    
    %%% remove functions, that where not calculated, from the list
    if (~isempty(f_e_i_l_e_t))
        f_y_l_t_n_a_m_s = fieldnames(MALER.WSP.Func);
        for y_n_d_e_k_s = 1 : 1 : length(f_y_l_t_n_a_m_s);
            f_y_l_d_t = f_y_l_t_n_a_m_s{y_n_d_e_k_s};
            f_y_l_d_t_v_a_l_u = MALER.WSP.Func.(f_y_l_d_t);
            if (~isempty(f_y_l_d_t_v_a_l_u))
                f_y_l_d_t_v_a_l_u(f_e_i_l_e_t) = [];
                MALER.WSP.Func.(f_y_l_d_t) = f_y_l_d_t_v_a_l_u;
            end
        end
    end    
    
    %%% check if function list is not empty
    if (isempty(MALER.WSP.Func.name))||(isempty(MALER.WSP.Func.name{1}))
        e_r_o_r_h_e_n_d_e_l = warndlg({'Warning: List of functions is empty!'});
        e_r_o_r_h_e_n_d_e_l = f_resize_dlgbox(e_r_o_r_h_e_n_d_e_l,MALER.GUI.Settings.fontsize);  
        uiwait(e_r_o_r_h_e_n_d_e_l);
        MALER.PPR.List.show = {};
        MALER.PPR.List.hide = {};
        f_update_listbox_funcs();
        return
    else
        %%% update PPR.List
        if ((isfield(MALER,'PPR'))&&(isfield(MALER.PPR,'List'))&&(isfield(MALER.PPR.List,'show'))&&(~isempty(MALER.PPR.List.show)))
            MALER.PPR.List.show = intersect(MALER.PPR.List.show,MALER.WSP.Func.name,'stable');        
        else
            MALER.PPR.List.show = MALER.WSP.Func.name;               
        end
        MALER.PPR.List.hide = setdiff(MALER.WSP.Func.name,MALER.PPR.List.show,'stable');        
    end
    
    %%% update listbox <show functions>
    f_update_listbox_funcs();
     
    %%% evaluate slice
    if ((isfield(MALER.PPR.Slc,'val')) && (~isempty(MALER.PPR.Slc.val{1})) && (~isempty(MALER.PPR.Slc.val{2})))
        f_eval_slice_funcs(f_u_u_n_k_I_n_D);
    end
    
    varargout{1} = 1;       % execution status
end    % end of the function <f_eval_all_funcs>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_find_func_dependences(varargin)
%%% find dependences between functions and parameters
    global MALER
    
    % get indexes of functions to find dependeces
    if (~isempty(varargin)),    IND = varargin{1};
    else                        IND = 1:1:numel(MALER.WSP.Func.name);
    end
    if (isempty(IND)),          return,         end,
    
    % delete too large indexes
    for j = length(IND) : -1 : 1
        if ( IND(j)) > numel(MALER.WSP.Func.name ),   IND(j) = [];    end,
    end
    if (isempty(IND)),          return,         end,
    
    % get list of functions to find dependences
    F = struct();
    F.name = MALER.WSP.Func.name(IND);    % func names to find dependences
    F.expr = MALER.WSP.Func.expr(IND);    % func expr to find dependences
    F.dep  = cell(size(F.name));          % dependences
    
    % create the list of all known quantities
    if ( (isfield(MALER.WSP,'Grid')) && (isfield(MALER.WSP.Grid,'name')) )
            varnames = MALER.WSP.Grid.name;
    else    varnames = {};
    end
    if ( (isfield(MALER.WSP,'Const')) && (isfield(MALER.WSP.Const,'name')) )
            constnames = MALER.WSP.Const.name;
    else    constnames = {};
    end
    if ( (isfield(MALER.WSP,'Param')) && (isfield(MALER.WSP.Param,'name')) )
            parnames = MALER.WSP.Param.name;
    else    parnames = {};
    end   
    funcnames = MALER.WSP.Func.name;
    
    qntlist = union(union(union(varnames,constnames),parnames),funcnames);
    
    before = {'','\W'};     % before = {'' , '(' , '+' , '-' , '*' , '/' , '\^'};
    after = before;         % after  = {'' , ')' , '+' , '-' , '*' , '/' , '\^' , '\.' , ',' , ';'};
    
    i1 = 1:1:length(before);
    i2 = 1:1:length(after);
    
    [I1,I2] = meshgrid(i1,i2);
    II = cat(2,I1.',I2.');
    II = reshape(II,[],2);      II(1,:) = [];   % remove case {'',''}     
    
    for j = 1:1:length(F.name)
        
        expr = F.expr{j};   
        expr = expr( ~isspace(expr) );      % remove spaces
        if (  (~strcmp(expr(end),';')) && (~strcmp(expr(end),','))  )
            expr = [expr,';'];
        end
        dep = {''};
        
        for k = 1:1:length(qntlist)
            for n = 1:1:size(II,1)
                patt = [before{II(n,1)} , qntlist{k} , after{II(n,2)}];
                if (~isempty(regexp(expr,patt,'once')))                   
                    dep = union(dep,qntlist(k));
                    break,
                end
            end
        end
        dep(1) = [];        % remove entry ''   
        F.dep{j} = dep;     
    end
    
    if (~isfield(MALER.WSP.Func,'dep')),    
        MALER.WSP.Func.dep = cell(size(MALER.WSP.Func.name));
    end
    MALER.WSP.Func.dep(IND) = F.dep;
    
end    % end of the function <f_find_func_dependences>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_find_dependlist_param(varargin)
%%% find denedence lists for parameters
    global MALER
    
    % get indexes of parameters to create dependence list
    if (~isempty(varargin)),    IND = varargin{1};
    else                        IND = 1:1:numel(MALER.WSP.Param.name);
    end
    if (isempty(IND)),          return,         end,
    
    % delete too large indexes
    for j = length(IND) : -1 : 1
        if ( IND(j)) > numel(MALER.WSP.Param.name ),   IND(j) = [];    end,
    end
    if (isempty(IND)),          return,         end,
    
    if (~isfield(MALER.WSP.Param,'dep'))
        MALER.WSP.Param.dep = cell(size(MALER.WSP.Param.name));
    end
    
    for j = IND
        parname = MALER.WSP.Param.name{j};      % name of parameter
        dep = {''};
        for k = 1:1:numel(MALER.WSP.Func.name)
            if (ismember(parname,MALER.WSP.Func.dep{k}))
                dep = union(dep,MALER.WSP.Func.name(k));
            end
        end
        dep(1) = [];        % remove entry '' 
        MALER.WSP.Param.dep{j} = dep;
    end    

end    % end of the function <f_find_dependlist_param>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_update_listbox_funcs()
%%% update listbox <select functionsto plot>
    global MALER
    
    MALER.GUI.Listbox.funcs.String = MALER.PPR.List.show;
    if (isfield(MALER.GUI.Listbox.funcs,'Value'))
        if (~isempty(MALER.GUI.Listbox.funcs.Value))
            for j = length(MALER.GUI.Listbox.funcs.Value) : -1 : 1
                if (MALER.GUI.Listbox.funcs.Value(j)) > length(MALER.PPR.List.show)
                    MALER.GUI.Listbox.funcs.Value(j) = [];
                else          break,
                end
            end
        else
            MALER.GUI.Listbox.funcs.Value = 1; 
        end
    else
        MALER.GUI.Listbox.funcs.Value = 1;
    end
    if (isempty(MALER.GUI.Listbox.funcs.Value))
        MALER.GUI.Listbox.funcs.Value = 1;
    end
    
    if ((isfield(MALER,'HDL'))&&(isfield(MALER.HDL,'Listbox'))&&(isfield(MALER.HDL.Listbox,'funcs')))
        set(MALER.HDL.Listbox.funcs,'String', MALER.GUI.Listbox.funcs.String);
        set(MALER.HDL.Listbox.funcs,'Value' , MALER.GUI.Listbox.funcs.Value);
    end
end     % end of the function <f_update_listbox_funcs>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotter_setup()
%%% setup plotter properties
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    switch MALER.DIM
        case 0            
            disp([mfilename,': FATAL ERROR! The problem dimension is zero!']);
            h = errordlg({'Run terminated. The problem dimension is zero!'});   
            h = f_resize_dlgbox(h,fsize);  uiwait(h);   return,
            
        case 1
            MALER.PPR.type = '1d';                  % {'1d','2d','slc'}
            MALER.PPR.func = 'plot';
            MALER.PPR.absciss = MALER.WSP.Grid.name{1};
            
        case 2    
            MALER.PPR.type = '2d';                  % {'1d','2d','slc'}                  
            MALER.PPR.func = 'contourf';
            MALER.PPR.absciss = '';
            
        otherwise
            disp([mfilename,': FATAL ERROR! Too high problem dimension!']);
            h = errordlg({'Run terminated. Too high problem dimension!'});   
            h = f_resize_dlgbox(h,fsize);  uiwait(h);   return,            
    end
end     % end of the function <f_plotter_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function col = f_color_setup(N,varargin)
%%% setup color scheme for N-curve plot
    N = int8( abs(N) );     col = cell(N,N);   
    if ( N < 1 ),           return,         end,
    
    folder = '';    infname = '';   LimVal = -1;   
    if (~isempty(varargin)),    folder = varargin{1};       end,
    % list the folder content, look for the file 'maler_colspec_*'
    listing = dir(folder);
    for j = 1:1:length(listing),   entry = listing(j).name;
        if (~isempty(regexp(entry,'maler_colspec','once')))
            infname = entry;       break,              
        end
    end
        
    if (~isempty(infname))     
        % try to get limiting number of lines from the input file name
        res = regexp(infname , '(?<limval>\d+)' , 'names');
        if (~isempty(res)),      LimVal = eval(res(end).limval);      end,
    end

    if ( (~isempty(infname)) && (exist(infname,'file')) && (N < LimVal+1) )
        % read color specifiers from the text file
        try
            A = cell(0,0);
            fid = fopen(infname,'r');  a = fscanf(fid,'%f');   fclose(fid);

            while (~isempty(a))   
                j = size(A,1) + 1;
                for k = 1:1:j
                    A{j,k}(1) = a(1);  A{j,k}(2) = a(2);  A{j,k}(3) = a(3);
                    a(1:3) = [];
                end   
            end
            col = A;
        catch
           disp([mfilename,': Error. Cannot read file ',infname]);
           col = {};
        end
    end
        
    if ((isempty(col)) || (isempty(col{1}))) % create color specifiers here
        for j = 1:1:N
            x = [1,2];          y = rand( numel(x), j);
            fh = figure('Visible','off');     hold on,   ph = plot(x,y);
            for k = 1:1:j,         col{j,k} = get(ph(k),'Color');      end,
            close(fh),
        end
    end
    
    %%% create color specification file
    if (isempty(infname))
        N = 25;     
        fname = ['maler_colspec_',num2str(N),'.txt'];
        f_create_colspec_file(N,fname);
    end
end     % end of the function <f_color_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_colspec_file(N,fname)
%%% create file with color specifications
    global MALER
    if (~exist(MALER.GUI.Settings.dir,'dir')),      return,     end,
    
    fname = fullfile(MALER.GUI.Settings.dir,fname);
    ColSpec = cell(N,N);
    
    for j = 1:1:N
        col = cell(1,j);           
        x = [1,2];          
        y = rand( numel(x), j);
        fh = figure('Visible','off');     hold on,   ph = plot(x,y);

        for k = 1:1:j,   col{1,k} = get(ph(k),'Color');   end,   close(fh),
        ColSpec(j,1:j) = col;
    end

    fid = fopen(fname,'w');  if (sign(fid)<0),     return,     end,   
    for j = 1:1:size(ColSpec,1)
        for k = 1:1:j
            fprintf(fid,'%.3f %.3f %.3f ',ColSpec{j,k}(1),ColSpec{j,k}(2),ColSpec{j,k}(3));
        end
        fprintf(fid,'\r\n');
    end   
    fclose(fid);

end     % end of the function <f_create_colspec_file>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gui_draw()
%%% create figure, axes object and upper toolbar pushbutton
    global MALER 

    %%% create figure for plotting profiles
    obj = 'Fig0';   
    MALER.HDL.(obj) = figure();                 f_create_gui_object(obj);
    
    %%% create empty upper buttons
    MALER.HDL.UserToolbar = findall(MALER.HDL.Fig0,'Type','uitoolbar');
    for j = 1 : 1 : 20
        MALER.HDL.deadbutt(j) = uitoggletool( MALER.HDL.UserToolbar,'Enable','Off',...
            'Visible','Off','Separator','On','HandleVisibility','Off');
    end
    
    %%% create active toolbar buttons including user-defined ones
    tbarbuttnames = fieldnames(MALER.GUI.TbarButt);  flag = true;
    for j = 1 : 1 : length(tbarbuttnames) 
        
        Name = tbarbuttnames{j};
        
        if (  (isempty(regexp(Name,'UserAct','once'))) && (flag)  )
            MALER.HDL.TbarButt.(Name) = uitoggletool(MALER.HDL.UserToolbar,...
                'State','off','Separator','On','HandleVisibility','On');        
            obj = ['TbarButt.',Name];       f_create_gui_object(obj);
            
        elseif (  (~isempty(regexp(Name,'UserAct','once'))) && (flag)  )
            %%% create empty upper buttons
            k0 = length(MALER.HDL.deadbutt);
            for k = k0+1 : 1 : k0+10
                MALER.HDL.deadbutt(k) = uitoggletool( MALER.HDL.UserToolbar,'Enable','Off',...
                    'Visible','Off','Separator','On','HandleVisibility','Off');
            end           
            flag = false;
            
            MALER.HDL.TbarButt.(Name) = uitoggletool(MALER.HDL.UserToolbar,...
                'State','off','Separator','On','HandleVisibility','On');          
            obj = ['TbarButt.',Name];       f_create_gui_object(obj);
            
        else
             MALER.HDL.TbarButt.(Name) = uitoggletool(MALER.HDL.UserToolbar,...
                'State','off','Separator','On','HandleVisibility','On');          
            obj = ['TbarButt.',Name];       f_create_gui_object(obj);           
        end
    end
    
    %%% create button <Apply>
    MALER.HDL.FigButt.Apply = uicontrol(MALER.HDL.Fig0,'Style','pushbutton');
    obj = 'FigButt.Apply';              f_create_gui_object(obj);
    
    %%% create containing axes object
    obj = 'Axe0';   
    MALER.HDL.(obj) = axes('Parent',MALER.HDL.Fig0);   f_create_gui_object(obj);
        
    %%% create visible axes object
    obj = 'AxeI';   
    MALER.HDL.(obj) = axes('Parent',MALER.HDL.Fig0);    f_create_gui_object(obj);
            
    %%% create grid variables sliders
    f_create_grid_sliders();
    
    %%% create button <choose absciss>
    MALER.HDL.FigButt.ChAbs = uicontrol(MALER.HDL.Fig0,'Style','pushbutton');
    obj = 'FigButt.ChAbs';          f_create_gui_object(obj);
    
    %%% create button <set grid parameters>
    MALER.HDL.FigButt.SetGrid = uicontrol(MALER.HDL.Fig0,'Style','pushbutton');
    obj = 'FigButt.SetGrid';        f_create_gui_object(obj); 
    
    %%% create button <set parameter values>
    MALER.HDL.FigButt.SetParam = uicontrol(MALER.HDL.Fig0,'Style','pushbutton');
    obj = 'FigButt.SetParam';       f_create_gui_object(obj);
    
    %%% create listbox <choose parameters to tune> and its label
    MALER.HDL.Listbox.param = uicontrol(MALER.HDL.Fig0,'Style','listbox');
    obj = 'Listbox.param';          f_create_gui_object(obj);
    MALER.HDL.TextParam = uicontrol('Parent',MALER.HDL.Fig0,'Style','text');
    obj = 'TextParam';              f_create_gui_object(obj);
    
    %%% create listbox <choose functions to plot> and its label
    MALER.HDL.Listbox.funcs = uicontrol(MALER.HDL.Fig0,'Style','listbox');
    obj = 'Listbox.funcs';          f_create_gui_object(obj);
    MALER.HDL.TextFuncs = uicontrol('Parent',MALER.HDL.Fig0,'Style','text');
    obj = 'TextFuncs';              f_create_gui_object(obj);
    
    %%% create listbox <hide/show functions>
    MALER.HDL.Listbox.hideshow = uicontrol(MALER.HDL.Fig0,'Style','listbox');
    obj = 'Listbox.hideshow';       f_create_gui_object(obj);
    
    %%% create parameter sliders and their labels
    f_create_param_sliders(); 
        
    %%% create listbox <choose absciss>
    MALER.HDL.Listbox.absciss = uicontrol(MALER.HDL.Fig0,'Style','listbox');
    obj = 'Listbox.absciss';               f_create_gui_object(obj);
    
    MALER.GUI.Exist = true;   
end    % end of the 'f_gui_draw' function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plot_selected()
%%% plot selected functions
    global MALER 
    fsize = MALER.GUI.Settings.fontsize;
    
    %%% remove old plots
    if (isfield(MALER.HDL,'Plot'))
        f_close_figobj(MALER.HDL.Plot);         
        MALER.HDL = rmfield(MALER.HDL,'Plot');
    end      
    
    % define plot-function
    if (strcmp(MALER.PPR.func,'cross-section'))
        MALER.PPR.func = 'plot';
    end
    pfunc = eval(['@',MALER.PPR.func,';']);             
    
    fselnumbs = get(MALER.HDL.Listbox.funcs,'Value');
    allfnames = get(MALER.HDL.Listbox.funcs,'String');
    
    if (numel(allfnames) < 1),              return,     end,
    if (fselnumbs > length(allfnames)),     return,     end,    
    fselnames = allfnames(fselnumbs);             % names of selected funcs
    
    % get color specifiers
    if ( length(fselnames) > size(MALER.PPR.linecolor,1) ) 
        MALER.PPR.linecolor = f_color_setup( length(allfnames) + 1 );
    end       
       
    Funcs = MALER.WSP.Func;     % all functions
    
    if (isempty(findall(MALER.HDL.Fig0,'type','axes','Tag','AxeI')))
        obj = 'AxeI';   MALER.HDL.(obj) = axes('Parent',MALER.HDL.Fig0);    
        f_create_gui_object(obj);
    end
    axes(MALER.HDL.AxeI);       % make AxeI active
    
    switch MALER.PPR.type
        
        case '1d' ,     f_plot_1d(pfunc,fselnames,Funcs);
            
        case 'slc',     f_plot_slice(pfunc,fselnames);
                        
        case '2d' ,
            
            x = MALER.WSP.Grid.val{1,2};    % 2d-array   
            y = MALER.WSP.Grid.val{2,2};    % 2d-array
            
            switch MALER.PPR.func
                
                case {'contour','contourf','mesh','meshc'} % separate plots
                    f_plot_2d_cont_mesh(x,y,fselnames,Funcs);
                    
                case {'surf','surfc'}                      % plot in 1 axes
                    f_plot_2d_surf(pfunc,x,y,fselnames,Funcs);
                    
                otherwise
                    disp([mfilename,': ERROR! ','Unknown plot function: ',MALER.PPR.func]);
                    h = errordlg(['Unknown plot function: ',MALER.PPR.func]);
                    h = f_resize_dlgbox(h,fsize);   uiwait(h);      return,                                                  
            end
                        
        otherwise
            disp([mfilename,': ERROR! ','Unknown plot type: ',MALER.PPR.type]);
            h = errordlg(['Unknown plot type: ',MALER.PPR.type]);
            h = f_resize_dlgbox(h,fsize);   uiwait(h);      return,            
    end
end    % end of the 'f_plot_selected' function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n1,n2,ax] = f_subplot(N,spec,fsize,fh)
%%% create N subplots
%%% N - total number of subplots
%%% spec - axes handel, or 'Position' vector in normalized units
%%% fsize - axes labels fontsize, fh - parent figure handle

    N = abs(N);   n1 = [];    n2 = [];    ax = [];  
    
    if ( N < 1 ),       return,     end,
  
    n2 = floor(sqrt(N));    % number of columns
    n1 = ceil(N/n2);        % number of rows
    
    if (numel(spec)==1),    Pos0 = get(spec,'Position');
    else                    Pos0 = spec;
    end
    
    hshift = Pos0(3)/10;    % horizontal shift between subplots
    vshift = Pos0(4)/10;    % vertical shift between subplots
    try
        units = get(fh,'Units');       set(fh,'Units','points');
        figpos = get(fh,'Position');   set(fh,'Units',units);
        
        coeff  = (3 - 0.3*(fsize-6)/26);
        hshift = coeff * fsize/figpos(3);
        vshift = coeff * fsize/figpos(4);  
    catch
    end
    
    width  = ( Pos0(3) - (n2-1)*hshift ) / n2;      % subplot width
    height = ( Pos0(4) - (n1-1)*vshift ) / n1;      % subplot height
    
    ii = 0;
    for j = 1:1:n1              % rows
        for k = 1:1:n2          % columns
            ii = ii+1;
            pos(1) = Pos0(1) + (k-1)*(hshift+width);                    
            pos(2) = Pos0(2) + Pos0(4) - (j-1)*vshift - j*height;       
            pos(3) = width;     pos(4) = height;
            ax(ii) = subplot('Position',pos);
            if (ii == N),    break,      end,
        end
        if (ii == N),    break,      end,
    end
    
    % resize last subplot to fill page in horizontal direction
    if (n1*n2 > N)                
        pos = get(ax(N),'Position');              
        pos(3) = Pos0(1) + Pos0(3) - pos(1);
        set(ax(N),'Position',pos);
    end
end     % end of the function <f_subplot>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  f_plot_1d(pfunc,fselnames,Funcs)
%%% plot 1d graphs.
%%% Arguments: plot function, selected function names, functions(structure)
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
            
    titstr = '';    legstr = {};
    if (strcmp(MALER.PPR.absciss,MALER.WSP.Grid.name{1}))
        x = MALER.WSP.Grid.val{1};
    else
        ind = find((strcmp(MALER.PPR.absciss,MALER.WSP.Func.name)),1,'first');
        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',MALER.PPR.absciss]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   return,                    
        end
        x = MALER.WSP.Func.val{ind};
    end

    for j = 1 : 1 : length(fselnames)

        fnc = fselnames{j};     % the name of one function
        ind = find(strcmp(Funcs.name,fnc),1,'first');

        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',fnc]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   continue,                    
        end

        f = Funcs.val{ind};        % the value of the function
        MALER.HDL.Plot.data(j) = pfunc(x,f,'LineWidth',2,'Color',MALER.PPR.linecolor{length(fselnames),j});
        titstr = [ titstr , [fnc,', '] ];
        legstr = [ legstr , fnc ];
    end             
    axis tight, grid on, view(2),

    MALER.HDL.Plot.xlabel = xlabel(MALER.PPR.absciss);
    if ( ~isempty(titstr) )
        titstr(end-1:end) = []; MALER.HDL.Plot.tit = title(titstr,'Interpreter','none'); 
    end
    if ( numel(legstr) > 1)
        MALER.HDL.Plot.lgd = legend(legstr,'Location','best');    
    end

end     % end of the function <f_plot_1d>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  f_plot_slice(pfunc,fselnames)
%%% plot 1d slices.
%%% Arguments: plot function, selected function names
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    titstr = '';    legstr = {};
    ind = find((strcmp(MALER.PPR.absciss,union(MALER.WSP.Func.name,MALER.WSP.Grid.name))),1,'first');
    if (isempty(ind))
        h = errordlg(['ERROR! Cannot find function ',MALER.PPR.absciss]); 
        h = f_resize_dlgbox(h,fsize);   

        MALER.GUI.Listbox.absciss.Value = 1;
        MALER.PPR.absciss = MALER.GUI.Listbox.absciss.String{MALER.GUI.Listbox.absciss.Value};
        set(MALER.HDL.Listbox.absciss,'Value',MALER.GUI.Listbox.absciss.Value);
    end            

    if     (strcmp(MALER.PPR.absciss,MALER.WSP.Grid.name{1}))
        x = MALER.PPR.Slc.val{1};
    elseif (strcmp(MALER.PPR.absciss,MALER.WSP.Grid.name{2}))
        x = MALER.PPR.Slc.val{2};
    else               
        ind = find((strcmp(MALER.PPR.absciss,MALER.PPR.Slc.fname)),1,'first');
        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',MALER.PPR.absciss]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   return,                    
        end                                
        x = MALER.PPR.Slc.fval{ind};                              
    end           

    for j = 1 : 1 : length(fselnames)

        fnc = fselnames{j};     % the name of one function
        ind = find(strcmp(MALER.PPR.Slc.fname,fnc),1,'first');

        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',fnc]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   continue,                    
        end

        f = MALER.PPR.Slc.fval{ind};    % the value of the function
        MALER.HDL.Plot.data(j) = pfunc(x,f,'LineWidth',2,'Color',MALER.PPR.linecolor{length(fselnames),j}); 
        if ( max(x)-min(x) < 10*eps ) && ( max(f)-min(f) < 10*eps )
            set(MALER.HDL.Plot.data(j),'Marker','*','MarkerSize',10);
        end                                 
        titstr = [ titstr , [fnc,', '] ];
        legstr = [ legstr , fnc ];
    end              
    axis tight, grid on, view(2),

    MALER.HDL.Plot.xlabel = xlabel(MALER.PPR.absciss);
    if ( ~isempty(titstr) )
        titstr(end-1:end) = []; MALER.HDL.Plot.tit = title(titstr,'Interpreter','none'); 
    end
    if ( numel(legstr) > 1)
        MALER.HDL.Plot.lgd = legend(legstr,'Location','best');    
    end

end     % end of the function <f_plot_slice>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plot_2d_cont_mesh(x,y,fselnames,Funcs)
%%% plot 2d: contour, contourf, mesh, meshc
%%% Arguments: grid, selected function names, functs(struct)
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
            
    titstr = '';
    N = length(fselnames); 
    
    [n1,n2,ax] = f_subplot(N,MALER.GUI.AxeI.Position,fsize+1,MALER.HDL.Fig0);
    MALER.HDL.Plot.subplots = ax;
    
    %releasecond = ( N < 2 || ~MALER.uset.releaseCond ); % for colorbar
    releasecond = ( ~MALER.uset.releaseCond ); % for colorbar
    
    contlinenumber = max(10,round(40-20*(N/10)));   % for contour(f)
    for j = 1 : 1 : N

        fnc = fselnames{j};         % the name of one function
        ind = find(strcmp(Funcs.name,fnc),1,'first');
        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',fnc]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   continue,                    
        end 
        pfunc = eval(['@',MALER.PPR.func,';']);  % renew each time

        f = Funcs.val{ind};         % numerical values of the function
        aH = ax(j);                 % current subplot
        pos = get(aH,'Position');

        if (max(max(real(f)))-min(min(real(f)))<10*eps), pfunc = @mesh;end,

        if (~isempty(regexp(char(pfunc),'contour','once')))                                        
            [~,h] = pfunc(aH, x, y, f, contlinenumber, 'LineWidth',1); 
            view(2),    grid off,
            if strcmp(char(pfunc),'contour')
                set(h,'LineWidth',2);
            end
        else
            h = pfunc(aH, x, y, f);
            view(3);    grid on,    colormap jet;                        
        end
        axis tight;

        set(aH,'FontSize',MALER.GUI.AxeI.FontSize,'FontWeight',MALER.GUI.AxeI.FontWeight);
        MALER.HDL.Plot.data(j) = h; 
        MALER.HDL.Plot.tit(j)  = title(aH,fnc,'Interpreter','none');

        if (~isempty(regexp(char(pfunc),'contour','once')))
            if (pos(2) < MALER.GUI.AxeI.Position(2)+0.1)
                MALER.HDL.Plot.xlabel(j) = xlabel(aH,MALER.WSP.Grid.name{1});
            end
            if (pos(1) < MALER.GUI.AxeI.Position(1)+0.1)
                MALER.HDL.Plot.ylabel(j) = ylabel(aH,MALER.WSP.Grid.name{2},'Rotation',0,'Tag','YLabel');
            end                        
        else
            MALER.HDL.Plot.xlabel(j) = xlabel(aH,MALER.WSP.Grid.name{1});
            MALER.HDL.Plot.ylabel(j) = ylabel(aH,MALER.WSP.Grid.name{2},'Rotation',0,'Tag','YLabel');
        end 

        % make colorbar(s): Matlab release dependent procedure                                        
        if (releasecond)
            MALER.HDL.Plot.cbr(j) = colorbar('peer',aH,'FontSize',MALER.GUI.AxeI.FontSize,'FontWeight',MALER.GUI.AxeI.FontWeight);
        end                    
    end   % end of the for <selected functions> - loop 

    % make multiple colorbars in specific Matlab releases
    if (~releasecond)
        for j = 1:1:N
            axes(ax(j)), MALER.HDL.Plot.cbr(j) = colorbar(gca,'FontSize',MALER.GUI.AxeI.FontSize,'FontWeight',MALER.GUI.AxeI.FontWeight);  %#ok<LAXES>
        end
    end

end     % end of the function <f_plot_2d_cont_mesh>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plot_2d_surf(pfunc,x,y,fselnames,Funcs)
%%% plot 2d: surf, surfc
%%% Arguments: plot function, grid, selected function names, functs(struct)
    global MALER
    fsize = MALER.GUI.Settings.fontsize;

    titstr = '';
    N = length(fselnames);

    aH = MALER.HDL.AxeI;  
    for j = 1 : 1 : N
        fnc = fselnames{j};       % the name of one function
        ind = find(strcmp(Funcs.name,fnc),1,'first');
        if (isempty(ind))
            h = errordlg(['ERROR! Cannot find function ',fnc]); 
            h = f_resize_dlgbox(h,fsize);   uiwait(h);   continue,                    
        end 

        f = Funcs.val{ind};       % numerical values of the function
        h = pfunc(aH, x, y, f);   hold on,                    
        view(3);   grid on,   colormap(autumn(70)),  shading interp,

        titstr = [ titstr , [fnc,', '] ]; 
        MALER.HDL.Plot.data(j) = h;
     end    % end of for-loop

     axis tight;
     set(aH,'FontSize',MALER.GUI.AxeI.FontSize,'FontWeight',MALER.GUI.AxeI.FontWeight);
     if (~isempty(titstr))
        titstr(end-1:end) = [];  MALER.HDL.Plot.tit = title(aH,titstr,'Interpreter','none');
     end
     MALER.HDL.Plot.xlabel = xlabel(aH,MALER.WSP.Grid.name{1});
     MALER.HDL.Plot.ylabel = ylabel(aH,MALER.WSP.Grid.name{2},'Rotation',0); 
     MALER.HDL.Plot.cbr = colorbar();
     set(MALER.HDL.Plot.cbr,'FontSize',MALER.GUI.AxeI.FontSize,'FontWeight',MALER.GUI.AxeI.FontWeight);    

end     % end of the function <f_plot_2d_surf>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_close_figobj(hObj)
%%% close figure objects

    if ( ~isstruct(hObj) )  
        for k = 1:numel(hObj)
            if ( ishandle(hObj(k)) )
                try     delete(hObj(k)),
                catch,  continue,
                end
            end
        end      
              
    else  
        fields = fieldnames(hObj);
        for j = 1:numel(fields)
            h = hObj.(fields{j});
            for k = 1:numel(h)
                if (ishandle(h(k)))  
                    try     f_close_figobj(h(k));  
                    catch,  continue,
                    end
                end
            end
        end
    end
end     % end of the function <f_close_figobj>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%   SETUP MALER.GUI   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gui_setup()
%%% setup gui
    global MALER 

    MALER.GUI.Exist = false;
    
    %%% get workspace variables
    DefaultDir = MALER.GUI.Defaults.dir;

    FontSize  = MALER.GUI.Settings.fontsize;
    CustFSize = f_uifsz(FontSize);

    UserButts = MALER.GUI.Settings.userbutts;
    
    %%% main figure propeties
    MALER.GUI.Fig0.Name        = MALER.uset.progVersion;     
    MALER.GUI.Fig0.NumberTitle = 'Off';
    MALER.GUI.Fig0.Units       = 'Normalized';     MALER.GUI.Fig0.Resize      = 'On';
    MALER.GUI.Fig0.NextPlot    = 'Add';            MALER.GUI.Fig0.Toolbar     = 'Figure';
    MALER.GUI.Fig0.MenuBar     = 'Figure';         MALER.GUI.Fig0.DeleteFcn   = 'delete *.asv';
    MALER.GUI.Fig0.Position    = MALER.GUI.Settings.position; 
    MALER.GUI.Fig0.CloseRequestFcn = @f_callback_fig0_closereqfcn;
    MALER.GUI.Fig0.ResizeFcn   = @f_callback_fig0_resizefcn; 
    MALER.GUI.Fig0.Color = [0.94, 0.9400, 0.9400]; %[1, 130, 87]/255;

    %%% containing axes object properties
    MALER.GUI.Axe0.Units      = MALER.GUI.Fig0.Units;
    MALER.GUI.Axe0.Position   = [0.11 0.13 0.72 0.85];
    MALER.GUI.Axe0.XTick      = [];       MALER.GUI.Axe0.XTickLabel = [];
    MALER.GUI.Axe0.YTick      = [];       MALER.GUI.Axe0.YTickLabel = [];
    MALER.GUI.Axe0.Box        = 'On';
    MALER.GUI.Axe0.Color      = [212 214 219]/255;
    MALER.GUI.Axe0.Tag        = 'Axe0';

    %%% visible axes object properties
    MALER.GUI.AxeI.Units = MALER.GUI.Fig0.Units;
    MALER.GUI.AxeI.Position = [MALER.GUI.Axe0.Position(1)+0.1*MALER.GUI.Axe0.Position(3) , ...
                         MALER.GUI.Axe0.Position(2)+0.1*MALER.GUI.Axe0.Position(4) , ...
                         0.8*MALER.GUI.Axe0.Position(3) , 0.8*MALER.GUI.Axe0.Position(4)];
    MALER.GUI.AxeI.FontSize   = FontSize+1;
    MALER.GUI.AxeI.FontWeight = 'Bold';
    MALER.GUI.AxeI.NextPlot   = 'add';
    if (strcmp(MALER.PPR.func,'plot'))
        MALER.GUI.AxeI.XGrid = 'On'; MALER.GUI.AxeI.YGrid = 'On'; MALER.GUI.AxeI.ZGrid = 'On';
    else
        MALER.GUI.AxeI.XGrid = 'Off'; MALER.GUI.AxeI.YGrid = 'Off'; MALER.GUI.AxeI.ZGrid = 'Off';
    end
    MALER.GUI.AxeI.GridLineStyle = '--';
    MALER.GUI.AxeI.LineWidth = 0.5;
    MALER.GUI.AxeI.Visible = 'On';
    MALER.GUI.AxeI.Tag = 'AxeI';

    %%% --------   setup parameters of the active toolbar buttons   -------
    ii = 0;
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'NewPar';
    TbbT_TxtStr   {ii} = 'Add/Remove Parameters';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_filtcoef.gif'); 
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'help','signal','csh','button_filtcoef.gif');
    TbbT_imgColor {ii} = [255, 51, 51];   TbbT_mapset{ii} = [];      % red
    %-----    
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'NewFunc';
    TbbT_TxtStr   {ii} = 'Add/Change/Remove Function';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_function.gif');
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','shared','dastudio','resources','Function.gif');
    TbbT_imgColor {ii} = [51, 153, 255]; TbbT_mapset{ii} = []; % light-blue
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'PlotType';
    TbbT_TxtStr   {ii} = 'Plot Type';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_box.gif');  
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','matlab','icons','HDF_rasterimage.gif'); 
    TbbT_imgColor {ii} = [5, 204, 51]; TbbT_mapset{ii} = [];  % light-green
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'Save';
    TbbT_TxtStr   {ii} = 'Save';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_save.gif');
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','shared','dastudio','resources','save.gif');
    TbbT_imgColor {ii} = [0, 255, 255];    TbbT_mapset{ii} = [];    % cyan
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'FontSize';
    TbbT_TxtStr   {ii} = 'Reset Font Size';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_tool_text.gif');
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','matlab','icons','tool_text.gif');
    TbbT_imgColor {ii} = [0, 0, 0];    TbbT_mapset{ii} = [];    % black    
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'Digits';
    TbbT_TxtStr   {ii} = 'Reset Sliders Display Accuracy';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_target.gif');
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','shared','dastudio','resources','Target.gif');
    TbbT_imgColor {ii} = [255, 51, 255];   TbbT_mapset{ii} = [];  % magenta    
    %----- 
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'Defaults';
    TbbT_TxtStr   {ii} = 'Reset Default Settings';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_reload.gif'); 
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'help','matlab','reload.gif');
    TbbT_imgColor {ii} = [255, 204, 0];   TbbT_mapset{ii} = [];  % orange    
    %-----
    
    ii = ii+1;
    TbbT_FieldName{ii} = 'Help';
    TbbT_TxtStr   {ii} = 'Open User Guide';
    TbbT_ImgPath1 {ii} = fullfile(DefaultDir,'img_query.gif');
    TbbT_ImgPath2 {ii} = fullfile(matlabroot,'toolbox','matlab','icons','csh_icon.gif');
    TbbT_imgColor {ii} = [153, 153, 153];   TbbT_mapset{ii} = [];  % grey    
    %%%--------------------------------------------------------------------
    
    for j = 1:1:ii
        %TbbT_TxtStr{j} = ['<html><b><font size=5>',TbbT_TxtStr{j},'</font></b></html>'];
        TbbT_TxtStr{j} = ['<html><font size=',num2str(CustFSize,'%d'),'>',TbbT_TxtStr{j},'</font></html>'];
    end

    %%% setup GUI.TbarButt
    for ii = 1 : 1 : numel(TbbT_FieldName)        
        try     [imgData, map] = imread(TbbT_ImgPath1{ii});
        catch  
            try   
                [imgData, map] = imread(TbbT_ImgPath2{ii});
                [status,msg,msgID] = copyfile(TbbT_ImgPath2{ii}, TbbT_ImgPath1{ii}, 'f'); 
            catch,    imgData = TbbT_imgColor{ii};   map = TbbT_mapset{ii};
            end
        end
        
        clear ButtIcon
        if ( isempty(map) )           
            ButtIcon(:,:,1) = ones(16,16)/255*imgData(1); 
            ButtIcon(:,:,2) = ones(16,16)/255*imgData(2);
            ButtIcon(:,:,3) = ones(16,16)/255*imgData(3);
        else            
            ButtIcon = ind2rgb(imgData,map);
        end 
        
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).Tag = lower(TbbT_FieldName{ii});        
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).CData = ButtIcon;
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).Enable = 'on';
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).Visible = 'on';
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).TooltipString = TbbT_TxtStr{ii};        
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).ClickedCallback = ...
            eval(['@f_callback_tbarbutt_',lower(TbbT_FieldName{ii}),';']);
        MALER.GUI.TbarButt.(TbbT_FieldName{ii}).Interruptible = 'off';

    end
    if (MALER.DIM < 2),  MALER.GUI.TbarButt.PlotType.Enable  = 'off';  end,
    %%%--------------------------------------------------------------------   

    %%% setup user defined GUI.TbarButt (works correctly up to 9)
    for iK = 1 : 1 : UserButts
        Img1 = fullfile(DefaultDir,['img_',num2str(iK),'.gif']);
        Img2 = fullfile(matlabroot,'sys','namespace','docbook',...
                'v4','xsl','images','callouts',[num2str(iK),'.gif']);
        try     [imgData, map] = imread( Img1 );
        catch
            try 
                [imgData, map] = imread( Img2 );
                [status,msg,msgID] = copyfile(Img2, Img1, 'f');
            catch
                map = [];
                switch iK
                    case 1,     imgData = [255, 0, 0];          % red
                    case 2,     imgData = [0, 255, 0];          % green
                    case 3,     imgData = [0, 0, 255];          % blue
                    case 4,     imgData = [255, 204, 0];        % orange
                    case 5,     imgData = [255, 51, 255];       % magenta
                    case 6,     imgData = [0, 0, 0];            % black
                    case 7,     imgData = [0, 255, 255];        % cyan
                    otherwise,  imgData = [255, 255, 255];      % white
                end
            end
        end
        
        clear ButtIcon
        if ( isempty(map) )
            ButtIcon(:,:,1) = ones(16,16)/255*imgData(1);
            ButtIcon(:,:,2) = ones(16,16)/255*imgData(2);
            ButtIcon(:,:,3) = ones(16,16)/255*imgData(3);
        else
            ButtIcon = ind2rgb(imgData,map);
        end
        
        TxtStr = eval(['MALER.GUI.Settings.tooltip',num2str(iK,'%02d')]);
        TxtStr = ['<html><font size=',num2str(CustFSize,'%d'),'>',TxtStr,'</font></html>'];   
        
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.CData = ButtIcon;']);
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.ClickedCallback = ',...
              '@maler_userbutt_callbacks;']);
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.Tag = ',...
              '[''action'',num2str(iK,''%02d'')];']); 
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.Enable = ''on'';']);
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.Visible = ''on'';']);
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.TooltipString = ',...
              '''',TxtStr,'''',';']);
        eval(['MALER.GUI.TbarButt.UserAct_',num2str(iK,'%02d'),'.Interruptible = ''off'';']);
    end
    %%%--------------------------------------------------------------------

    %%% ---------   define grid slider and label parameters   -------------
    MALER.GUI.SliderGrid(1,1).Units = MALER.GUI.Fig0.Units;
    MALER.GUI.SliderGrid(1,1).BackgroundColor = [212 214 219]/255; % slider color
    MALER.GUI.SliderGrid(1,1).Callback = @f_callback_gridslider;
    MALER.GUI.SliderGrid(1,1).Tag = '11';
    MALER.GUI.SliderGrid(1,1).String = MALER.GUI.SliderGrid(1,1).Tag;
    MALER.GUI.SliderGrid(1,1).Interruptible = 'off';
        
    %%%--parameters of horizontal sliders and their labels-----------------
    MALER.GUI.sldgrid.hor.shift  = 0.04;         % horizontal slider horizontal shift
    MALER.GUI.sldgrid.hor.width  = (MALER.GUI.Axe0.Position(3)-2*MALER.GUI.sldgrid.hor.shift)/3; % horizontal slider width
    MALER.GUI.sldgrid.hor.height = MALER.GUI.Axe0.Position(1)/6;    % horizontal slider height
    if ( MALER.GUI.SCR.Prop > 0 )    
        MALER.GUI.sldgrid.hor.height = MALER.GUI.sldgrid.hor.height * sqrt(MALER.GUI.SCR.Width/MALER.GUI.SCR.Height);
    end
    MALER.GUI.sldgrid.hor.vbase = 0.05;            % horizontal slider vertical base

    MALER.GUI.SliderGrid(1,1).Position = [MALER.GUI.Axe0.Position(1), MALER.GUI.sldgrid.hor.vbase, MALER.GUI.sldgrid.hor.width, MALER.GUI.sldgrid.hor.height];
    MALER.GUI.SliderGrid(1,2) = MALER.GUI.SliderGrid(1,1);
    MALER.GUI.SliderGrid(1,2).Position(1) = MALER.GUI.SliderGrid(1,1).Position(1) + MALER.GUI.sldgrid.hor.width +MALER.GUI.sldgrid.hor.shift;
    MALER.GUI.SliderGrid(1,2).Tag = '12';
    MALER.GUI.SliderGrid(1,2).String = MALER.GUI.SliderGrid(1,2).Tag;
    
    MALER.GUI.SliderGrid(1,3) = MALER.GUI.SliderGrid(1,2);
    MALER.GUI.SliderGrid(1,3).Position(1) = MALER.GUI.SliderGrid(1,2).Position(1) + MALER.GUI.sldgrid.hor.width + MALER.GUI.sldgrid.hor.shift;
    MALER.GUI.SliderGrid(1,3).Tag = '13';
    MALER.GUI.SliderGrid(1,3).String = MALER.GUI.SliderGrid(1,3).Tag;
    
    MALER.GUI.LabelGrid(1,1).Units = MALER.GUI.Fig0.Units;
    MALER.GUI.LabelGrid(1,1).BackgroundColor = MALER.GUI.Fig0.Color;   % slider label color
    MALER.GUI.LabelGrid(1,1).FontSize = FontSize;
    MALER.GUI.LabelGrid(1,1).FontWeight = 'Bold';

    MALER.GUI.txtgrid.hor.width = MALER.GUI.sldgrid.hor.width;         % text box width
    MALER.GUI.txtgrid.hor.height = 2*FontSize/MALER.GUI.SCR.Height;    % text box height
    MALER.GUI.txtgrid.hor.vbase = MALER.GUI.sldgrid.hor.vbase + MALER.GUI.sldgrid.hor.height + 0.2*MALER.GUI.txtgrid.hor.height;      % text box vertical base

    MALER.GUI.LabelGrid(1,1).Position = [MALER.GUI.Axe0.Position(1) + (MALER.GUI.sldgrid.hor.width - MALER.GUI.txtgrid.hor.width)/2, MALER.GUI.txtgrid.hor.vbase, MALER.GUI.txtgrid.hor.width, MALER.GUI.txtgrid.hor.height];
    MALER.GUI.LabelGrid(1,2) = MALER.GUI.LabelGrid(1,1);
    MALER.GUI.LabelGrid(1,2).Position(1) = MALER.GUI.LabelGrid(1,1).Position(1) + MALER.GUI.sldgrid.hor.width + MALER.GUI.sldgrid.hor.shift;
    MALER.GUI.LabelGrid(1,3) = MALER.GUI.LabelGrid(1,2);
    MALER.GUI.LabelGrid(1,3).Position(1) = MALER.GUI.LabelGrid(1,2).Position(1) + MALER.GUI.sldgrid.hor.width + MALER.GUI.sldgrid.hor.shift;
    
    %%%--parameters of vertical sliders and their labels-------------------
    MALER.GUI.txtgrid.ver.width  = MALER.GUI.Axe0.Position(1)/1.2;   % text box width
    MALER.GUI.txtgrid.ver.height = 1.6*MALER.GUI.txtgrid.hor.height;       %text box height
    MALER.GUI.sldgrid.ver.shift = 0.03*min(FontSize/14,1); % vertical shift between vertical sliders
    MALER.GUI.sldgrid.ver.width = MALER.GUI.sldgrid.hor.height;   % vertical slider width
    if ( MALER.GUI.SCR.Prop > 0 )
        MALER.GUI.sldgrid.ver.width = MALER.GUI.sldgrid.ver.width/MALER.GUI.SCR.Width*MALER.GUI.SCR.Height;       % vertical slider width
    end
    MALER.GUI.sldgrid.ver.height = (MALER.GUI.Axe0.Position(4)-3*MALER.GUI.txtgrid.ver.height-2*MALER.GUI.sldgrid.ver.shift)/3;  %slider height
    MALER.GUI.sldgrid.ver.hbase  = (MALER.GUI.Axe0.Position(1) - MALER.GUI.sldgrid.ver.width)/2; % horiz. base
    MALER.GUI.txtgrid.ver.hbase  = MALER.GUI.Axe0.Position(1)/12;   % box horizontal base

    MALER.GUI.SliderGrid(2,1) = MALER.GUI.SliderGrid(1,1);
    MALER.GUI.SliderGrid(2,1).Position = [MALER.GUI.sldgrid.ver.hbase, MALER.GUI.Axe0.Position(2), MALER.GUI.sldgrid.ver.width, MALER.GUI.sldgrid.ver.height];
    MALER.GUI.SliderGrid(2,1).Tag = '21';
    MALER.GUI.SliderGrid(2,1).String = MALER.GUI.SliderGrid(2,1).Tag;
    
    MALER.GUI.SliderGrid(2,2) = MALER.GUI.SliderGrid(2,1);
    MALER.GUI.SliderGrid(2,2).Position(2) = MALER.GUI.SliderGrid(2,1).Position(2) + MALER.GUI.sldgrid.ver.height + MALER.GUI.txtgrid.ver.height + MALER.GUI.sldgrid.ver.shift;
    MALER.GUI.SliderGrid(2,2).Tag = '22';
    MALER.GUI.SliderGrid(2,2).String = MALER.GUI.SliderGrid(2,2).Tag;    
    
    MALER.GUI.SliderGrid(2,3) = MALER.GUI.SliderGrid(2,2);
    MALER.GUI.SliderGrid(2,3).Position(2) = MALER.GUI.SliderGrid(2,2).Position(2) + MALER.GUI.sldgrid.ver.height + MALER.GUI.txtgrid.ver.height + MALER.GUI.sldgrid.ver.shift;
    MALER.GUI.SliderGrid(2,3).Tag = '23';
    MALER.GUI.SliderGrid(2,3).String = MALER.GUI.SliderGrid(2,3).Tag;
    
    MALER.GUI.LabelGrid(2,1) = MALER.GUI.LabelGrid(1,1);
    MALER.GUI.LabelGrid(2,1).Position = [MALER.GUI.txtgrid.ver.hbase, (MALER.GUI.Axe0.Position(2)+MALER.GUI.sldgrid.ver.height), MALER.GUI.txtgrid.ver.width, MALER.GUI.txtgrid.ver.height];
    MALER.GUI.LabelGrid(2,2) = MALER.GUI.LabelGrid(2,1);
    MALER.GUI.LabelGrid(2,2).Position(2) = MALER.GUI.LabelGrid(2,1).Position(2) + MALER.GUI.sldgrid.ver.height + MALER.GUI.txtgrid.ver.height + MALER.GUI.sldgrid.ver.shift;
    MALER.GUI.LabelGrid(2,3) = MALER.GUI.LabelGrid(2,2);
    MALER.GUI.LabelGrid(2,3).Position(2) = MALER.GUI.LabelGrid(2,2).Position(2) + MALER.GUI.sldgrid.ver.height + MALER.GUI.txtgrid.ver.height + MALER.GUI.sldgrid.ver.shift;

    f_gridsliders_setup();  % setup slider values and text labels
    %%%--------------------------------------------------------------------
    
    %%%-parameters of the button 'choose absciss'--------------------------
    MALER.GUI.FigButt.ChAbs.Units = MALER.GUI.Fig0.Units;
    MALER.GUI.FigButt.ChAbs.Position = [MALER.GUI.LabelGrid(2,1).Position(1), MALER.GUI.sldgrid.hor.vbase, ...
        0.9*(MALER.GUI.Axe0.Position(1)-MALER.GUI.LabelGrid(2,1).Position(1)),0.65*(MALER.GUI.Axe0.Position(2)-MALER.GUI.sldgrid.hor.vbase)];
    MALER.GUI.FigButt.ChAbs.FontSize = MALER.GUI.LabelGrid(1,1).FontSize;
    MALER.GUI.FigButt.ChAbs.FontWeight = MALER.GUI.LabelGrid(1,1).FontWeight;
    MALER.GUI.FigButt.ChAbs.String = {'Abscissa'};
    MALER.GUI.FigButt.ChAbs.Callback = @f_callback_pushbutt_absciss;
    MALER.GUI.FigButt.ChAbs.tag = 'chabs';
    MALER.GUI.FigButt.ChAbs.Backgroundcolor = 'blue';   % button color
    MALER.GUI.FigButt.ChAbs.Foregroundcolor = 'white';  % button label color
    MALER.GUI.FigButt.ChAbs.Interruptible = 'off';
    if (MALER.DIM>1)
            MALER.GUI.FigButt.ChAbs.Visible = 'on';     MALER.GUI.FigButt.ChAbs.Enable = 'off';   
    else    MALER.GUI.FigButt.ChAbs.Visible = 'on';     MALER.GUI.FigButt.ChAbs.Enable = 'on';    
    end
    %%%--------------------------------------------------------------------

    %%%-parameters of the button 'setup grid parameters'-------------------
    MALER.GUI.rightpanel.width = 0.9*(1 - MALER.GUI.Axe0.Position(1) - MALER.GUI.Axe0.Position(3));
    MALER.GUI.rightpanel.hbase = MALER.GUI.Axe0.Position(1) + MALER.GUI.Axe0.Position(3) + ...
        0.5*(1 - MALER.GUI.Axe0.Position(1) - MALER.GUI.Axe0.Position(3)) - ...
        0.5*MALER.GUI.rightpanel.width;

    MALER.GUI.FigButt.SetGrid = MALER.GUI.FigButt.ChAbs;
    MALER.GUI.FigButt.SetGrid.Position(1) = MALER.GUI.rightpanel.hbase;
    MALER.GUI.FigButt.SetGrid.Position(3) = MALER.GUI.rightpanel.width;
    MALER.GUI.FigButt.SetGrid.String = {'Setup Grid Variables'};
    MALER.GUI.FigButt.SetGrid.Enable = 'on';
    MALER.GUI.FigButt.SetGrid.Visible = 'on';
    MALER.GUI.FigButt.SetGrid.Callback = @f_callback_pushbutt_setgrid;
    MALER.GUI.FigButt.SetGrid.Interruptible = 'off';
    %%%--------------------------------------------------------------------

    %%%-parameters of the button 'setup parameter values'------------------
    MALER.GUI.FigButt.SetParam = MALER.GUI.FigButt.SetGrid;
    MALER.GUI.FigButt.SetParam.Position(2) = MALER.GUI.FigButt.SetGrid.Position(2) + 1.1*MALER.GUI.FigButt.SetGrid.Position(4);
    MALER.GUI.FigButt.SetParam.String = {'Setup Parameters'};
    MALER.GUI.FigButt.SetParam.Callback = @f_callback_pushbutt_setparam;
    if (isempty(MALER.WSP.Param.name))   
            MALER.GUI.FigButt.SetParam.Enable = 'off';    MALER.GUI.FigButt.SetParam.Visible = 'off';
    else    MALER.GUI.FigButt.SetParam.Enable = 'on';     MALER.GUI.FigButt.SetParam.Visible = 'on';
    end
    %%%--------------------------------------------------------------------

    %%% parameters of the text label of the functions listbox
    MALER.GUI.TextFuncs = MALER.GUI.LabelGrid(2,3);
    MALER.GUI.TextFuncs.Position(1) = MALER.GUI.rightpanel.hbase;
    MALER.GUI.TextFuncs.Position(3) = MALER.GUI.rightpanel.width;
    MALER.GUI.TextFuncs.String = 'plot funcs';
    MALER.GUI.TextFuncs.Enable = 'on';
    MALER.GUI.TextFuncs.Visible = 'on';
    MALER.GUI.TextFuncs.Callback = '';
    MALER.GUI.TextFuncs.ButtonDownFcn = @f_callback_textfunc_bdfcn;
    %MALER.GUI.TextFuncs.TooltipString = ['<html><b><font size=5>','Rightclick to open hide/show menu','</font></b></html>'];
    MALER.GUI.TextFuncs.Selected = 'off';
    
    %%%-parameters of the listbox 'choose functions to plot'---------------
    MALER.GUI.Listbox.funcs.Units = MALER.GUI.Fig0.Units;
    MALER.GUI.Listbox.funcs.FontSize = FontSize;
    MALER.GUI.Listbox.funcs.FontWeight = 'Bold';
    MALER.GUI.Listbox.funcs.Backgroundcolor = MALER.GUI.FigButt.SetGrid.Backgroundcolor;
    MALER.GUI.Listbox.funcs.Foregroundcolor = MALER.GUI.FigButt.SetGrid.Foregroundcolor;
    MALER.GUI.Listbox.funcs.Position = MALER.GUI.TextFuncs.Position;
    if (isempty(MALER.WSP.Param.name)),     maxheight = 0.4;
    else                                    maxheight = 0.2;
    end
    MALER.GUI.Listbox.funcs.Position(4) = max(0.1,min(maxheight,2*length(MALER.PPR.List.show)*FontSize/MALER.GUI.SCR.Size(4)));  
    MALER.GUI.Listbox.funcs.Position(2) = MALER.GUI.Listbox.funcs.Position(2)-MALER.GUI.Listbox.funcs.Position(4);
    MALER.GUI.Listbox.funcs.String = MALER.PPR.List.show;
    MALER.GUI.Listbox.funcs.Callback = @f_callback_listbox_funcs;     
    MALER.GUI.Listbox.funcs.Min = 0;
    MALER.GUI.Listbox.funcs.Max = max(1,numel(MALER.GUI.Listbox.funcs.String)+1);
    MALER.GUI.Listbox.funcs.Enable  = 'on';     
    MALER.GUI.Listbox.funcs.Visible = 'on';
    MALER.GUI.Listbox.funcs.Interruptible = 'off';
    
    %%%-parameters of the listbox 'hide/show functions'--------------------
    MALER.GUI.Listbox.hideshow = MALER.GUI.Listbox.funcs;
    MALER.GUI.Listbox.hideshow.Enable = 'off';
    MALER.GUI.Listbox.hideshow.Visible = 'off';
    MALER.GUI.Listbox.hideshow.Backgroundcolor = [1, 1, 0];
    MALER.GUI.Listbox.hideshow.Foregroundcolor = [0, 0, 0];
    MALER.GUI.Listbox.hideshow.Callback = @f_callback_listbox_hideshow;
    MALER.GUI.Listbox.hideshow.Interruptible = 'off';

    %%% parameters of the text label of the parameters listbox
    MALER.GUI.TextParam = MALER.GUI.TextFuncs;
    MALER.GUI.TextParam.Position(2) = MALER.GUI.Listbox.funcs.Position(2) - 1.4*MALER.GUI.TextParam.Position(4);
    MALER.GUI.TextParam.String = 'tune params';
    MALER.GUI.TextParam = rmfield(MALER.GUI.TextParam,'ButtonDownFcn');
    MALER.GUI.TextParam.TooltipString = '';
    MALER.GUI.TextParam.Selected = 'off';

    %%%-parameters of the listbox 'choose parameters to tune'--------------
    MALER.GUI.Listbox.param = MALER.GUI.Listbox.funcs;
    MALER.GUI.Listbox.param.FontSize = FontSize;
    if (isempty(MALER.WSP.Func.name)),     maxheight = 0.3;
    else                                   maxheight = 0.15;
    end
    MALER.GUI.Listbox.param.Position(4) = max(0.1,min(maxheight,2*length(MALER.WSP.Param.name)*FontSize/MALER.GUI.SCR.Size(4)));
    MALER.GUI.Listbox.param.Position(2) = MALER.GUI.TextParam.Position(2) - MALER.GUI.Listbox.param.Position(4);
    MALER.GUI.Listbox.param.String = MALER.WSP.Param.name;
    MALER.GUI.Listbox.param.Callback = @f_callback_listbox_param;     
    MALER.GUI.Listbox.param.Min = 0;
    MALER.GUI.Listbox.param.Max = 3;

    switch (numel(MALER.WSP.Param.name))
        case 0
            MALER.GUI.Listbox.param.Enable = 'off';    MALER.GUI.Listbox.param.Visible = 'off';
            MALER.GUI.TextParam.Enable = 'off';        MALER.GUI.TextParam.Visible = 'off';
            MALER.GUI.Listbox.param.Value = [];
        case 1,     MALER.GUI.Listbox.param.Value = 1;
        case 2,     MALER.GUI.Listbox.param.Value = [1,2];
        otherwise,  MALER.GUI.Listbox.param.Value = [1,2,3];
    end

    switch (numel(MALER.WSP.Func.name))
        case 0
            MALER.GUI.Listbox.funcs.Enable = 'off';    MALER.GUI.Listbox.funcs.Visible = 'off';
            MALER.GUI.TextFuncs.Enable = 'off';        MALER.GUI.TextFuncs.Visible = 'off';
            MALER.GUI.Listbox.funcs.Value = [];
        otherwise,  MALER.GUI.Listbox.funcs.Value = 1;
    end
    %%%--------------------------------------------------------------------

    %-setup <tune parameters> sliders
    MALER.GUI.sldparampanel.vshift = (MALER.GUI.FigButt.SetParam.Position(2) - MALER.GUI.FigButt.SetGrid.Position(2) - MALER.GUI.FigButt.SetGrid.Position(4));
    MALER.GUI.sldparampanel.hbase = MALER.GUI.FigButt.SetGrid.Position(1);
    MALER.GUI.sldparampanel.width = MALER.GUI.FigButt.SetGrid.Position(3);
    MALER.GUI.sldparampanel.vbase = 2*MALER.GUI.FigButt.SetParam.Position(2) - MALER.GUI.FigButt.SetGrid.Position(2);
    MALER.GUI.sldparampanel.height = MALER.GUI.Listbox.param.Position(2) - MALER.GUI.sldparampanel.vbase - MALER.GUI.sldparampanel.vshift;
    
    f_paramsliders_setup();
    %%%--------------------------------------------------------------------
    
    %%%-parameters of the listbox 'choose absciss'-------------------------
    MALER.GUI.Listbox.absciss = MALER.GUI.Listbox.funcs;
    MALER.GUI.Listbox.absciss.String = union(MALER.WSP.Grid.name{1},MALER.WSP.Func.name,'stable');
    MALER.GUI.Listbox.absciss.Position(1) = MALER.GUI.FigButt.ChAbs.Position(1);
    MALER.GUI.Listbox.absciss.Position(3) = MALER.GUI.FigButt.ChAbs.Position(3);
    MALER.GUI.Listbox.absciss.Position(4) = min(0.85,2*(1+length(MALER.GUI.Listbox.absciss.String))*FontSize/MALER.GUI.SCR.Size(4));
    MALER.GUI.Listbox.absciss.Position(2) = MALER.GUI.FigButt.ChAbs.Position(2) + 1.05*MALER.GUI.FigButt.ChAbs.Position(4);        
    MALER.GUI.Listbox.absciss.Value = 1;                
    MALER.GUI.Listbox.absciss.Min = 1;
    MALER.GUI.Listbox.absciss.Max = 1;    
    MALER.GUI.Listbox.absciss.Callback = @f_callback_listbox_absciss; 
    MALER.GUI.Listbox.absciss.Enable  = 'off';
    MALER.GUI.Listbox.absciss.Visible = 'off';
    %%%--------------------------------------------------------------------
    
    %%%-parameters of the button 'Apply'-----------------------------------
    MALER.GUI.FigButt.Apply.Units = MALER.GUI.Fig0.Units;
    MALER.GUI.FigButt.Apply.Position = MALER.GUI.TextFuncs.Position;
    MALER.GUI.FigButt.Apply.Position(1) = MALER.GUI.FigButt.Apply.Position(1)+0.2*MALER.GUI.FigButt.Apply.Position(3);
    MALER.GUI.FigButt.Apply.Position(3) = 0.6*MALER.GUI.FigButt.Apply.Position(3);
    MALER.GUI.FigButt.Apply.String = 'Apply'; 
    MALER.GUI.FigButt.Apply.FontSize = FontSize;
    MALER.GUI.FigButt.Apply.FontWeight = 'Bold';
    MALER.GUI.FigButt.Apply.Backgroundcolor = MALER.GUI.Axe0.Color;
    MALER.GUI.FigButt.Apply.Foregroundcolor = [0, 0, 0];
    MALER.GUI.FigButt.Apply.Callback = @f_callback_figbutt_apply;
    MALER.GUI.FigButt.Apply.Enable = 'off';
    MALER.GUI.FigButt.Apply.Visible = 'off';
    MALER.GUI.FigButt.Apply.Interruptible = 'off';
    %%%--------------------------------------------------------------------
           
end     % end of the function <f_gui_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gridsliders_setup()
%%% setup grid sliders and their labels
    global MALER
          
    x = MALER.WSP.Grid.val{1,1};       
    vname = MALER.WSP.Grid.name{1};

    [value,vmin,vmax] = f_slider_update(min(x),[],1);
    MALER.GUI.SliderGrid(1,1).Value = value;   MALER.GUI.SliderGrid(1,1).Min = vmin;   MALER.GUI.SliderGrid(1,1).Max = vmax;

    [value,vmin,vmax] = f_slider_update( (max(x)-min(x))/max(1,(length(x)-1)) , [] , 2 );
    MALER.GUI.SliderGrid(1,2).Value = value;   MALER.GUI.SliderGrid(1,2).Min = vmin;   MALER.GUI.SliderGrid(1,2).Max = vmax;

    [value,vmin,vmax] = f_slider_update(max(x),[],3);
    MALER.GUI.SliderGrid(1,3).Value = value;   MALER.GUI.SliderGrid(1,3).Min = vmin;   MALER.GUI.SliderGrid(1,3).Max = vmax;
    
    numspec = MALER.uset.numSpec;
    
    MALER.GUI.LabelGrid(1,1).String = ['min ',vname,' = ',num2str( MALER.GUI.SliderGrid(1,1).Value, numspec )];
    MALER.GUI.LabelGrid(1,2).String = ['d',   vname,' = ',num2str( MALER.GUI.SliderGrid(1,2).Value, numspec )];
    MALER.GUI.LabelGrid(1,3).String = ['max ',vname,' = ',num2str( MALER.GUI.SliderGrid(1,3).Value, numspec )];

    if (~strcmp(MALER.PPR.func,'plot')),      
        y = MALER.WSP.Grid.val{2,1};       
        vname = MALER.WSP.Grid.name{2};

        [value,vmin,vmax] = f_slider_update(min(y),[],1);
        MALER.GUI.SliderGrid(2,1).Value = value;   MALER.GUI.SliderGrid(2,1).Min = vmin;   MALER.GUI.SliderGrid(2,1).Max = vmax;

        [value,vmin,vmax] = f_slider_update( (max(y)-min(y))/max(1,(length(y)-1)) , [] , 2 );
        MALER.GUI.SliderGrid(2,2).Value = value;   MALER.GUI.SliderGrid(2,2).Min = vmin;   MALER.GUI.SliderGrid(2,2).Max = vmax;

        [value,vmin,vmax] = f_slider_update(max(y),[],3);
        MALER.GUI.SliderGrid(2,3).Value = value;   MALER.GUI.SliderGrid(2,3).Min = vmin;   MALER.GUI.SliderGrid(2,3).Max = vmax;

        MALER.GUI.LabelGrid(2,1).String = ['min ',vname,' = ',num2str( MALER.GUI.SliderGrid(2,1).Value, numspec )];
        MALER.GUI.LabelGrid(2,2).String = ['d',   vname,' = ',num2str( MALER.GUI.SliderGrid(2,2).Value, numspec )];
        MALER.GUI.LabelGrid(2,3).String = ['max ',vname,' = ',num2str( MALER.GUI.SliderGrid(2,3).Value, numspec )];
    end
end     % end of the function <f_gridsliders_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,vmin,vmax] = f_slider_update(value,range,varargin)
%%% update slider parameters
%%% <value> is the value to be set, range = [vmin,vmax] or []
    if ( numel(varargin) > 0 )     
        i2 = varargin{1};       % second index of grid slider
    else
        i2 = -1;
    end
    switch numel(range)
        case 0
            if ( i2 == 2 ),     
                coeff = 4;  
                vmin = value/coeff;         vmax = value*coeff;
            else
                delta = max(1,abs(value));
                vmin = value - delta;       vmax = value + delta;
            end
        case 1
            vmin = min(range,value-range);  vmax = max(value,value+range);
            if ( abs(vmax-vmin)<10*eps ),   vmax = vmin + 1;    end,
        otherwise
            vmin = min(range(1:2));         vmax = max(range(1:2));
            if ( abs(vmax-vmin)<10*eps ),   vmax = vmin + 1;    end,
    end
end    % end of the 'f_slider_update' function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = f_create_gui_object(obj , varargin)
%%% create gui element <obj> (type string)
global MALER                                                   
    answer = regexp(obj,'\.','split');      str = '';
    for j = 1:1:length(answer)
        str = [str,'.(''',answer{j},''')'];
    end
    source = eval(['MALER.GUI',str,';']);
    target = eval(['MALER.HDL',str,';']);
    pnames = fieldnames(source);        % get settings' names
    switch numel(varargin)
        case 0
            for j = 1:1:length(pnames)
                try  
                    set( target , pnames{j} , source.(pnames{j}));                      
                catch
                end
            end
        case 1
            handle = varargin{1}; 
            for j = 1:1:length(pnames)
                try      
                    set( handle , pnames{j} , MALER.GUI.(obj).(pnames{j}));
                catch
                end
            end
            varargout{1} = handle;    
        case 2
            handle = varargin{1};   i1 = varargin{2}; 
            for j = 1:1:length(pnames)
                try  
                    set( handle , pnames{j} , source(i1).(pnames{j}));                        
                catch
                end
            end
            varargout{1} = handle;    
        case 3
            handle = varargin{1};   i1 = varargin{2};   i2 = varargin{3}; 
            for j = 1:1:length(pnames)
                try      set( handle , pnames{j} , MALER.GUI.(obj)(i1,i2).(pnames{j}));
                catch
                end
            end
            varargout{1} = handle;
        otherwise
    end
end     % end of the function <f_create_gui_object>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_grid_sliders()
%%% create grid-variable sliders
    global MALER 
    
    for i1 = 1:1:min(2,MALER.DIM)
        for i2 = 1:1:3
            ind = [num2str(i1),num2str(i2)];
            MALER.HDL.SliderGrid(i1,i2) = uicontrol('Parent',MALER.HDL.Fig0,'Style','slider');
            obj = 'SliderGrid';
            MALER.HDL.SliderGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.SliderGrid(i1,i2),i1,i2);
            
            MALER.HDL.LabelGrid(i1,i2) = uicontrol('Parent',MALER.HDL.Fig0,'Style','text');
            obj = 'LabelGrid';
            MALER.HDL.LabelGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.LabelGrid(i1,i2),i1,i2);
        end
    end
end     % end of the function <f_create_grid_sliders>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_paramsliders_setup(varargin)
%%% setup parameters of the <tune> parameter sliders
    if (~isempty(varargin)),    mindrange = true;   
    else                        mindrange = false;
    end
    global MALER
    
    %%% remove old sliders and labels if exist
    if (isfield(MALER.HDL,'SliderParam'))
        delete(MALER.HDL.SliderParam),      delete(MALER.HDL.LabelParam);
        MALER.HDL = rmfield(MALER.HDL,{'SliderParam','LabelParam'});
        MALER.GUI = rmfield(MALER.GUI,{'SliderParam','LabelParam'});
    end
    
    numspec = MALER.uset.numSpec;
    select  = MALER.GUI.Listbox.param.Value;
    Lth = numel(select);
    
    if (select < 1),   return,     end,
    
    if (Lth > 0)
        MALER.GUI.LabelParam(1) = MALER.GUI.TextParam;
        MALER.GUI.SliderParam(1) = MALER.GUI.SliderGrid(1,1);
        ii = 1; vname = MALER.GUI.Listbox.param.String{select(ii)};
        f_assign_sliderparam_values(ii,vname,select,mindrange,numspec);

        MALER.GUI.LabelParam(1).Position(1) = MALER.GUI.sldparampanel.hbase;
        MALER.GUI.LabelParam(1).Position(3) = MALER.GUI.sldparampanel.width;
        MALER.GUI.LabelParam(1).Position(2) = MALER.GUI.sldparampanel.vbase + MALER.GUI.sldparampanel.height - 1.2*MALER.GUI.LabelParam(1).Position(4); 
       
        MALER.GUI.SliderParam(1).Position(3) = 0.1*MALER.GUI.sldparampanel.width;
        MALER.GUI.SliderParam(1).Position(1) = MALER.GUI.sldparampanel.hbase+0.5*(MALER.GUI.sldparampanel.width - MALER.GUI.SliderParam(1).Position(3));
        MALER.GUI.SliderParam(1).Position(2) = MALER.GUI.sldparampanel.vbase;
        MALER.GUI.SliderParam(1).Position(4) = MALER.GUI.sldparampanel.height - MALER.GUI.LabelParam(1).Position(4) - 2*MALER.GUI.sldparampanel.vshift;
        MALER.GUI.SliderParam(1).Callback = @f_callback_paramslider;
        MALER.GUI.SliderParam(1).String = vname;
        MALER.GUI.SliderParam(1).Tag    = '31';
        MALER.GUI.SliderParam(1).Interruptible = 'off';
        MALER.GUI.SliderParam(1).Enable = MALER.GUI.LabelParam(1).Enable;
    end
    
    if (Lth > 1)
        MALER.GUI.LabelParam(1).Position(3) = (2/3)*MALER.GUI.sldparampanel.width-MALER.GUI.SliderParam(1).Position(3);
        MALER.GUI.LabelParam(1).Position(1) = MALER.GUI.LabelParam(1).Position(1)+(1/2)*MALER.GUI.SliderParam(1).Position(3);
        MALER.GUI.LabelParam(1).Position(2) = MALER.GUI.LabelParam(1).Position(2) - MALER.GUI.LabelParam(1).Position(4);
        
        MALER.GUI.SliderParam(1).Position(1) = MALER.GUI.SliderParam(1).Position(1) - (1/2-1/3)*MALER.GUI.sldparampanel.width;
        MALER.GUI.SliderParam(1).Position(4) = MALER.GUI.SliderParam(1).Position(4) - MALER.GUI.LabelParam(1).Position(4);
        %-----------------------------------
        
        MALER.GUI.LabelParam(2) = MALER.GUI.LabelParam(1);
        MALER.GUI.SliderParam(2) = MALER.GUI.SliderParam(1);
        ii = 2; vname = MALER.GUI.Listbox.param.String{select(ii)};
        f_assign_sliderparam_values(ii,vname,select,mindrange,numspec);
        
        MALER.GUI.LabelParam(2).Position(1) = MALER.GUI.LabelParam(1).Position(1) + (1/3)*MALER.GUI.sldparampanel.width;
        MALER.GUI.LabelParam(2).Position(2) = MALER.GUI.LabelParam(1).Position(2) + MALER.GUI.LabelParam(1).Position(4);
                
        MALER.GUI.SliderParam(2).Position(1) = MALER.GUI.SliderParam(2).Position(1) + (1/3)*MALER.GUI.sldparampanel.width;
        MALER.GUI.SliderParam(2).Position(4) = MALER.GUI.SliderParam(1).Position(4) + MALER.GUI.LabelParam(2).Position(4);
        MALER.GUI.SliderParam(2).String = vname;
        MALER.GUI.SliderParam(2).Tag    = '32';
        MALER.GUI.SliderParam(2).Interruptible = 'off';
        MALER.GUI.SliderParam(2).Enable = MALER.GUI.LabelParam(2).Enable;
    end
    
    if (Lth>2)
        
        MALER.GUI.LabelParam(1).Position(3) = (1/2)*MALER.GUI.sldparampanel.width - MALER.GUI.SliderParam(1).Position(3);
        MALER.GUI.LabelParam(1).Position(1) = MALER.GUI.sldparampanel.hbase+(1/2)*MALER.GUI.SliderParam(1).Position(3);
        
        MALER.GUI.SliderParam(1).Position(1) = MALER.GUI.sldparampanel.hbase + (1/4)*MALER.GUI.sldparampanel.width - (1/2)*MALER.GUI.SliderParam(1).Position(3);
        %-----------------------------------
        
        MALER.GUI.SliderParam(2).Position(1) = MALER.GUI.SliderParam(1).Position(1) + (1/4)*MALER.GUI.sldparampanel.width;
        
        MALER.GUI.LabelParam(2).Position(3) = MALER.GUI.LabelParam(1).Position(3);
        MALER.GUI.LabelParam(2).Position(1) = MALER.GUI.SliderParam(1).Position(1) + MALER.GUI.SliderParam(1).Position(3);
        %-----------------------------------
        
        MALER.GUI.LabelParam(3) = MALER.GUI.LabelParam(2);
        MALER.GUI.SliderParam(3) = MALER.GUI.SliderParam(2);
        ii = 3; vname = MALER.GUI.Listbox.param.String{select(ii)};
        f_assign_sliderparam_values(ii,vname,select,mindrange,numspec);

        MALER.GUI.LabelParam(3).Position(1) = MALER.GUI.SliderParam(2).Position(1) + MALER.GUI.SliderParam(2).Position(3);
        MALER.GUI.LabelParam(3).Position(2) = MALER.GUI.LabelParam(2).Position(2) - MALER.GUI.LabelParam(2).Position(4);
       
        MALER.GUI.SliderParam(3).Position(1) = MALER.GUI.SliderParam(2).Position(1) + (1/4)*MALER.GUI.sldparampanel.width;
        MALER.GUI.SliderParam(3).Position(4) = MALER.GUI.SliderParam(2).Position(4) - MALER.GUI.LabelParam(1).Position(4);
        MALER.GUI.SliderParam(3).String = vname;
        MALER.GUI.SliderParam(3).Tag    = '33';
        MALER.GUI.SliderParam(3).Interruptible = 'off';
        MALER.GUI.SliderParam(3).Enable = MALER.GUI.LabelParam(3).Enable;
    end    
end     % end of the function <f_paramsliders_setup>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_assign_sliderparam_values(ii,vname,select,mindrange,numspec)
%%% setup [Value, Min, Max] for slider parameter N ii
    global MALER

    try
        value = MALER.WSP.Param.val{select(ii)}(1);
        if (mindrange)
            vmin = MALER.WSP.Param.val{select(ii)}(2);
            vmax = MALER.WSP.Param.val{select(ii)}(3);
        else
            [value,vmin,vmax] = f_slider_update(value,[]);
        end 
        if (~isnan(value))
            MALER.GUI.SliderParam(ii).Value = value;   MALER.GUI.SliderParam(ii).Min = vmin;   MALER.GUI.SliderParam(ii).Max = vmax;
            MALER.GUI.LabelParam(ii).Enable = 'on';
        else
            nan_or_empty_param_value();            
        end
    catch
        nan_or_empty_param_value();
    end
    
    MALER.GUI.LabelParam(ii).String = [vname,' = ', num2str(value,numspec)];
    
    function nan_or_empty_param_value()
        value =  0;     vmin  = -1;       vmax  = +1;    
        MALER.GUI.SliderParam(ii).Value = value;   MALER.GUI.SliderParam(ii).Min = vmin;   MALER.GUI.SliderParam(ii).Max = vmax;
        MALER.GUI.LabelParam(ii).Enable = 'off';
        MALER.WSP.Param.val{select(ii)}(1) = NaN;
        MALER.WSP.Param.val{select(ii)}(2) = NaN;
        MALER.WSP.Param.val{select(ii)}(3) = NaN;            
    end
    
end   % end of the function <f_assign_sliderparam_values>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_create_param_sliders()
%%% create grid-variable sliders
    global MALER 
    
    if ((isempty(MALER.GUI.Listbox.param.Value))||(MALER.GUI.Listbox.param.Value(1) < 1))
        return,
    end
    
    for i1 = 1:1:min(3,numel(MALER.GUI.Listbox.param.Value)) 
        
        MALER.HDL.SliderParam(i1) = uicontrol('Parent',MALER.HDL.Fig0,'Style','slider');           
        obj = 'SliderParam';
        MALER.HDL.SliderParam(i1) = f_create_gui_object(obj,MALER.HDL.SliderParam(i1),i1);

        MALER.HDL.LabelParam(i1) = uicontrol('Parent',MALER.HDL.Fig0,'Style','text','Visible','on','Enable','on');
        obj = 'LabelParam';
        MALER.HDL.LabelParam(i1) = f_create_gui_object(obj,MALER.HDL.LabelParam(i1),i1);
    end
end     % end of the function <f_create_param_sliders>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotter_update()
%%% update plotter settings
    global MALER    
    if (strcmp(MALER.PPR.type,'slc')),   MALER.PPR.func = 'plot';     end,
    
    if (strcmp(MALER.PPR.type,'slc') || strcmp(MALER.PPR.type,'1d'))
        MALER.GUI.FigButt.ChAbs.Enable  = 'on';
        set(MALER.HDL.FigButt.ChAbs,'Enable',MALER.GUI.FigButt.ChAbs.Enable);
        
        MALER.GUI.FigButt.ChAbs.Visible = 'on';
        set(MALER.HDL.FigButt.ChAbs,'Visible',MALER.GUI.FigButt.ChAbs.Visible);
    else
        MALER.GUI.FigButt.ChAbs.Enable  = 'off';
        set(MALER.HDL.FigButt.ChAbs,'Enable',MALER.GUI.FigButt.ChAbs.Enable);
        
%         MALER.GUI.FigButt.ChAbs.Visible = 'off';
%         set(MALER.HDL.FigButt.ChAbs,'Visible',MALER.GUI.FigButt.ChAbs.Visible);        
    end
end     % end of the function <f_plotter_update>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%   GUI Callback functions   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_fig0_closereqfcn(varargin)
%%% close figure callback
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    try
        delete(MALER.HDL.Fig0),
        MALER = rmfield(MALER,'HDL');       % remove handles
        MALER = rmfield(MALER,'uset');      % remove utility settings
        MALER.GUI.Exist = false;            % clear existance flag
        assignin('base','MALER',MALER);
    catch,     delete(gcf),
    end
end     % end of the function <f_callback_figI_closereqfcn>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_fig0_resizefcn(varargin)
%%% resize figure callback
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    if (exist('MALER','var'))&&(isfield(MALER,'GUI'))&&(isfield(MALER.GUI,'Exist'))&&(MALER.GUI.Exist)    
        MALER.GUI.Fig0.Position = get(MALER.HDL.Fig0,'Position');
        f_plot_selected();
    end
end     % end of the function <f_callback_fig0_resizefcn>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_listbox_param(varargin)
%%% callback to litbox parameters
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    MALER.GUI.Listbox.param.Value = get(MALER.HDL.Listbox.param,'Value');
    f_paramsliders_setup();
    f_create_param_sliders();
end     % end of the function <f_callback_listbox_param>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_listbox_funcs(varargin)
%%% callback to litbox parameters
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
        
    MALER.GUI.Listbox.funcs.Value = get(MALER.HDL.Listbox.funcs,'Value');
    f_plot_selected();
end     % end of the function <f_callback_listbox_funcs>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_gridslider(varargin)
%%% callback function for all grid sliders
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    tag = get(hObject,'Tag');   i1 = eval(tag(1));      i2 = eval(tag(2));
    MALER.GUI.SliderGrid(i1,i2).Value = get(MALER.HDL.SliderGrid(i1,i2),'Value');
    
    % update one particular grid slider and its label
    f_gridslider_update(i1,i2);     % update GUI settings       
    obj = 'SliderGrid';
    MALER.HDL.SliderGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.SliderGrid(i1,i2),i1,i2);
    
    obj = 'LabelGrid';
    MALER.HDL.LabelGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.LabelGrid(i1,i2),i1,i2); 
    
    f_update_grid(i1);      % update computational grid and <GUI.Grid>
    f_eval_all_funcs();     % recalculate all functions
    f_plot_selected();      % plot selected functions   
end     % end of the function <f_callback_gridslider>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_gridslider_update(i1,i2)
%%% update one grid slider and its label
    global MALER
        
    vname = MALER.WSP.Grid.name{i1};
    value = MALER.GUI.SliderGrid(i1,i2).Value;
    vmin  = MALER.GUI.SliderGrid(i1,i2).Min;
    vmax  = MALER.GUI.SliderGrid(i1,i2).Max;
    numspec = MALER.uset.numSpec;
    fsize = MALER.GUI.Settings.fontsize;
    
    range = [vmin,vmax];    
    if ( abs(value-vmin) < 10*eps ),    range    = [];     end,
    if ( abs(value-vmax) < 10*eps ),    range(1) = [];     end,
    
    [value,vmin,vmax] = f_slider_update(value,range,i2);
    MALER.GUI.SliderGrid(i1,i2).Value = value;   MALER.GUI.SliderGrid(i1,i2).Min = vmin;   MALER.GUI.SliderGrid(i1,i2).Max = vmax;
    
    % assign GUI parameters with checking consistency   
    umin = MALER.GUI.SliderGrid(i1,1).Value;   
    du   = MALER.GUI.SliderGrid(i1,2).Value;
    umax = MALER.GUI.SliderGrid(i1,3).Value;
    
    switch int8(i2)
        case 1
            if (  ( umin > umax ) || ( du > (umax-umin) )  )     
                umin = umax - du;
                MALER.GUI.SliderGrid(i1,i2).Value = umin;
                MALER.GUI.SliderGrid(i1,i2).Min = min(MALER.GUI.SliderGrid(i1,i2).Min,umin);
                MALER.GUI.SliderGrid(i1,i2).Max = min(MALER.GUI.SliderGrid(i1,i2).Max,umax);
            end            
            MALER.GUI.LabelGrid(i1,i2).String = ['min ',vname,' = ',num2str( MALER.GUI.SliderGrid(i1,i2).Value, numspec )];
        
        case 2
            if ( du > (umax-umin) )
                du = max(10*eps,abs(umax-umin)/2);
                MALER.GUI.SliderGrid(i1,i2).Value = du;
                MALER.GUI.SliderGrid(i1,i2).Min = max(5*eps,min(MALER.GUI.SliderGrid(i1,i2).Min,du/4));
                MALER.GUI.SliderGrid(i1,i2).Max = max(abs(umax-umin),min(MALER.GUI.SliderGrid(i1,i2).Max,abs(umax-umin)));
                     
            elseif ( du < 10*eps )
                du = max(10*eps,abs(umax-umin)/100);
                MALER.GUI.SliderGrid(i1,i2).Value = du;
                MALER.GUI.SliderGrid(i1,i2).Min = max(10*eps,min(MALER.GUI.SliderGrid(i1,i2).Min,du/4));
                MALER.GUI.SliderGrid(i1,i2).Max = max(abs(umax-umin)/2,min(MALER.GUI.SliderGrid(i1,i2).Max,abs(umax-umin)));                                  
            end    
            MALER.GUI.LabelGrid(i1,i2).String = ['d',   vname,' = ',num2str( MALER.GUI.SliderGrid(i1,i2).Value, numspec )];
        
        case 3
            if (  ( umin > umax ) || ( du > (umax-umin) )  )
                umax = umin + du;
                MALER.GUI.SliderGrid(i1,i2).Value = umax;
                MALER.GUI.SliderGrid(i1,i2).Min = max(MALER.GUI.SliderGrid(i1,i2).Min,umin);
                MALER.GUI.SliderGrid(i1,i2).Max = max(MALER.GUI.SliderGrid(i1,i2).Max,umax);
            end
            MALER.GUI.LabelGrid(i1,i2).String = ['max ',vname,' = ',num2str( MALER.GUI.SliderGrid(i1,i2).Value, numspec )];
        
        otherwise
            disp([mfilename,': ERROR! Unknowm slider number i2 = ',num2str(i2)]);
            h = errordlg(['ERROR! Unknowm slider number i2 = ',num2str(i2)]);
            h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
    end   
end     % end of the function <f_gridslider_update>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_update_grid(i1)
%%% update computational grid and <GUI.Grid>
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    ii = setdiff([1,2],i1);         % ii = <not i1>
    
    umin = MALER.GUI.SliderGrid(i1,1).Value;    
    du   = MALER.GUI.SliderGrid(i1,2).Value;
    umax = MALER.GUI.SliderGrid(i1,3).Value;
    u = umin : du : umax;    
    
    if ( int8(i1) == 2 )            % recalculate y-componet      
        
        x = MALER.WSP.Grid.val{ii,1};                    
        [X,Y] = meshgrid(x,u);
        
        MALER.WSP.Grid.val{ii,2} = X;
        MALER.WSP.Grid.val{i1,1} = u;
        MALER.WSP.Grid.val{i1,2} = Y;
        
    elseif ( int8(i1) == 1 )        % recalculate x - component
        
        MALER.WSP.Grid.val{i1,1} = u;
        if ( MALER.DIM < 2),    return,     end,
        
        y = MALER.WSP.Grid.val{ii,1};
        [X,Y] = meshgrid(u,y);
        
        MALER.WSP.Grid.val{i1,2} = X;
        MALER.WSP.Grid.val{ii,2} = Y;
        
    else
        h = errordlg(['ERROR! Unknowm slider number i1 = ',num2str(i1)]);
        h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
    end

    MALER.WSP.Grid.par{i1} = [umin,umax,length(u)];
    if (ismember('WORKSPACE',evalin('base','who')))
        wsp = evalin('base','WORKSPACE');   field = 'GRIDPARAM';
        wsp.(field){i1} = MALER.WSP.Grid.par{i1};
        assignin('base','WORKSPACE',wsp);
    end
    
    % update grid variables at the slice if exist
    if ( (isfield(MALER.PPR.Slc,'val'))&&(~isempty(MALER.PPR.Slc.val{1}))&&(strcmp(MALER.PPR.Slc.name,MALER.WSP.Grid.name{i1})) )
        f_eval_slice_grid();
    end    
end     % end of the function <f_update_grid>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_pushbutt_setgrid(varargin)
%%% callback function for the button 'Setup Grid'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    delimiters = MALER.uset.inputDelims;
    
    dlg_name = 'SETUP GRID PARAMETERS: VALUE, MIN, MAX';
    dlg_lines  = [1,round(length(dlg_name)*fsize/deffsize)];
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          
    vname = {};
    
    for i1 = 1:1:MALER.DIM 
        vname = [vname, MALER.WSP.Grid.name{i1}];             

        i2 = 1;
        dlg_prompt(3*(i1-1)+i2) = {['setup min(',vname{i1},')']};        
        dlg_defans(3*(i1-1)+i2) = {[ ...
            num2str(MALER.GUI.SliderGrid(i1,i2).Value),',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Min)  ,',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Max) ]};
        
        i2 = 2;
        dlg_prompt(3*(i1-1)+i2) = {['setup grid step d',vname{i1}]};        
        dlg_defans(3*(i1-1)+i2) = {[ ...
            num2str(MALER.GUI.SliderGrid(i1,i2).Value),',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Min)  ,',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Max) ]}; 
        
        i2 = 3;
        dlg_prompt(3*(i1-1)+i2) = {['setup max(',vname{i1},')']};        
        dlg_defans(3*(i1-1)+i2) = {[ ...
            num2str(MALER.GUI.SliderGrid(i1,i2).Value),',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Min)  ,',  ',...
            num2str(MALER.GUI.SliderGrid(i1,i2).Max) ]};             
    end
    
    circleflag = true;
    
    while (circleflag)
        
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch   
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if (isempty(answer)),        return,         end,
        
        circleflag = false;
    
        for j = 1:1:length(answer)
        
            tline = answer{j};        
            qnt = f_extract_keyboard_input(tline,delimiters); % type char
            qnt = f_char2num(qnt);                            % type double                     
           
            i1 = floor( (j-1)/3 ) + 1;
            i2 = mod( (j-1) , 3 ) + 1;
            
            if ( numel(qnt) > 2 )
                value = qnt(1);        vmin = qnt(2);       vmax = qnt(3);
                % check input consistency
                if (( abs(vmax-vmin) < 10*eps ) || ( value < vmin ) || ( value > vmax ))
                    h = warndlg(['Inconsistent parameters in line ',num2str(j),'. Try again']);
                    h = f_resize_dlgbox(h,fsize);       uiwait(h);
                    circleflag = true;                  continue,
                end 
                
            elseif ( numel(qnt) > 1 )     
                value = qnt(1);        vmin  = qnt(2);
                % check input consistency
                if ( value < vmin ) 
                    h = warndlg(['Inconsistent parameters in line ',num2str(j),'. Try again']);
                    h = f_resize_dlgbox(h,fsize);       uiwait(h);
                    circleflag = true;                  continue,
                end
                vmax = MALER.GUI.SliderGrid(i1,i2).Max;
                if (( abs(vmax-vmin) < 10*eps ) || ( value > vmax ))
                    vmax = max(vmin,value)+1;
                end
                
            elseif ( numel(qnt) > 0 )
                value = qnt(1);
                vmin  = MALER.GUI.SliderGrid(i1,i2).Min;
                vmax  = MALER.GUI.SliderGrid(i1,i2).Max;
                if ( value < vmin ),        vmin = value - 1;       end,
                if ( value > vmax ),        vmax = value + 1;       end,
                if ( abs(vmax-vmin) < 10*eps ),   vmax = vmax +1;   end,
                
            else               
                value = MALER.GUI.SliderGrid(i1,i2).Value;
                vmin  = MALER.GUI.SliderGrid(i1,i2).Min;
                vmax  = MALER.GUI.SliderGrid(i1,i2).Max;                             
            end
            
            MALER.GUI.SliderGrid(i1,i2).Value = value;
            MALER.GUI.SliderGrid(i1,i2).Min   = vmin;
            MALER.GUI.SliderGrid(i1,i2).Max   = vmax;
            
            dlg_defans(3*(i1-1)+i2) = {[ ...
                num2str(MALER.GUI.SliderGrid(i1,i2).Value),',  ',...
                num2str(MALER.GUI.SliderGrid(i1,i2).Min)  ,',  ',...
                num2str(MALER.GUI.SliderGrid(i1,i2).Max) ]};   %%#ok<*AGROW>
        end 
        
        % check consistency of the grid parameters and update WSP.Grid.par
        for i1 = 1:MALER.DIM
            umin = MALER.GUI.SliderGrid(i1,1).Value;
            du   = MALER.GUI.SliderGrid(i1,2).Value;
            umax = MALER.GUI.SliderGrid(i1,3).Value;          
            
            if (( (umax-umin) < 10*eps ) || ( du > (umax-umin) ))
                h = warndlg(['Inconsistent parameters of variable ',vname{i1},'. Try again']);
                h = f_resize_dlgbox(h,fsize);       uiwait(h);
                circleflag = true;                  continue,                
            end
            
            % update grid parameters and 'WORKSPACE'
            u = umin:du:umax;
            MALER.WSP.Grid.par{i1} = [umin,umax,length(u)];
            
            if (ismember('WORKSPACE',evalin('base','who')))
                wsp = evalin('base','WORKSPACE');   field = 'GRIDPARAM';
                wsp.(field){i1} = MALER.WSP.Grid.par{i1};
                assignin('base','WORKSPACE',wsp);
            end                         
        end
    end   

    % update grid sliders and their labels
    for i1 = 1:1:MALER.DIM
        for i2 = 1:1:3
            f_gridslider_update(i1,i2);     % update GUI settings       
            obj = 'SliderGrid';
            MALER.HDL.SliderGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.SliderGrid(i1,i2),i1,i2);

            obj = 'LabelGrid';
            MALER.HDL.LabelGrid(i1,i2) = f_create_gui_object(obj,MALER.HDL.LabelGrid(i1,i2),i1,i2); 
        end
    end
    
    % update computational grid
    MALER.WSP.Grid.val = f_create_grid(MALER.WSP.Grid.name,MALER.WSP.Grid.par,0,fsize);
    
    % update grid variables at the slice if exist
    if ( (isfield(MALER.PPR.Slc,'val'))&&(~isempty(MALER.PPR.Slc.val{1})) )
        f_eval_slice_grid;
    end
    
    f_eval_all_funcs();     % recalculate all functions
    f_plot_selected();      % plot selected functions      
end     % end of the function <f_callback_pushbutt_setgrid>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f_char2num(qnt)
%%% convert cell array of chars like this: {'1','1.2'} to vector
    val = NaN*zeros(size(qnt));
    for j = 1:1:length(qnt)
        val(j) = eval(qnt{j});
    end
end     % end of the function <f_char2num>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_paramslider(varargin)
%%% callback function for parameter sliders
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    vname = get(hObject,'String');      % name of parameter
    tag = get(hObject,'Tag');           i1 = eval(tag(2));
    MALER.GUI.SliderParam(i1).Value = get(MALER.HDL.SliderParam(i1),'Value');
    
    % update parameter slider, its label, WSP.Param, and 'WORKSPACE'
    f_paramslider_update(i1);             
    obj = 'SliderParam';
    MALER.HDL.SliderParam(i1) = f_create_gui_object(obj,MALER.HDL.SliderParam(i1),i1);
    
    obj = 'LabelParam';
    MALER.HDL.LabelParam(i1) = f_create_gui_object(obj,MALER.HDL.LabelParam(i1),i1); 
       
    n = f_recalc_depend_func({vname});
    if (n > 0),     f_plot_selected();        end,
    
end     % end of the function <f_callback_paramslider>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_paramslider_update(i1)
%%% update one parameter slider and its label
    global MALER
        
    vname = MALER.GUI.SliderParam(i1).String;    
    value = MALER.GUI.SliderParam(i1).Value;
    vmin  = MALER.GUI.SliderParam(i1).Min;
    vmax  = MALER.GUI.SliderParam(i1).Max;
    numspec = MALER.uset.numSpec;
    fsize = MALER.GUI.Settings.fontsize;
    
    range = [vmin,vmax];    
    if ( abs(value-vmin) < 10*eps ),    range    = [];     end,
    if ( abs(value-vmax) < 10*eps ),    range(1) = [];     end,
    
    [value,vmin,vmax] = f_slider_update(value,range);
    MALER.GUI.SliderParam(i1).Value = value;   MALER.GUI.SliderParam(i1).Min = vmin;   MALER.GUI.SliderParam(i1).Max = vmax;
    
    MALER.GUI.LabelParam(i1).String = [vname,' = ',num2str( MALER.GUI.SliderParam(i1).Value, numspec )];
    
    % update WSP.Param
    parnames = MALER.WSP.Param.name;
    ii = find(strcmp(parnames,vname),1,'first');    % index of current par
    if (~isempty(ii))
        MALER.WSP.Param.val{ii} = [ MALER.GUI.SliderParam(i1).Value,...
            MALER.GUI.SliderParam(i1).Min, MALER.GUI.SliderParam(i1).Max ];
    end
    
    % transfer new value to <WORKSPACE.PARSAMPVAL>
    if (ismember('WORKSPACE',evalin('base','who')))
        wsp = evalin('base','WORKSPACE');       field = 'PARAMETERS'; 
        parnames = wsp.(field);                      % names of parameters
        ii = find(strcmp(parnames,vname),1,'first'); % index of current par
        if (~isempty(ii))
            field = 'PARSAMPVAL';   parvals = wsp.(field);  % param values 
            parvals{ii} = [ MALER.GUI.SliderParam(i1).Value,...
              MALER.GUI.SliderParam(i1).Min, MALER.GUI.SliderParam(i1).Max];             
            wsp.(field) = parvals;                          % reset field
            assignin('base','WORKSPACE',wsp);
        end                     
    end
    
    % update grid variables at the slice if exist
    if ( (isfield(MALER.PPR.Slc,'val'))&&(~isempty(MALER.PPR.Slc.val{1})) )
        f_eval_slice_grid();
    end 
    
    % remove CS center data
    if (isfield(MALER,'CScent')),  MALER = rmfield(MALER,'CScent');    end,
end     % end of the function <f_paramslider_update>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_pushbutt_setparam(varargin)
%%% callback function for the button 'Setup Parameters'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    if (numel(MALER.WSP.Param.name) < 1),   return,     end,
    
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    delimiters = MALER.uset.inputDelims;
    numspec = MALER.uset.numSpec;
    
    dlg_name = 'SETUP PARAMETER VALUES: VALUE, MIN, MAX';
    dlg_lines  = [1,round(length(dlg_name)*fsize/deffsize)];
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          
    vname = {};
    
    for i1 = 1:1:length(MALER.WSP.Param.name) 
        vname = [vname, MALER.WSP.Param.name{i1}];             

        dlg_prompt(i1) = {['setup parameter ',vname{i1}]};        
        dlg_defans(i1) = {[ ...
            num2str(MALER.WSP.Param.val{i1}(1)),',  ',...
            num2str(MALER.WSP.Param.val{i1}(2))  ,',  ',...
            num2str(MALER.WSP.Param.val{i1}(3)) ]};                   
    end
    
    circleflag = true;
    
    while (circleflag)
        
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if (isempty(answer)),        return,         end,
        
        circleflag = false;
        changepar  = {};
        for i1 = 1:1:length(answer)
        
            tline = answer{i1};        
            qnt = f_extract_keyboard_input(tline,delimiters); % type char
            qnt = f_char2num(qnt);                            % type double
            
            if ( numel(qnt) > 2 )
                value = qnt(1);        vmin = qnt(2);       vmax = qnt(3);
                % check input consistency
                if (( abs(vmax-vmin) < 10*eps ) || ( value < vmin ) || ( value > vmax ))
                    h = warndlg({['Inconsistent settings of parameter ',num2str(vname{i1})];'Try again'});
                    h = f_resize_dlgbox(h,fsize);       uiwait(h);
                    circleflag = true;                  continue,
                end 
                
            elseif ( numel(qnt) > 1 )     
                value = qnt(1);        vmin  = qnt(2);
                % check input consistency
                if ( value < vmin ) 
                    h = warndlg({'Inconsistent settings of parameter ',num2str(vname{i1});'Try again'});
                    h = f_resize_dlgbox(h,fsize);       uiwait(h);
                    circleflag = true;                  continue,
                end
                vmax = MALER.WSP.Param.val{i1}(3);
                if (( abs(vmax-vmin) < 10*eps ) || ( value > vmax ))
                    vmax = max(vmin,value)+1;
                end
                
            elseif ( numel(qnt) > 0 )
                value = qnt(1);
                vmin  = MALER.WSP.Param.val{i1}(2);
                vmax  = MALER.WSP.Param.val{i1}(3);
                if ( value < vmin ),        vmin = value - 1;       end,
                if ( value > vmax ),        vmax = value + 1;       end,
                if ( abs(vmax-vmin) < 10*eps ),   vmax = vmax +1;   end,
                
            else               
                value = MALER.WSP.Param.val{i1}(1);
                vmin  = MALER.WSP.Param.val{i1}(2);
                vmax  = MALER.WSP.Param.val{i1}(3);                             
            end
            
            if ( value ~= MALER.WSP.Param.val{i1}(1) )
                changepar = [changepar,MALER.WSP.Param.name{i1}];
            end
            
            MALER.WSP.Param.val{i1}(1) = value;
            MALER.WSP.Param.val{i1}(2) = vmin;
            MALER.WSP.Param.val{i1}(3) = vmax;
            
            dlg_defans(i1) = {[ ...
                num2str(MALER.WSP.Param.val{i1}(1)),',  ',...
                num2str(MALER.WSP.Param.val{i1}(2))  ,',  ',...
                num2str(MALER.WSP.Param.val{i1}(3)) ]};            
        end 
    end  % end of while-loop
    changepar = unique(changepar);

    %%% update grid sliders and their labels
    sldparind = MALER.GUI.Listbox.param.Value;  % indexes of slider parames
    if (~isempty(sldparind))
        sldparnames = {};            % names of slider parameters
        for j = 1:1:length(sldparind)
            sldparnames = [sldparnames,MALER.GUI.SliderParam(j).String];
        end
    end
    
    if (~isempty(sldparnames))
        
        for i1 = 1:1:length(MALER.WSP.Param.name)
            vname = MALER.WSP.Param.name{i1};   ind = [];
            for k = 1:1:length(sldparnames)
                if (~isempty(regexp(sldparnames{k},vname,'once')))
                    ind = k;        break,
                end
            end
            
            if (~isempty(ind))
        
                value = MALER.WSP.Param.val{i1}(1);
                vmin  = MALER.WSP.Param.val{i1}(2);
                vmax  = MALER.WSP.Param.val{i1}(3);               
              
                range = [vmin,vmax];    
                if ( abs(value-vmin) < 10*eps ),    range    = [];     end,
                if ( abs(value-vmax) < 10*eps ),    range(1) = [];     end,
    
                [value,vmin,vmax] = f_slider_update(value,range);
                MALER.GUI.SliderParam(ind).Value = value;   
                MALER.GUI.SliderParam(ind).Min = vmin;   
                MALER.GUI.SliderParam(ind).Max = vmax;
                MALER.GUI.SliderParam(ind).Enable = 'on';
    
                MALER.GUI.LabelParam(ind).String = [vname,' = ',num2str( value, numspec )];
                MALER.GUI.LabelParam(ind).Enable = 'on';
                
                obj = 'SliderParam';
                MALER.HDL.SliderParam(ind) = f_create_gui_object(obj,MALER.HDL.SliderParam(ind),ind);

                obj = 'LabelParam';
                MALER.HDL.LabelParam(ind) = f_create_gui_object(obj,MALER.HDL.LabelParam(ind),ind);               
            end
        end         
    end
         
    %%% update 'WORKSPACE'
    if (ismember('WORKSPACE',evalin('base','who')))
        wsp = evalin('base','WORKSPACE');   field = 'PARSAMPVAL';
        for i1 = 1:1:length(MALER.WSP.Param.name)        
            wsp.(field){i1} = MALER.WSP.Param.val{i1};
        end
        assignin('base','WORKSPACE',wsp);
    end
    
    % update grid variables at the slice if exist
    if ( (isfield(MALER.PPR.Slc,'val'))&&(~isempty(MALER.PPR.Slc.val{1})) )
        f_eval_slice_grid();
    end 
    
    n = f_recalc_depend_func(changepar);
    if (n > 0),       f_plot_selected();        end,
end     % end of the function <f_callback_pushbutt_setparam>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N = f_recalc_depend_func(changepar)
%%% recalculate functions depending on specific parametrs
    global MALER
    N = int8(0);      % number of recalculated functions
    
    % recalculate functions depending on the changed parameters
    if (~MALER.uset.makeSym)
        
        f_eval_all_funcs();         % recalculate all functions
        N = numel(MALER.WSP.Func.name);
        
    else
        % get list of functions depending on changed parameters
        if (isempty(changepar)),    return,         end,
        
        fnc = {''};
        for j = 1:1:length(changepar)
            ind = find(strcmp(MALER.WSP.Param.name,changepar{j}),1);
            if (isempty(ind)),      continue,       end,
            fnc = [fnc,MALER.WSP.Param.dep{ind}];
        end
        fnc(1) = [];    % remove entry ''

        if (isempty(fnc)),         return,          end,
        
        ind = find(ismember(MALER.WSP.Func.name,fnc));
        if (isempty(ind)),         return,          end,
        
        sz = size(ind);
        if (sz(1) > sz(2)),       ind = ind.';      end,            
                
        f_eval_all_funcs(ind);    % recalculate depending functions
        N = length(ind);                             
    end 
    if (isempty(N)),    N = int8(0);        end,
end    % end of the function <f_recalc_depend_func>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_newpar(varargin)              
%%% callback function for the toolbar button 'New Parameter'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    makesym = MALER.uset.makeSym;
    delimiters = MALER.uset.inputDelims;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.NewPar.State = 'off';
    set(MALER.HDL.TbarButt.NewPar,'State',MALER.GUI.TbarButt.NewPar.State);
        
    %%% keyboar input of parameters, create symbolic parameters    
    parnames = f_keyboard_input_params(deffsize,fsize,delimiters,false,true);
    if (isempty(parnames)),     return,     end,
    
    %%% remove variable's names from the new parameter's list
    [C,ia,ib] = intersect(parnames,MALER.WSP.Grid.name);
    if (~isempty(C))
        h = errordlg({['Error! The notation ',parnames{ia(1)},' is occupied by variable.']});
        h = f_resize_dlgbox(h,fsize);   uiwait(h),      return,
    end
 
    %%% remove function's names from the new parameter's list
    [C,ia,ib] = intersect(parnames,MALER.WSP.Func.name);
    if (~isempty(C))
        h = errordlg({['Error! The notation ',parnames{ia(1)},' is occupied by function.']});
        h = f_resize_dlgbox(h,fsize);   uiwait(h),      return,
    end 
    
    %%% clear known parameters if they are listed as new ones
    parnames = sort(parnames);      clearind = [];
    [knownpar,ia,ib] = intersect(parnames,MALER.WSP.Param.name,'stable');
    if (~isempty(knownpar))
        for j = length(knownpar) : - 1 : 1
            choice = questdlg(['Clear parameter ',knownpar{j},'?'], ...
                    'Clear parameter', 'Yes', 'No', 'Yes');
            if strcmp(choice,'Yes')
                clearind = [clearind, ib(j)];
                evalin('base',['clear ',knownpar{j}]);  % clear in base ws
            end
            parnames(ia(j)) = [];    % remove parameter from input list
        end
    end
    if (~isempty(clearind))
        
        % find depending functions
        if (makesym)
            cleardep = {};
            for k = 1:1:length(clearind)
                cleardep = union(cleardep,MALER.WSP.Param.dep{clearind(k)});
            end        
            [C,ia,ib] = intersect(MALER.WSP.Func.name,cleardep,'stable');
        else
            ia = 1:1:numel(MALER.WSP.Func.name);
        end
        
        MALER.WSP.Param.name(clearind) = [];  % clear parameter name
        MALER.WSP.Param.val (clearind) = [];  % clear parameter value        
        MALER.WSP.Param.dep (clearind) = [];  % clear parameter depend-list
        
        if ( numel(ia) > 0 )
            f_eval_all_funcs(ia);   % remove depending functions
            f_plot_selected();
        end
    end
    
    %%% update the list of parameters in MALER.WSP
    if (~isempty(parnames))
        if (makesym),       f_create_symvars(parnames,'base');      end,
        MALER.WSP.Param.name = union(MALER.WSP.Param.name,parnames,'stable');
        MALER.WSP.Param.val(end+1:end+length(parnames)) = {[0, -1, 1]};
        if (~isfield(MALER.WSP.Param,'dep')),   MALER.WSP.Param.dep = {};   end,
        MALER.WSP.Param.dep(end+1:end+length(parnames)) = {''};
    end    
    
    %%% update the list of parameters and functions in WORKSPACE
    if (ismember('WORKSPACE',evalin('base','who')))
        wsp = evalin('base','WORKSPACE');
    end
    field = 'PARAMETERS';   wsp.(field) = MALER.WSP.Param.name;
    field = 'PARSAMPVAL';   wsp.(field) = MALER.WSP.Param.val;
    field = 'FUNCTIONS';
    if (  ( isfield(wsp,field) ) && ( isfield(wsp.(field),'name') )  )
        [C,ia] = setdiff(wsp.(field).name,MALER.WSP.Func.name,'stable');
        if (~isempty(ia))
            finames = fieldnames(wsp.(field));
            for j = 1:1:numel(finames)
                wsp.(field).(finames{j})(ia) = [];
            end
        end
    end
    assignin('base','WORKSPACE',wsp);
        
    %%% update gui listbox <parameters to tune>
    MALER.GUI.Listbox.param.String = MALER.WSP.Param.name;
    MALER.GUI.Listbox.param.Value = min(1,numel(MALER.WSP.Param.name));
    if (MALER.GUI.Listbox.param.Value > 0)
        MALER.GUI.Listbox.param.Enable = 'on';   MALER.GUI.TextParam.Enable = 'on';
        MALER.GUI.Listbox.param.Visible = 'on';  MALER.GUI.TextParam.Visible = 'on';
    else
        MALER.GUI.Listbox.param.Enable = 'off';   MALER.GUI.TextParam.Enable = 'off';
        MALER.GUI.Listbox.param.Visible = 'off';  MALER.GUI.TextParam.Visible = 'off';
    end
    set(MALER.HDL.Listbox.param,'String',MALER.GUI.Listbox.param.String);
    set(MALER.HDL.Listbox.param,'Value',MALER.GUI.Listbox.param.Value);
    set(MALER.HDL.Listbox.param,'Enable',MALER.GUI.Listbox.param.Enable);
    set(MALER.HDL.Listbox.param,'Visible',MALER.GUI.Listbox.param.Visible);
    set(MALER.HDL.TextParam,'Enable',MALER.GUI.TextParam.Enable);
    set(MALER.HDL.TextParam,'Visible',MALER.GUI.TextParam.Visible);
    
    %%% update gui sliders <tune parameters>
    f_paramsliders_setup(1);
    f_create_param_sliders();    
    f_callback_pushbutt_setparam(); 
    
    %%% show figure button <Set parameters>
    if (numel(MALER.WSP.Param.name) > 0)
        MALER.GUI.FigButt.SetParam.Enable = 'on';     
        MALER.GUI.FigButt.SetParam.Visible = 'on';
    else
        MALER.GUI.FigButt.SetParam.Enable = 'off';     
        MALER.GUI.FigButt.SetParam.Visible = 'off';
    end
    set(MALER.HDL.FigButt.SetParam,'Enable',MALER.GUI.FigButt.SetParam.Enable);
    set(MALER.HDL.FigButt.SetParam,'Visible',MALER.GUI.FigButt.SetParam.Visible);
end     % end of the function <f_callback_tbarbutt_newpar>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_newfunc(varargin)             
%%% callback function for the toolbar button 'New Function'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    makesym = MALER.uset.makeSym;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.NewFunc.State = 'off';
    set(MALER.HDL.TbarButt.NewFunc,'State',MALER.GUI.TbarButt.NewFunc.State);
    
    %%% keyboard input of the new function
    newfuncnames = f_keyboard_input_funcs(makesym);
    
    %%% remove deleted functions from the show-list
    MALER.PPR.List.show = intersect(MALER.WSP.Func.name,MALER.PPR.List.show,'stable');
    
    %%% add new functions to the show-list
    MALER.PPR.List.show = union(MALER.PPR.List.show,newfuncnames,'stable');
    MALER.PPR.List.hide = setdiff(MALER.WSP.Func.name,MALER.PPR.List.show,'stable');  
       
    %%% evaluate new functions
    if (~isempty(newfuncnames))
        
        %%% get the list of functions in WORKSPACE
        if (ismember('WORKSPACE',evalin('base','who')))
                wsp = evalin('base','WORKSPACE');
        else    wsp = struct();
        end
        field = 'FUNCTIONS';
        if (~isfield(wsp,field)),           wsp.(field) = struct();    end,
        if (~isfield(wsp.(field),'name')),  wsp.(field).name = {};     end,
        
        for j = 1:1:length(newfuncnames)
            fname = newfuncnames{j};
            ind = find(strcmp(MALER.WSP.Func.name,fname),1,'first');
            
            if (isempty(ind))
                h = warndlg({['Error! Cannot find function ',fname]});
                h = f_resize_dlgbox(h,fsize);   uiwait(h);  continue,
            end
            
            % transfer new/changed function in base workspace
            fnc.name{1} = MALER.WSP.Func.name{ind};  
            fnc.expr{1} = MALER.WSP.Func.expr{ind};
            
            ii = find(strcmp(wsp.(field).name,fnc.name{1}),1,'first');
            if (~isempty(ii))
                wsp.(field).expr{ii} = fnc.expr{1};
            else
                wsp.(field).name{end+1} = fnc.name{1};
                wsp.(field).expr{end+1} = fnc.expr{1};
            end
        
            % correct function expressions (atan-->atan2, operands)
            fnc = f_correct_functions(fnc);
            MALER.WSP.Func.expr{ind} = fnc.expr{1};
            
            % calculate new/updated functions
            if (~makesym),  f_eval_all_funcs();
            else            f_eval_all_funcs(ind);
            end
            f_find_func_dependences(ind);   % find functional dependences            
        end
        f_find_dependlist_param();    % find denedence lists for parameters
        
        assignin('base','WORKSPACE',wsp);
    end
    
    f_update_listbox_funcs();
        
    MALER.GUI.Listbox.absciss.String = union(MALER.WSP.Grid.name,MALER.WSP.Func.name,'stable');
    set(MALER.HDL.Listbox.absciss,'String',MALER.GUI.Listbox.absciss.String);
        
    f_plot_selected();
end     % end of the function <f_callback_tbarbutt_newfunc>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_plottype(varargin)            
%%% callback function for the toolbar button 'Plot Type'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    MALER.GUI.TbarButt.PlotType.State = 'off';
    set(MALER.HDL.TbarButt.PlotType,'State',MALER.GUI.TbarButt.PlotType.State);
    
     %%% create context menu
    if ((~isfield(MALER.HDL,'PlotMenu'))||(isempty( get(MALER.HDL.PlotMenu,'Tag'))))    
        c = uicontextmenu('Tag','PlotMenu');     
        m1 = uimenu(c,'Callback',@f_uimenu_plotype,'Interruptible','off');
        m2 = uimenu(c,'Callback',@f_uimenu_plotype,'Interruptible','off');
        m3 = uimenu(c,'Callback',@f_uimenu_plotype,'Interruptible','off');
        m4 = uimenu(c,'Callback',@f_uimenu_plotype,'Interruptible','off');
        m5 = uimenu(c,'Callback',@f_uimenu_plotype,'Interruptible','off');
        custfsz = f_uifsz(fsize);
        set(m1,'label',['<html><b><font size=',num2str(custfsz),'>contour</font></b></html>']);    
        set(m2,'label',['<html><b><font size=',num2str(custfsz),'>contourf</font></b></html>']);     
        set(m3,'label',['<html><b><font size=',num2str(custfsz),'>mesh</font></b></html>']);
        set(m4,'label',['<html><b><font size=',num2str(custfsz),'>surf</font></b></html>']);
        set(m5,'label',['<html><b><font size=',num2str(custfsz),'>cross-section</font></b></html>']);
        MALER.HDL.PlotMenu = c;                
    end
    set(MALER.HDL.PlotMenu,'Position',f_get_pointer_pos(MALER.HDL.Fig0));
    set(MALER.HDL.PlotMenu,'Visible','on');      

end     % end of the function <f_callback_tbarbutt_plottype>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_uimenu_plotype(varargin)
%%% callback function for uimenu items of the <Plot Type> toolbar button
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER    
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    
    Label = get(hObject,'Label');            % plot function 
    
    if (  ~isempty( regexp(Label,'cross-section','once') )  )
        MALER.PPR.func = 'cross-section';
    elseif (  ~isempty( regexp(Label,'contourf','once') )  )
        MALER.PPR.func = 'contourf';
    elseif (  ~isempty( regexp(Label,'contour','once') )  )
        MALER.PPR.func = 'contour';        
    elseif (  ~isempty( regexp(Label,'meshc','once') )  )
        MALER.PPR.func = 'meshc';
    elseif (  ~isempty( regexp(Label,'mesh','once') )  )
        MALER.PPR.func = 'mesh';        
    elseif (  ~isempty( regexp(Label,'surfl','once') )  )
        MALER.PPR.func = 'surfl'; 
    elseif (  ~isempty( regexp(Label,'surf','once') )  )
        MALER.PPR.func = 'surf'; 
    else
        h = errordlg('Unknown plot function');
        h = f_resize_dlgbox(h,fsize);      uiwait(h);      return,
    end
    
    % evaluate selected plot type 
    if (regexp(MALER.PPR.func,'cross-section'))
        
        % create input dialog box 'cross-section' 
        dlg_name = 'MAKE CROSS-SECTION AT';
        dlg_prompt = 'template: y = x^2';
        dlg_lines  = [1,round( max(length(dlg_name),length(dlg_prompt)) * fsize/deffsize )];     
        dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
        dlg_userfsize = fsize;          dlg_defans = {''};
        if ( (isfield(MALER.PPR.Slc,'lhside')) && (~isempty(MALER.PPR.Slc.lhside)) )
            dlg_defans = {[MALER.PPR.Slc.lhside,' = ',MALER.PPR.Slc.rhside]};
        end
        fname = '';     fstr = '';

        while (isempty(fname))
            try
                answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
            catch
                answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
            end
            if (isempty(answer)),       break,       end,  % input canceled
            if (isempty(answer{1})),    break,       end,  % input canceled

            tline = answer{1};          inqnt = regexp(tline,'=','split');
            if ( length(inqnt) < 2 )           
                h = errordlg({'The equality sign is lost! try again'}); 
                h = f_resize_dlgbox(h,fsize);     uiwait(h);      continue,
            elseif ( length(inqnt) == 2 )   
                fname = strtrim(inqnt{1});      % left-hand side 
                fstr  = strtrim(inqnt{2});      % right-hand side
            else
                h = errordlg({'Too many equality signs! try again'}); 
                h = f_resize_dlgbox(h,fsize);     uiwait(h);      continue,
            end
            
            if (~strcmp(MALER.WSP.Grid.name,fname))
                h = errordlg({'Unknown left-hand side! Try again'}); 
                h = f_resize_dlgbox(h,fsize);     uiwait(h);          
                fname = '';     fstr = '';
                continue,           
            end                     
        end    % end of the <while> loop            
        
        % if input was canceled, then quit        
        if ( (isempty(fname))||(isempty(fstr)) ) 
            
            if ( (~isfield(MALER.PPR,'type'))||(isempty(MALER.PPR.type)) )
                MALER.PPR.type = [num2str(MALER.DIM),'d'];
                if ( MALER.DIM > 1.5 ),     MALER.PPR.func = 'contourf';
                else                        MALER.PPR.func = 'plot';
                end                
            else
                if     (strcmp(MALER.PPR.type,'2d')),   MALER.PPR.func = 'contourf';
                elseif (strcmp(MALER.PPR.type,'1d')),   MALER.PPR.func = 'plot'; 
                end
            end    
            f_plotter_update();
            return,
        end
        
        % correct input expression
        [runstat,outstr] = f_replace_atan(fstr); 
        if (sign(runstat) < 0)
            disp([mfilename,': Warning! Unable to replace <atan> by <atan2> in ',fname]);
            outstr = fstr;
        elseif (~strcmp(fstr,outstr))
            disp([mfilename,' Message: <atan> is replaced by <atan2> in ',fname]);
        end
        fstr = f_matrix_operands(outstr);
        
        ind = regexp(fstr,'\s');
        if (~isempty(ind)),     fstr(ind) = [];        end,
        
        % assign values to PPR.Slc
        MALER.PPR.type = 'slc';
        MALER.PPR.Slc.lhside = fname;
        MALER.PPR.Slc.rhside = fstr;   
        
        % evaluate grid variables at the slice
        f_eval_slice_grid();
        if ( (isempty(MALER.PPR.Slc.val{1})) || (isempty(MALER.PPR.Slc.val{2})) )
            h = errordlg({'Unable evaluate cross-section at ';[fname,' =',fstr]}); 
            h = f_resize_dlgbox(h,fsize);     uiwait(h);
            MALER.PPR.type = [num2str(MALER.DIM),'d'];      return,
        end
        
        % add slice expression in the list of functions
        if (~strcmp(MALER.WSP.Func.name{1},MALER.WSP.Grid.name))
            MALER.WSP.Func.name = union(MALER.PPR.Slc.lhside,MALER.WSP.Func.name,'stable');
            MALER.WSP.Func.expr(2:end+1) = MALER.WSP.Func.expr(1:end);
            MALER.WSP.Func.expr{1} = MALER.PPR.Slc.rhside;
            MALER.WSP.Func.val(2:end+1) = MALER.WSP.Func.val(1:end);
            f_eval_all_funcs(1);    % evaluate slice expression in 2d
            
            MALER.PPR.List.show = union(MALER.PPR.Slc.lhside,MALER.PPR.List.show,'stable');
            
            MALER.GUI.Listbox.funcs.Value = MALER.GUI.Listbox.funcs.Value + 1;
            MALER.GUI.Listbox.funcs.Max   = MALER.GUI.Listbox.funcs.Max + 1;
        else
            MALER.WSP.Func.name{1} = MALER.PPR.Slc.lhside;
            MALER.WSP.Func.expr{1} = MALER.PPR.Slc.rhside;
            MALER.PPR.List.show = union(MALER.PPR.Slc.lhside,MALER.PPR.List.show,'stable');
            f_eval_all_funcs(1);    % evaluate slice expression in 2d
        end     
        
        % reset <listbox functions>
        MALER.GUI.Listbox.funcs.String = MALER.PPR.List.show;
        set(MALER.HDL.Listbox.funcs, 'String', MALER.GUI.Listbox.funcs.String);
        set(MALER.HDL.Listbox.funcs, 'Max'   , MALER.GUI.Listbox.funcs.Max);
        set(MALER.HDL.Listbox.funcs, 'Value' , MALER.GUI.Listbox.funcs.Value);
        
        % reset <choose absciss>
        MALER.GUI.Listbox.absciss.String = union(setdiff(MALER.WSP.Grid.name,MALER.WSP.Func.name),MALER.WSP.Func.name,'stable');
        set(MALER.HDL.Listbox.absciss,'String',MALER.GUI.Listbox.absciss.String);
        
        MALER.PPR.absciss = MALER.PPR.Slc.name;
        
        % evaluate functions at the slice
        f_eval_slice_funcs();
        
    elseif (strcmp(MALER.PPR.func,'plot')),     MALER.PPR.type = '1d';        
    else                                        MALER.PPR.type = '2d';
    end
    
    f_plotter_update();
    f_plot_selected();    
       
end     % end of the function <f_uimenu_plotype>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_eval_slice_grid()
%%% evaluate grid variables at the given cross-section
    global MALER
    MALER.PPR.Slc.val = cell(1,2);
    
    %%% create numerical parameters
    if ( ~isempty(MALER.WSP.Param.name) )&&( ~isempty(MALER.WSP.Param.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Param.name)
            eval( [MALER.WSP.Param.name{y_n_d_e_k_s},' = MALER.WSP.Param.val{y_n_d_e_k_s}(1);'] );
        end
    end
    
    %%% create constants
    if ( ~isempty(MALER.WSP.Const.name) )&&( ~isempty(MALER.WSP.Const.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Const.name)
            eval( [MALER.WSP.Const.name{y_n_d_e_k_s},' = MALER.WSP.Const.val{y_n_d_e_k_s};'] );
        end
    end    
    
    %%% evaluate grid variables
    for y_n_d_e_k_s = 1:1:2    % index of the right-hand side variable
        
        e_n_d_e_k_s = setdiff([1,2],y_n_d_e_k_s);  % index of the left-hand side variable
    
        if (~strcmp(MALER.PPR.Slc.lhside,MALER.WSP.Grid.name{e_n_d_e_k_s}))      
            continue,
        else
            MALER.PPR.Slc.name = MALER.WSP.Grid.name{y_n_d_e_k_s}; 
            eval( [MALER.WSP.Grid.name{y_n_d_e_k_s},' = MALER.WSP.Grid.val{y_n_d_e_k_s,1};'] );           
            MALER.PPR.Slc.val{y_n_d_e_k_s} = MALER.WSP.Grid.val{y_n_d_e_k_s,1};
            try                
                MALER.PPR.Slc.val{e_n_d_e_k_s} = eval([MALER.PPR.Slc.rhside,';']);

                if (~isnumeric(MALER.PPR.Slc.val{e_n_d_e_k_s}(1)))
                    MALER.PPR.Slc.val{e_n_d_e_k_s} = eval(MALER.PPR.Slc.val{e_n_d_e_k_s});
                end
                if (numel(MALER.PPR.Slc.val{e_n_d_e_k_s}) == 1)
                    MALER.PPR.Slc.val{e_n_d_e_k_s} = MALER.PPR.Slc.val{e_n_d_e_k_s}*ones(size(MALER.PPR.Slc.val{y_n_d_e_k_s}));
                end
            catch
                return
            end
        end        
    end
end     % end of the function <f_eval_slice_grid>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_eval_slice_funcs(varargin)
%%% evaluate functions at the given cross-section
    global MALER
    
    if (~isempty(varargin)),    f_u_u_n_k_I_n_D = varargin{1};
    else                        f_u_u_n_k_I_n_D = 1:1:numel(MALER.WSP.Func.name);
    end    
    
    if ( (isempty(MALER.PPR.Slc.val{1})) || (isempty(MALER.PPR.Slc.val{2})) )
        h = errordlg({'Cannot evaluate functions: Cross-section is not defined'}); 
        h = f_resize_dlgbox(h,MALER.GUI.Settings.fontsize);     uiwait(h);
        MALER.PPR.type = [num2str(MALER.DIM),'d'];      return,
    end
    
    %%% create numerical parameters
    if ( ~isempty(MALER.WSP.Param.name) )&&( ~isempty(MALER.WSP.Param.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Param.name)
            eval( [MALER.WSP.Param.name{y_n_d_e_k_s},' = MALER.WSP.Param.val{y_n_d_e_k_s}(1);'] );
        end
    end
    
    %%% create constants
    if ( ~isempty(MALER.WSP.Const.name) )&&( ~isempty(MALER.WSP.Const.name{1}) )
        for y_n_d_e_k_s = 1:1:length(MALER.WSP.Const.name)
            eval( [MALER.WSP.Const.name{y_n_d_e_k_s},' = MALER.WSP.Const.val{y_n_d_e_k_s};'] );
        end
    end
    
    % create numerical grid 
    for y_n_d_e_k_s = 1:1:length(MALER.WSP.Grid.name)
        eval( [MALER.WSP.Grid.name{y_n_d_e_k_s} , ' = MALER.PPR.Slc.val{y_n_d_e_k_s};'] );    
    end 
    
    %%% evaluate functions
    f_u_u_n_k_s = MALER.WSP.Func;             % struct
    f_e_i_l_e_t = [];
    for y_n_d_e_k_s = f_u_u_n_k_I_n_D        
        f_u_u_n_k_c.name = f_u_u_n_k_s.name{y_n_d_e_k_s};
        f_u_u_n_k_c.expr = f_u_u_n_k_s.expr{y_n_d_e_k_s};
        if (~strcmpi(f_u_u_n_k_c.expr(end),';')),  f_u_u_n_k_c.expr = [f_u_u_n_k_c.expr,';'];  end,
        try     
            f_u_u_n_k_v_a_l_u = eval(f_u_u_n_k_c.expr);
            if (~isnumeric(f_u_u_n_k_v_a_l_u(1))),   f_u_u_n_k_v_a_l_u = eval(f_u_u_n_k_v_a_l_u);    end,
            if (numel(f_u_u_n_k_v_a_l_u) < 2),       f_u_u_n_k_v_a_l_u = f_u_u_n_k_v_a_l_u*ones(size(MALER.PPR.Slc.val{1}));    end,
        catch
            e_r_o_r_h_e_n_d_e_l = errordlg({['ERROR: Unable evaluate function ',f_u_u_n_k_c.name]});
            e_r_o_r_h_e_n_d_e_l = f_resize_dlgbox(e_r_o_r_h_e_n_d_e_l,MALER.GUI.Settings.fontsize);  
            uiwait(e_r_o_r_h_e_n_d_e_l);
            f_e_i_l_e_t = [f_e_i_l_e_t , y_n_d_e_k_s];
            f_u_u_n_k_v_a_l_u = NaN*ones(size(MALER.PPR.Slc.val{1}));
        end
        
%         if (strcmpi(f_u_u_n_k_c.name,'theta'))
%             f_u_u_n_k_v_a_l_u = f_u_u_n_k_v_a_l_u + (f_u_u_n_k_v_a_l_u < 0)*2*pi;
%         end

        MALER.PPR.Slc.fname{y_n_d_e_k_s} = f_u_u_n_k_c.name;
        MALER.PPR.Slc.fval {y_n_d_e_k_s} = f_u_u_n_k_v_a_l_u; 
        
        if (~MALER.uset.makeSym)
            eval([MALER.PPR.Slc.fname{y_n_d_e_k_s},' = MALER.PPR.Slc.fval{y_n_d_e_k_s};']);
        end        
    end 
end     % end of the function <f_eval_slice_funcs>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_save(varargin)                     
%%% callback function for the toolbar button 'Save'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.Save.State = 'off';
    set(MALER.HDL.TbarButt.Save,'State',MALER.GUI.TbarButt.Save.State);    
        
    %%% create context menu
    if ((~isfield(MALER.HDL,'SaveMenu'))||(isempty( get(MALER.HDL.SaveMenu,'Tag'))))    
        c = uicontextmenu('Tag','SaveMenu');     
        m1 = uimenu(c,'Callback',@f_uimenu_save_figure,'Interruptible','off');
        m2 = uimenu(c,'Callback',@f_save_problem,'Interruptible','off');
        m3 = uimenu(c,'Callback',@f_save_ws,'Interruptible','off');
        custfsz = f_uifsz(fsize);
        set(m1,'label',['<html><b><i><font size=',num2str(custfsz),'>Figure</font></i></b></html>']);    
        set(m2,'label',['<html><b><i><font size=',num2str(custfsz),'>Problem</font></i></b></html>']);     
        set(m3,'label',['<html><b><i><font size=',num2str(custfsz),'>Workspace</font></i></b></html>']);         
        MALER.HDL.SaveMenu = c;               
    end
    set(MALER.HDL.SaveMenu,'Position',f_get_pointer_pos(MALER.HDL.Fig0));
    set(MALER.HDL.SaveMenu,'Visible','on');
    
    % make prefix and suffix for save-file name
    [pref,suff] = f_make_pref_suff();
    MALER.uset.saveFname.prefix = pref;
    MALER.uset.saveFname.suffix = suff;

end     % end of the function <f_callback_tbarbutt_save>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_uimenu_save_figure(varargin)                     
%%% callback function for the uimenu <Figure> of the toolbar button <Save>
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.Save.State = 'off';
    set(MALER.HDL.TbarButt.Save,'State',MALER.GUI.TbarButt.Save.State);    
        
    %%% create context menu
    if ((~isfield(MALER.HDL,'SaveFigMenu'))||(isempty( get(MALER.HDL.SaveFigMenu,'Tag'))))    
        c = uicontextmenu('Tag','SaveFigMenu');
        m1 = uimenu(c,'Callback',@f_save_fig,'Interruptible','off');
        m2 = uimenu(c,'Callback',@f_specify_save_fig_type,'Interruptible','off');        
        custfsz = f_uifsz(fsize);
        set(m1,'label',['<html><b><i><font size=',num2str(custfsz),'>Save</font></i></b></html>']); 
        set(m2,'label',['<html><b><i><font size=',num2str(custfsz),'>Specify format</font></i></b></html>']);                
        MALER.HDL.SaveFigMenu = c;                
    end
    set(MALER.HDL.SaveFigMenu,'Position',f_get_pointer_pos(MALER.HDL.Fig0)); 
    set(MALER.HDL.SaveFigMenu,'Visible','on');

end     % end of the function <f_uimenu_save_figure>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_specify_save_fig_type(varargin)                     
%%% callback function for the uimenu item <choose figure saving format>
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    figsavetype = MALER.GUI.Settings.figsavetype;
    
    dlg_name   = 'ENTER FILE FORMAT(S)';
    dlg_prompt = {'EXAMPLE: ''eps'', ''pdf'',...'};
    dlg_lines  = [1,round(length(dlg_prompt{1})*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {''};
    
    if (~isempty(figsavetype))
        for j = 1:1:numel(figsavetype)
            dlg_defans{1} = [dlg_defans{1},'''',figsavetype{j},'''',', '];
        end
        dlg_defans{1}(end-1:end) = [];
    end
    
    qnt = [];
    while 1
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end                        
        if (isempty(answer)),       return,     end,  % cancel button press
        
        tline = ['{',answer{1},'}'];
        try    
            qnt = eval(tline);
        catch
            h = errordlg({['Error! Cannot evaluate input:',answer{1}];'Try again'});
            h = f_resize_dlgbox(h,fsize);   uiwait(h);      continue,
        end
        if (~iscell(qnt))
            h = errordlg({'Error! Input expression is not cell array.';'Try again'});
            h = f_resize_dlgbox(h,fsize);   uiwait(h);      continue,
        end
        break,
    end
    
    MALER.GUI.Settings.figsavetype = qnt;     
    f_save_fig();
end     % end of the function <f_specify_save_fig_type>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_save_fig(varargin)             
%%% callback function for the context menu item 'Save Figure'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    figsavetype = MALER.GUI.Settings.figsavetype;
    pref = MALER.uset.saveFname.prefix;
    suff = MALER.uset.saveFname.suffix;
    folder = MALER.uset.outputDir; 
    
    % save initial figure position settings
    fig0units = get(MALER.HDL.Fig0,'Units');
    if (~strcmpi(fig0units,'normalized'))
        set(MALER.HDL.Fig0,'Units','normalized');
    end
    fig0normpos = get(MALER.HDL.Fig0,'Position'); % position in norm. units
        
    %%% create file name: slice specification          
    fname = '';
    if (  (strcmp(MALER.PPR.func,'plot'))&&(MALER.DIM > 1)  )
        fname = ['slc','_',MALER.PPR.Slc.lhside,'=',MALER.PPR.Slc.rhside];
        fname = strrep(fname,'.','');   fname = strrep(fname,'/','_dev_');      
        fname = strrep(fname,'*','');   fname = strrep(fname,'^','_pow_');
        fname = [fname,'__'];
    end
    fname = [fname,MALER.PPR.func];

    % get selected function names
    fselnumbs = get(MALER.HDL.Listbox.funcs,'Value');
    allfnames = get(MALER.HDL.Listbox.funcs,'String');
    fselnames = allfnames(fselnumbs);  
    
    % create file name: list of the plotted functions
    for j = 1:1:length(fselnames),  fname = [fname,'_',fselnames{j}];  end,
    fname = [pref,fname,suff];          % default file name
    
    fname = f_keyboard_input_fname(fname,deffsize,fsize);
    if (isnan(fname)),      return,         end,    
    fname = fullfile(folder,fname);     % user approved full file name 
    %%%--------------------------------------------------------------------
    
    %%% save in .fig format
    nfh = f_make_fig_copy();
    
    set(nfh,'Visible','on');
    try     saveas(nfh,[fname,'.fig'],'fig');
    catch
        h = errordlg('Error! Unable save figure. Check file name.');
        h = f_resize_dlgbox(h,fsize);   uiwait(h);      return,
    end
    
    % the case figsavetype = {}: leave figure open and return
    if (  ( isempty(figsavetype) ) || ( isempty(figsavetype{1}) )  )
        h = msgbox({' ';'Operation completed' ; ['Output file is <',fname,'>']});
        h = f_resize_dlgbox(h,fsize);   uiwait(h);        
        return,
    end
    
    close(nfh),
    figsavetype = setdiff(figsavetype,'fig');
    
    % the case figsavetype = {'fig'}: return
    if (  ( isempty(figsavetype) ) || ( isempty(figsavetype{1}) )  )            
        h = msgbox({' ';'Operation completed' ; ['Output file is <',fname,'>']});
        h = f_resize_dlgbox(h,fsize);   uiwait(h);                
        return,       
    end
    %%%--------------------------------------------------------------------
        
    %%% general case: specify paper orientation
    deforient = 'landscape';          % default paper orientation
    orientation = f_specify_orientation(length(fselnames),MALER.PPR.func,deforient,1);
    
    % resize GUI figure and try to place it in the center of the screen
    if (strcmp(orientation,'landscape')),       FigPos = [0 0 29.7 21];            
    else                                        FigPos = [0 0 21 29.7];             
    end
    
    try
        set(0,'Units','centimeters');  ScrPos = get(0,'ScreenSize');
        
        if (ScrPos(3) > FigPos(3)),    FigPos(1) = (ScrPos(3) - FigPos(3))/2.0;
        else
            coeff = ScrPos(3)/FigPos(3);
            FigPos(3) = coeff*FigPos(3);
            FigPos(4) = coeff*FigPos(4);
        end
        
        if (ScrPos(4) > FigPos(4)),    FigPos(2) = (ScrPos(4) - FigPos(4))/2.0;
        else
            coeff = ScrPos(4)/FigPos(4);
            FigPos(4) = coeff*FigPos(4);
            FigPos(3) = coeff*FigPos(3);
            FigPos(1) = (ScrPos(3) - FigPos(3))/2.0;
        end           
    catch
    end
    if (strcmp(orientation,'landscape')),   FigPos(2) = FigPos(2)-1;   end,
    
    set(MALER.HDL.Fig0,'units','centimeters','Position',FigPos);
    set(MALER.HDL.Fig0,'Visible','off');    
    drawnow,
    nfh = f_make_fig_copy();     % make figure copy     
    
    %%% save figure in non-matlab formats
    status = f_save_figure_nonfig(nfh,fname,orientation,figsavetype,fsize);
    
    if (status == 2)   
        h = msgbox({' ';'Operation completed' ; ['Output file is <',fname,'>']});
        h = f_resize_dlgbox(h,fsize);   uiwait(h);
    elseif (status == 0)
        h = errordlg({' ';'Operation failed';' '});
        h = f_resize_dlgbox(h,fsize);   uiwait(h);
    end
    
    % return Fig0 to initial state
    set(MALER.HDL.Fig0,'Units','normalized','Position',fig0normpos);
    if (~strcmpi(fig0units,'normalized'))
        set(MALER.HDL.Fig0,'Units',fig0units);
    end
    drawnow,     set(MALER.HDL.Fig0,'Visible','on');
end     % end of the function <f_save_fig>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nfh = f_make_fig_copy()
%%% make copy of the gui figure without gui elements
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
       
    %%% find all axes objects in GUI
    GUI_fig_children = get(MALER.HDL.Fig0,'children');
    Fig_Axes = findall(GUI_fig_children,'type','Axes');
    Fig_Cbr  = findall(GUI_fig_children,'type','ColorBar');
    if (isempty(Fig_Axes))&&(isempty(Fig_Cbr)),        return,         end,
    
    %%% create new figure of the A4 size and get axes position
    set(MALER.HDL.Fig0,'Units','centimeters');
    FigPos = get(MALER.HDL.Fig0,'Position');        % is to be A4!
    set(MALER.HDL.Fig0,'Units','Normalized');
    
    nfh = figure('Visible','off','Units','centimeters','Position',FigPos);
    set(nfh,'Visible','off','Units','Normalized', 'toolbar', 'figure');       
    nah = axes('Parent',nfh);   POS = get(nah,'Position');   delete(nah);
    
    % setup image size in normalized units
    p1 = 1;     p2 = 1;     p3 = 0;     p4 = 0;     remind = [];
    
    %%% get position of axes (and colorbar) objects
    N1 = length(Fig_Axes);      N2 = length(Fig_Cbr);
    
    for j = 1:1:N1

        tag = get(Fig_Axes(j),'Tag');
        pos{j,1} = get(Fig_Axes(j),'Position');

        if (strcmpi(tag,'Axe0')),        remind = j;   % do not copy 'Axe0'    
        else
            [p1,p2,p3,p4] = f_get_axes_cbr_min_max_pos(pos{j,1},p1,p2,p3,p4);
        end                                    
    end    
    %%% skip 'Axe0'
    if (~isempty(remind))  
        Fig_Axes(remind) = [];  pos(remind) = [];   N1 = N1 - 1;  
    end
    if (isempty(Fig_Axes))
        disp('Nothing to copy: no any axes object detected.');
        return,     
    end
    
    if (N1 == N2) % make 2d cell array for positions
        
        for j = 1 : 1 : N1
            pos{j,2} = get(Fig_Cbr(j),'Position');
            [p1,p2,p3,p4] = f_get_axes_cbr_min_max_pos(pos{j,2},p1,p2,p3,p4);
        end 
        
    elseif ( N2 > 0)    % continue 1d cell array for positions
        
        for j = N1+1 : 1 : N1+N2
            pos{j} = get(Fig_Cbr(j),'Position');
            [p1,p2,p3,p4] = f_get_axes_cbr_min_max_pos(pos{j},p1,p2,p3,p4);
        end         
    end
    
    %%% calculate image shifts and scales
    dx = (p1-POS(1));           % shift leftward
    dy = (p2-POS(2));           % shift downward
    fx = POS(3) / (p3-p1);      % scaling in the horizontal direction
    fy = POS(4) / (p4-p2);      % scaling in the vertical direction
    
    if (N1 == N2)
        width = 1;
        for j = 1 : 1 : N1           
            new_handle(2*j-1 : 2*j) = copyobj([Fig_Axes(j),Fig_Cbr(j)],nfh);             
            for k = 1:1:2               
                pos{j,k} = f_shift_axes_cbr_pos(pos{j,k},dx,dy,fx,fy,POS);                       
            end
                       
            width = min(width , pos{j,2}(3));  % minimal colorbar width
        end
        for j = 1 : 1 : N1
            pos{j,2}(3) = width;               % correct colorbar width bug
            set(new_handle(2*j-1),'Position',pos{j,1}); 
            set(new_handle(2*j)  ,'Position',pos{j,2});           
        end       
    else       
        new_handle = copyobj(union(Fig_Axes,Fig_Cbr,'stable'),nfh);     
        for j = 1 : 1 : N1+N2
            pos{j} = f_shift_axes_cbr_pos(pos{j},dx,dy,fx,fy,POS);
            set(new_handle(j),'Position',pos{j});
        end
    end 
     
    % correct ylabel rotation if any
    obj = findall(new_handle,'Tag','YLabel');
    if (~isempty(obj)) 
        for k = 1:1:numel(obj),     set(obj(k),'Rotation',0);      end,
    end
    drawnow,
    
    function [p1,p2,p3,p4] = f_get_axes_cbr_min_max_pos(u,p1,p2,p3,p4)
        p1 = min(p1,u(1));                  % min x
        p2 = min(p2,u(2));                  % min y
        p3 = max(p3,u(1)+u(3));             % max x
        p4 = max(p4,u(2)+u(4));             % max y                                                    
    end
    
    function u = f_shift_axes_cbr_pos(u,dx,dy,fx,fy,POS)
        u(1) = u(1) - dx;                   % shift leftward 
        u(2) = u(2) - dy;                   % shift downward           

        u(1) = POS(1) + (u(1)-POS(1))*fx;   % shift rightward
        u(3) = u(3) * fx;                   % scale            

        u(2) = POS(2) + (u(2)-POS(2))*fy;   % shift upward
        u(4) = u(4) * fy;                   % scale        
    end
    
end     % end of the function <f_make_fig_copy>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function orientation = f_specify_orientation(N,plotfunc,deforient,flag)
%%% specify plot orientation for N selected functions
    
    orientation = deforient;          % default paper orientation
    if flag   
        
        if (strcmp(plotfunc,'plot')) || (~isempty(regexp(plotfunc,'surf','once')))       
            orientation = 'landscape'; 
            
        elseif (strcmp(plotfunc,'mesh')) || (~isempty(regexp(plotfunc,'contour','once')))       
            n = 0:1:10;            
            switch N
                case num2cell(n.^2)
                    orientation = 'landscape';                       
                otherwise
                    orientation = 'portrait';
            end       
        else
            disp([mfilename,': Warning! Unknown plot type']);
            orientation = deforient;
        end
    end    
end    % end of the function <f_specify_orientation>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = f_keyboard_input_fname(deffname,deffsize,fsize)
%%% keyboar input of the save file name
    fname = deffname;

    dlg_name = 'ENTER FILE NAME (WITHOUT EXTENSION)';
    dlg_prompt = {'PRESS  CANCEL  FOR  BREAK,  OR  OK  FOR  DEFAULT:'};
    dlg_lines  = [1,round(max(length(dlg_name),length(dlg_prompt{1}))*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = min(18,fsize);          %dlg_defans = {deffname};
    
    if (length(deffname) < length(dlg_prompt{1})), dlg_defans = {deffname};
    else    dlg_defans = {[deffname(1:length(dlg_prompt{1})),' ...']};
    end
    defans_ini = dlg_defans{1};
    
    while 1
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if (isempty(answer)),         fname = NaN;        break,     end,
        
        if (isempty(answer{1})),      fname = deffname;   break,     end,
        
        tline = answer{1};
        if strcmp(tline,defans_ini),  fname = deffname;   break,     end,
        
        ind = regexp(tline, '[/\\*:?"<>|\s]');
        if (~isempty(ind))     
            tline(ind) = [];
            h = errordlg('Error! Enter valid file name.');
            h = f_resize_dlgbox(h,fsize);   uiwait(h);
            dlg_defans{1} = tline;
            continue,
        end
        
        fname = tline;
        break,        
    end
end     % end of the function <f_keyboard_input_fname>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status = f_save_figure_nonfig(fh,fname,ori,format,fsize)
%%% save figure in non - .fig format. Arguments: handle <fh>, file name <fname>, 
%%% paper orientation <ori>, <format> = {'fig',...}, and font size <fsize>
    status = int8(0);     
    
    set(fh,'Visible','off','PaperUnits','default','PaperType','default','PaperSize','default');       
    set(fh,'PaperPositionMode','auto','PaperOrientation',ori);
    
    ext = 'cpbd';   % save in clipboard (Windows only)
    if (  ( ismember(ext,format) ) && ( ispc )  ) 
        try     print(fh,'-dmeta','-r300');            
        catch,  disp('Error! Cannot save figure to clipboard.');
        end            
    end
    format = setdiff(format,ext);
       
    ext = 'emf';    % Windows only
    if (  ( ismember(ext,format) ) && (ispc)  )
        try
            print(fh,'-dmeta','-r300',[fname,'.',ext]);
            format = setdiff(format,ext);
        catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
        end 
    else
        format = setdiff(format,ext);  % remove from the list on non-Windows systems
    end    
    
    ext = 'pdf';    % Windows and MAC
    if (  ( ismember(ext,format) ) && ( ispc || ismac )  ) 
        try
            print(fh,'-dpdf',[fname,'.',ext]);  %saveas(fh,[fname,'.',ext],'pdf');            
            format = setdiff(format,ext);
        catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
        end
    else
        format = setdiff(format,ext);  % remove from the list on Unix systems        
    end
       
    ext = 'eps';
    if ( ismember(ext,format) )             
        try
            print(fh,'-depsc','-r300', [fname,'.',ext]);
            format = setdiff(format,ext);
        catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
        end
    end
    
    ext = 'jpg';    % print in screen resolution
    if ( ismember(ext,format) ) 
        try
            print(fh,'-djpeg',[fname,'.',ext]);
            format = setdiff(format,ext);
        catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
        end            
    end
        
    ext = 'tif';
    if ( ismember(ext,format) )
        try
            print(fh,'-dtiff','-r300', [fname,'.',ext]);
            format = setdiff(format,ext);
        catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
        end
    end
    
    % ext = {'bmp','hdf','jpeg','pbm','pgm','png','ppm','ps','svg','tiff'}
    if ( ~isempty(format) )
        remind = [];
        for j = 1:1:numel(format)
            ext = format{j}; 
            try
                print(fh,['-d',ext],'-r300',[fname,'.',ext]);
                remind = [remind,j];
            catch,  disp(['Error! Cannot save figure in format ','''',ext,'''.']);
            end
        end
        if ( ~isempty(remind) ),       format(remind) = [];       end,
    end
    
    % check if we have saved picture in all specified formats
    if (isempty(format))
        
        status = int8(2);   close(fh), 
        
    else
        status = int8(1);
        set(fh,'Visible','on');
        
        msg = {' '};
        mess = 'Saving in format(s) <';
        for j = 1:1:numel(format)
            mess = [mess,'''',format{j},'''',', '];
        end
        mess(end-1:end) = [];   mess = [mess,'> is not implemented'];
        msg = [msg; mess; ' '; 'Use Matlab saving options'];
                
        h = msgbox(msg);     h = f_resize_dlgbox(h,fsize);     uiwait(h);        
        return
    end
end     % end of the function <f_make_fig_copy>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_save_problem(varargin)
%%% callback function for the context menu item 'Save Problem'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    
    makesym = MALER.uset.makeSym;
    folder  = MALER.uset.outputDir;
    
    fname = 'problem';
    fname = f_keyboard_input_fname(fname,deffsize,fsize);
    if (isnan(fname)),      return,     end,
    
    fname = fullfile(folder,fname);
    
    if (makesym),       fname = [fname,'.m'];
    else                fname = [fname,'.txt'];
    end
    
    fid = fopen(fname,'w');
    if (fid < 0)
        disp([mfilename, ': Error! Cannot open file ',fname]);
        h = errordlg(['Error! Cannot open file ',fname]);
        h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
    end    
    wsp = evalin('base','WORKSPACE');
    vars = wsp.VARIABLES;           varlim = wsp.GRIDPARAM;
    params = wsp.PARAMETERS;        parlim = wsp.PARSAMPVAL;
    if (isfield(wsp,'FUNCTIONS')),  funcs  = wsp.FUNCTIONS;
    else                            funcs  = MALER.WSP.Func;
    end 
    if (isfield(wsp,'CONSTANTS')),  const  = wsp.CONSTANTS;
    else                            const  = [];
    end    
    
    tline = 'WORKSPACE.VARIABLES = {  ';
    for j = 1:1:length(vars),   tline = [tline,'''',vars{j},'''',' , '];   end,
    tline(end-1:end) = [];      tline = [tline,'};          % variables'];
    fprintf(fid,'%s\r\n',tline);
    
    tline = 'WORKSPACE.GRIDPARAM = {  ';
    for j = 1:1:length(varlim), tline = [tline,'[', ...
        num2str(varlim{j}(1)),' , ',num2str(varlim{j}(2)),' , ',num2str(varlim{j}(3)),...
        ']',' , '];       
    end
    tline(end-1:end) = [];  tline = [tline,' };     % domain and the vector length'];
    fprintf(fid,'%s\r\n\n',tline);
    
    tline = 'WORKSPACE.PARAMETERS = {  ';
    for j = 1:1:length(params),   tline = [tline,'''',params{j},'''',' , '];   end,
    tline(end-1:end) = [];        tline = [tline,'};     % parameters'];
    fprintf(fid,'%s\r\n',tline);
    
    tline = 'WORKSPACE.PARSAMPVAL = {  ';
    for j = 1:1:length(parlim),   tline = [tline,num2str(parlim{j}(1)),' , '];   end,
    tline(end-1:end) = [];        tline = [tline,' };     % sampling values of parameters'];
    fprintf(fid,'%s\r\n\n',tline);
    
    if ( ~isempty(const) && ~isempty(const.name) )
        for j = 1:1:length(const.name)
            tline = [const.name{j},' = ',num2str(const.val{j}),';'];
            fprintf(fid,'%s\r\n',tline);
        end
        fprintf(fid,'\n');
    end   
    
    if (makesym)
        tline = 'syms ';
        for j = 1:1:length(vars),   tline = [tline,vars{j},' '];   end,
        tline(end) = ',';
        fprintf(fid,'%s\r\n',tline);
        
        if (~isempty(params))
            tline = 'syms ';
            for j = 1:1:length(params),   tline = [tline,params{j},' '];   end,
            tline(end) = ',';
            fprintf(fid,'%s\r\n\n',tline);            
        end
        
        if (~isempty(funcs.name))
            for j = 1:1:length(funcs.name)
                try
                    tline = [funcs.name{j},' = ',char(evalin('base',funcs.name{j})),';'];
                catch
                    tline = [funcs.name{j},' = ',funcs.expr{j},';'];
                end
                fprintf(fid,'%s\r\n\n',tline);
            end
        end       
        
    else     % make text file        
        if (~isempty(funcs.name))
            for j = 1:1:length(funcs.name)
                tline = [funcs.name{j},' = ','''',funcs.expr{j},'''',';']; 
                fprintf(fid,'%s\r\n\n',tline);
            end
        end        
    end  
    fclose(fid); 
    
    h = msgbox({' ';'Operation Completed' ; ['Output file is <',fname,'>']}); 
    h = f_resize_dlgbox(h,fsize);
end     % end of the function <f_save_problem>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_save_ws(varargin)
%%% callback function for the context menu item 'Save Workspace'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    pref = MALER.uset.saveFname.prefix;     % parameter specification
    if ( numel(pref) > 1 ),     pref(end-1:end) = [];       end,
    
    folder = MALER.uset.outputDir;   
    fname  = 'ws';
    fname  = [fname,'__',pref];
    fname  = f_keyboard_input_fname(fname,deffsize,fsize);
    if (isnan(fname)),      return,     end,
    
    fname  = fullfile(folder,[fname,'.mat']);
    
    if (ismember('WORKSPACE',evalin('base','who')))
        WORKSPACE = evalin('base','WORKSPACE');
    else
        WORKSPACE = [];
    end
    
    M = MALER;    
    MALER = rmfield(MALER,'HDL');       % remove handles
    MALER.GUI.Exist = false;            % clear existance flag
    save(fname,'MALER','WORKSPACE');
    MALER = M;
    
    h = msgbox({' ';'Operation Completed' ; ['Output file is <',fname,'>']}); 
    h = f_resize_dlgbox(h,fsize);    
end     % end of the function <f_save_ws>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pref,suff] = f_make_pref_suff()
%%% make prefix and suffix for file names
    global MALER
    Param = MALER.WSP.Param;
    Grid  = MALER.WSP.Grid;
    nspc  = ['%+.',num2str(MALER.GUI.Settings.dignumber),'f'];
    pref = '';      suff = '';
    
    % make file name prefix    
    if (~isempty(Param.name))
        for j = 1:1:numel(Param.name)
            pref = [pref,Param.name{j},'_',num2str(Param.val{j}(1),nspc),'_'];
        end
        pref = [pref,'_'];
    end
    
    % make file name suffix   
    if (~strcmp(MALER.PPR.type,'slc'))
        for j = 1:1:MALER.DIM
            suff = [suff,'__',Grid.name{j},'_',num2str(Grid.par{j}(1),nspc),'_',...
                num2str(Grid.par{j}(2),nspc),'_',num2str(Grid.par{j}(3))];
        end
    else        
        if     (strcmp(MALER.PPR.Slc.lhside,Grid.name{1}))     
            j = 2;
            suff = [suff,'__',Grid.name{j},'_',num2str(Grid.par{j}(1),nspc),'_',...
                num2str(Grid.par{j}(2),nspc),'_',num2str(Grid.par{j}(3))];
            
        elseif (strcmp(MALER.PPR.Slc.lhside,Grid.name{2}))     
            j = 1;
            suff = [suff,'__',Grid.name{j},'_',num2str(Grid.par{j}(1),nspc),'_',...
                num2str(Grid.par{j}(2),nspc),'_',num2str(Grid.par{j}(3))];            
        end                               
    end    
end    % end of the <f_make_pref_suff> function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_fontsize(varargin)             
%%% callback function for the toolbar button 'FontSize'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.FontSize.State = 'off';
    set(MALER.HDL.TbarButt.FontSize,'State',MALER.GUI.TbarButt.FontSize.State);

    dlg_name = 'CHANGE FONT SIZE';
    dlg_prompt = {'ENTER NEW SIZE'};
    dlg_lines  = [1,round(max(length(dlg_name),length(dlg_prompt{1}))*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {num2str(fsize)};

    newsize = []; 
    while (isempty(newsize))
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if ((isempty(answer))||(isempty(answer{1}))),     break,      end,
        
        try         
            newsize = eval(answer{1});
        catch
            h = errordlg('Error! Cannot evaluate input. Try again');
            h = f_resize_dlgbox(h,fsize);       uiwait(h);      continue,
        end
    end
    if (isempty(newsize)),    return,     end,
    
    newsize = abs(newsize);
    newsize = min(max(newsize,6),32);
    
    %%% reset fontsize in gui handles       
    h1 = findall(MALER.HDL.Fig0,'FontSize',fsize);
    h2 = findall(MALER.HDL.Fig0,'FontSize',fsize+1); 
    for j = 1:1:length(h1),    set(h1(j),'FontSize',newsize);       end,       
    for j = 1:1:length(h2),    set(h2(j),'FontSize',newsize+1);     end,
    
    %%% remove contextmenu tag
    h = findobj('Type','uicontextmenu');
    if (~isempty(h))
        for j = 1:1:length(h)
            Tag = get(h(j),'Tag');
            if (isfield(MALER.HDL,Tag))
                set(MALER.HDL.(Tag),'Tag','');
            end
        end
    end   
    
    %%% reset GUI
    MALER.GUI = f_reset_field(MALER.GUI,'FontSize',false,newsize);
    
    f_plot_selected();      % update plots
    
end     % end of the function <f_callback_tbarbutt_fontsize>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = f_reset_field(S,targetname,casesens,newval)
%%% reset value of the field <targetname> of the structure <S>
%%% <casesens> = {true, false} is the case sensitivity

    if ( (~isstruct(S)) || (isempty(S)) ),     return,     end,
        
    allnames = fieldnames(S);
    if (isempty(allnames)),         return,         end,
    
    if (casesens),      compfunc = @regexp;
    else                compfunc = @regexpi;
    end    
    
    for j = 1:1:length(allnames)
        currentname = allnames{j};      % current field name
        s = S.(currentname);            % current quantity
        if (isstruct(s))
            for k = 1:1:numel(s)
                s(k) = f_reset_field(s(k),targetname,casesens,newval); 
                S.(currentname)(k) = s(k);
            end
        else
            res = compfunc(currentname,targetname);
            if (~isempty(res))
                for k = 1:1:numel(s)
                    S.(currentname)(k) = newval; 
                end
            end
        end       
    end
end     % end of the function <f_reset_field>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_digits(varargin)              %
%%% callback function for the toolbar button 'Digit Number'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    deffsize = MALER.uset.defGraphFsize;
    fsize = MALER.GUI.Settings.fontsize;
    DigNum = MALER.GUI.Settings.dignumber;
    
    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.Digits.State = 'off';
    set(MALER.HDL.TbarButt.Digits,'State',MALER.GUI.TbarButt.Digits.State);

    dlg_name = 'CHANGE LABEL ACCURACY';
    dlg_prompt = {'ENTER NUMBER OF DIGITS'};
    dlg_lines  = [1,round(length(dlg_prompt{1})*fsize/deffsize)];     
    dlg_options.Resize = 'on';      dlg_options.WindowStyle = 'normal';
    dlg_userfsize = fsize;          dlg_defans = {num2str(DigNum)};
    
    dignum = []; 
    while (isempty(dignum))
        try
            answer = maler_indlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options,dlg_userfsize);
        catch
            answer = inputdlg(dlg_prompt,dlg_name,dlg_lines,dlg_defans,dlg_options);
        end
        if ((isempty(answer))||(isempty(answer{1}))),     break,      end,
        
        try         
            dignum = eval(answer{1});
        catch
            h = errordlg('Error! Cannot evaluate input. Try again');
            h = f_resize_dlgbox(h,fsize);       uiwait(h);      continue,
        end
    end
    if (isempty(dignum)),    return,     end, 
    
    dignum = abs(dignum);
    dignum = min(max(dignum,1),12);
    
    %%% reset GUI
    MALER.GUI.Settings.dignumber = dignum;
    MALER.uset.numSpec = ['%.',num2str(MALER.GUI.Settings.dignumber),'g'];
    
    %%% update slider labels
    f_paramsliders_setup();
    f_create_param_sliders();  
    
    obj = 'LabelGrid';
    for i1 = 1:1:MALER.DIM
        for i2 = 1:1:3
            f_gridslider_update(i1,i2);            
            MALER.HDL.(obj)(i1,i2) = ...
                f_create_gui_object(obj,MALER.HDL.(obj)(i1,i2),i1,i2);
        end
    end   
end     % end of the function <f_callback_tbarbutt_digits>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_defaults(varargin)            %
%%% callback function for the toolbar button 'Default Settings'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER

    %%% switch the button visual form to the normal state
    MALER.GUI.TbarButt.Defaults.State = 'off';
    set(MALER.HDL.TbarButt.Defaults,'State',MALER.GUI.TbarButt.Defaults.State);
    
    MALER.GUI.Settings.position = get(MALER.HDL.Fig0,'position');
    f_gui_defaults_update(true);
end     % end of the function <f_callback_tbarbutt_defaults>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_tbarbutt_help(varargin)                %
%%% callback function for the toolbar button 'Help'
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    help_file_name = 'maler_user_guide.pdf';
    help_file_name = fullfile(MALER.GUI.Settings.dir,help_file_name);
    
    if (~exist(MALER.GUI.Settings.dir,'dir'))||(~exist(help_file_name,'file'))
        h = errordlg('Error! Help file is not found');
        h = f_resize_dlgbox(h,fsize);   uiwait(h);
        return,
    end
    
    if (ispc)
        try     open(help_file_name);
        catch
            h = errordlg('Error! Unable open Help file.');
            h = f_resize_dlgbox(h,fsize);   uiwait(h);
            return,        
        end
    else
        h = msgbox(['Open ' help_file_name 'in the operating system'],'Product Help');
        h = f_resize_dlgbox(h,fsize);       uiwait(h);      return,
    end

end     % end of the function <f_callback_tbarbutt_help>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_pushbutt_absciss(varargin)
%%% callback function for the pushbutton <choose absciss>
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    MALER.GUI.Listbox.absciss.String = union(MALER.WSP.Grid.name,MALER.WSP.Func.name,'stable');
    set(MALER.HDL.Listbox.absciss,'String',MALER.GUI.Listbox.absciss.String);   
    
    if (strcmpi(get(MALER.HDL.Listbox.absciss,'Enable'),'off'))
        MALER.GUI.Listbox.absciss.Enable = 'on';            
        MALER.GUI.Listbox.absciss.Visible = 'on';        
    else
        MALER.GUI.Listbox.absciss.Enable = 'off';
        MALER.GUI.Listbox.absciss.Visible = 'off';
    end
    set(MALER.HDL.Listbox.absciss,'Enable',MALER.GUI.Listbox.absciss.Enable);
    set(MALER.HDL.Listbox.absciss,'Visible',MALER.GUI.Listbox.absciss.Visible);
end     % end of the function <f_callback_pushbutt_absciss>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_listbox_absciss(varargin)
%%% callback function for the listbox <choose absciss>
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    
    MALER.GUI.Listbox.absciss.Value = get(MALER.HDL.Listbox.absciss,'Value');
    MALER.PPR.absciss = MALER.GUI.Listbox.absciss.String{MALER.GUI.Listbox.absciss.Value};    
    
    MALER.GUI.Listbox.absciss.Enable = 'off';
    set(MALER.HDL.Listbox.absciss,'Enable',MALER.GUI.Listbox.absciss.Enable);
    MALER.GUI.Listbox.absciss.Visible = 'off';
    set(MALER.HDL.Listbox.absciss,'Visible',MALER.GUI.Listbox.absciss.Visible);
    
    f_plot_selected();
end     % end of the function <f_callback_listbox_absciss>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_textfunc_bdfcn(varargin)
%%% callback function for the text label <plot functions>
    if (~isempty(varargin)),hObject = varargin{1};eventdata = varargin{2};end,
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    %%% create context menu
    if ((~isfield(MALER.HDL,'HideShowMenu'))||(isempty(get(MALER.HDL.HideShowMenu,'Tag'))))    
        c = uicontextmenu('Tag','HideShowMenu');     
        m1 = uimenu(c,'Callback',@f_callback_showhide,'Interruptible','off');
        m2 = uimenu(c,'Callback',@f_callback_showhide,'Interruptible','off');
        custfsz = f_uifsz(fsize);
        set(m1,'label',['<html><b><i><font size=',num2str(custfsz),'>Hide</font></i></b></html>']);    
        set(m2,'label',['<html><b><i><font size=',num2str(custfsz),'>Show</font></i></b></html>']);     
        MALER.HDL.HideShowMenu = c;
    end
    set(MALER.HDL.HideShowMenu,'Position',f_get_pointer_pos(MALER.HDL.Fig0));
    set(MALER.HDL.HideShowMenu,'Visible','on');  
end     % end of the function <f_callback_textfunc_bdfcn>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_showhide(source,callbackdata)               %#ok<INUSD>
%%% callback function for the <show-hide> context menu
    global MALER
    fsize = MALER.GUI.Settings.fontsize;    
    Label = get(source,'Label');
    
    if (~isempty(regexp(lower(Label),'hide','once')))
        
        if (numel(MALER.PPR.List.show) < 2),    return,     end,
        
        MALER.GUI.Listbox.hideshow.String = union(MALER.PPR.List.show,'cancel','stable');
        MALER.GUI.Listbox.hideshow.Min = 0;
        MALER.GUI.Listbox.hideshow.Max = numel(MALER.PPR.List.show)-1; 
        MALER.GUI.Listbox.hideshow.Value = 1;
        MALER.GUI.Listbox.hideshow.Tag = 'hide';
        
    elseif (~isempty(regexp(lower(Label),'show','once')))
        
        if (numel(MALER.PPR.List.hide) < 1),    return,     end,
        
        MALER.GUI.Listbox.hideshow.String = union('all',MALER.PPR.List.hide,'stable');
        MALER.GUI.Listbox.hideshow.String = union(MALER.GUI.Listbox.hideshow.String,'cancel','stable');
        MALER.GUI.Listbox.hideshow.Min = 0;
        MALER.GUI.Listbox.hideshow.Max = numel(MALER.PPR.List.hide)+1;
        MALER.GUI.Listbox.hideshow.Value = 1;
        MALER.GUI.Listbox.hideshow.Tag = 'show';
        
    else
        h = errordlg('Error! Unknown context menu entry');
        h = f_resize_dlgbox(h,fsize);   uiwait(h);  return,
    end 
    
    set(MALER.HDL.Listbox.hideshow,'String',MALER.GUI.Listbox.hideshow.String);
    set(MALER.HDL.Listbox.hideshow,'Min',MALER.GUI.Listbox.hideshow.Min);
    set(MALER.HDL.Listbox.hideshow,'Max',MALER.GUI.Listbox.hideshow.Max);
    set(MALER.HDL.Listbox.hideshow,'Value',MALER.GUI.Listbox.hideshow.Value);
    set(MALER.HDL.Listbox.hideshow,'Tag',MALER.GUI.Listbox.hideshow.Tag);
    set(MALER.HDL.Listbox.hideshow,'Enable','on','Visible','on');
    
    set(MALER.HDL.TextFuncs,'Enable','off','Visible','off');    
    set(MALER.HDL.FigButt.Apply,'Enable','on','Visible','on');    
end     % end of the function <f_callback_showhide>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_listbox_hideshow(varargin)
%%% callback for the listbox <hide/show functions>
    global MALER
    fsize = MALER.GUI.Settings.fontsize;
    
    % get listbox parameters
    MALER.GUI.Listbox.hideshow.Tag = get(MALER.HDL.Listbox.hideshow,'Tag');             
    MALER.GUI.Listbox.hideshow.Value = get(MALER.HDL.Listbox.hideshow,'Value');
    MALER.GUI.Listbox.hideshow.String = get(MALER.HDL.Listbox.hideshow,'String');
    
    % if <cancel> selected then break
    if (MALER.GUI.Listbox.hideshow.Value(end) == numel(MALER.GUI.Listbox.hideshow.String))  
        set(MALER.HDL.Listbox.hideshow,'Enable','off','Visible','off');
        set(MALER.HDL.FigButt.Apply,'Enable','off','Visible','off');
        set(MALER.HDL.TextFuncs,'Enable','on','Visible','on');       
        return,
    end
    
    % if <all> selected, apply it immediately
    if strcmp(MALER.GUI.Listbox.hideshow.String{MALER.GUI.Listbox.hideshow.Value(1)},'all')
        set(MALER.HDL.Listbox.hideshow,'Enable','off','Visible','off');
        set(MALER.HDL.FigButt.Apply,'Enable','off','Visible','off');
        set(MALER.HDL.TextFuncs,'Enable','on','Visible','on');
        
        MALER.PPR.List.show = MALER.WSP.Func.name;
        MALER.PPR.List.hide = setdiff(MALER.WSP.Func.name,MALER.PPR.List.show,'stable');
        f_update_listbox_funcs();
        f_plot_selected();        
        return,        
    end    
end     % end of the function <f_callback_listbox_hideshow>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_callback_figbutt_apply(varargin)
%%% callback for the figure button  <Apply>
    global MALER
    
    set(MALER.HDL.Listbox.hideshow,'Enable','off','Visible','off');
    set(MALER.HDL.FigButt.Apply,'Enable','off','Visible','off');
    set(MALER.HDL.TextFuncs,'Enable','on','Visible','on');
    
    Tag = MALER.GUI.Listbox.hideshow.Tag;
    Value = MALER.GUI.Listbox.hideshow.Value;
    String = MALER.GUI.Listbox.hideshow.String;
    
    switch lower(Tag)
        case 'hide'

            Hidelist = String(Value);
            MALER.PPR.List.show = setdiff(MALER.PPR.List.show,Hidelist,'stable');           
            
        case 'show'
            Showlist = String(Value);
            if (strcmp(Showlist{1},'all'))
                MALER.PPR.List.show = MALER.WSP.Func.name;                
            else                
                MALER.PPR.List.show = union(MALER.PPR.List.show,Showlist,'stable');               
            end
            
        otherwise
            h = errordlg(['Error! Unknown action ',Tag]);
            h = f_resize_dlgbox(h,fsize);   uiwait(h);      return,
    end
    MALER.PPR.List.hide = setdiff(MALER.WSP.Func.name,MALER.PPR.List.show,'stable');
    
    f_update_listbox_funcs();
    f_plot_selected();    
end     % end of the function <f_callback_figbutt_apply>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function custfsz = f_uifsz(fsize)
%%% define font size for uimenu
    custfsz = min( max( int8(sqrt(fsize)+1.3) , 4 ) , 7 );
end    % end of the function <f_uifsz>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = f_get_pointer_pos(fh)
%%% get pointer position respective to the figure
    set(0,'Units','pixels');   pos0 = get(0,'PointerLocation');   

    un = get(fh,'Units');
    set(fh,'Units','pixels');   pos1 = get(fh,'Position'); 
    set(fh,'Units',un);
               
    pos = [pos0(1) - pos1(1) , pos0(2) - pos1(2)];
end     % end of the function <f_get_pointer_pos>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = f_resize_dlgbox(h,fsize)
%%% resize errordlg and message dialog boxes
    fsize = max(6,min(fsize,20));        % limiter
    htext = findobj(h, 'Type', 'Text');  %find text control in dialog
    set(htext,'FontSize',fsize);
    set(h,'Resize','on');           pos = get(h,'Position');
    try
        deffsize = get(0,'factoryUicontrolFontSize');
    catch
        deffsize = 8;
    end
    set(h,'Position',[pos(1)-(pos(3)*(fsize/deffsize-1)/2), pos(2), pos(3)*fsize/deffsize, pos(4)]);           
end    % end of the function <f_resize_dlgbox>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start the function group <maler_indlg>
function Answer = maler_indlg(Prompt, Title, NumLines, DefAns, Resize, UserFontSize)
if (nargin < 6), UserFontSize = get(0,'FactoryUicontrolFontSize'); end,
UserFontSize = max(6,min(UserFontSize,20));    % limiter
%  ANSWER = INPUTDLG(PROMPT) creates a modal dialog box that returns user
%  input for multiple prompts in the cell array ANSWER. PROMPT is a cell
%  array containing the PROMPT strings.
%
%  INPUTDLG uses UIWAIT to suspend execution until the user responds.
%
%  ANSWER = INPUTDLG(PROMPT,NAME) specifies the title for the dialog.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES) specifies the number of lines for
%  each answer in NUMLINES. NUMLINES may be a constant value or a column
%  vector having one element per PROMPT that specifies how many lines per
%  input field. NUMLINES may also be a matrix where the first column
%  specifies how many rows for the input field and the second column
%  specifies how many columns wide the input field should be.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER) specifies the
%  default answer to display for each PROMPT. DEFAULTANSWER must contain
%  the same number of elements as PROMPT and must be a cell array of
%  strings.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER,OPTIONS) specifies
%  additional options. If OPTIONS is the string 'on', the dialog is made
%  resizable. If OPTIONS is a structure, the fields Resize, WindowStyle, and
%  Interpreter are recognized. Resize can be either 'on' or
%  'off'. WindowStyle can be either 'normal' or 'modal'. Interpreter can be
%  either 'none' or 'tex'. If Interpreter is 'tex', the prompt strings are
%  rendered using LaTeX.
%
%  Examples:
%
%  prompt={'Enter the matrix size for x^2:','Enter the colormap name:'};
%  name='Input for Peaks function';
%  numlines=1;
%  defaultanswer={'20','hsv'};
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer);
%
%  options.Resize='on';
%  options.WindowStyle='normal';
%  options.Interpreter='tex';
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%
%  See also DIALOG, ERRORDLG, HELPDLG, LISTDLG, MSGBOX,
%    QUESTDLG, TEXTWRAP, UIWAIT, WARNDLG .
 
%  Copyright 1994-2010 The MathWorks, Inc.
%  $Revision: 1.58.4.21 $
 
%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
try     error(nargchk(0,6,nargin));             %#ok<NCHKN>
catch,  narginchk(0,6);
end
try     error(nargoutchk(0,1,nargout));         %#ok<NCHKE>
catch,  nargoutchk(0,1);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle Input Args %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
  Prompt=getString(message('MATLAB:uistring:popupdialogs:InputDlgInput'));
end
if ~iscell(Prompt)
  Prompt={Prompt};
end
NumQuest=numel(Prompt);
 
 
if nargin<2,
  Title=' ';
end
 
if nargin<3
  NumLines=1;
end
 
if nargin<4
  DefAns=cell(NumQuest,1);
  for lp=1:NumQuest
    DefAns{lp}='';
  end
end
 
if nargin<5
  Resize = 'off';
end
WindowStyle='modal';
Interpreter='none';
 
Options = struct([]); %%#ok
if nargin>4 && isstruct(Resize)
  Options = Resize;
  Resize  = 'off';
  if isfield(Options,'Resize'),      Resize=Options.Resize;           end
  if isfield(Options,'WindowStyle'), WindowStyle=Options.WindowStyle; end
  if isfield(Options,'Interpreter'), Interpreter=Options.Interpreter; end
end
 
[rw,cl]=size(NumLines);
OneVect = ones(NumQuest,1);
if (rw == 1 & cl == 2) %#ok Handle []
  NumLines=NumLines(OneVect,:);
elseif (rw == 1 & cl == 1) %#ok
  NumLines=NumLines(OneVect);
elseif (rw == 1 & cl == NumQuest) %#ok
  NumLines = NumLines';
elseif (rw ~= NumQuest | cl > 2) %#ok
  error(message('MATLAB:inputdlg:IncorrectSize'))
end
 
if ~iscell(DefAns),
  error(message('MATLAB:inputdlg:InvalidDefaultAnswer'));
end
 
%%%%%%%%%%%%%%%%%%%%%%%
%%% Create InputFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigWidth=175;
FigHeight=100;
FigPos(3:4)=[FigWidth FigHeight];  %%#ok
FigColor=get(0,'DefaultUicontrolBackgroundColor');
 
InputFig=dialog(                     ...
  'Visible'          ,'off'      , ...
  'KeyPressFcn'      ,@doFigureKeyPress, ...
  'Name'             ,Title      , ...
  'Pointer'          ,'arrow'    , ...
  'Units'            ,'pixels'   , ...
  'UserData'         ,'Cancel'   , ...
  'Tag'              ,Title      , ...
  'HandleVisibility' ,'callback' , ...
  'Color'            ,FigColor   , ...
  'NextPlot'         ,'add'      , ...
  'WindowStyle'      ,WindowStyle, ...
  'Resize'           ,Resize       ...
  );
 
 
%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset    = 5;
DefBtnWidth  = 53;
DefBtnHeight = 23;
 
TextInfo.Units              = 'pixels'   ;
TextInfo.FontSize           = UserFontSize; % get(0,'FactoryUicontrolFontSize');
TextInfo.FontWeight         = get(InputFig,'DefaultTextFontWeight');
TextInfo.HorizontalAlignment= 'left'     ;
TextInfo.HandleVisibility   = 'callback' ;
 
StInfo=TextInfo;
StInfo.Style              = 'text'  ;
StInfo.BackgroundColor    = FigColor;
 
 
EdInfo=StInfo;
EdInfo.FontWeight      = get(InputFig,'DefaultUicontrolFontWeight');
EdInfo.Style           = 'edit';
EdInfo.BackgroundColor = 'white';
 
BtnInfo=StInfo;
BtnInfo.FontWeight          = get(InputFig,'DefaultUicontrolFontWeight');
BtnInfo.Style               = 'pushbutton';
BtnInfo.HorizontalAlignment = 'center';
 
% Add VerticalAlignment here as it is not applicable to the above.
TextInfo.VerticalAlignment  = 'bottom';
TextInfo.Color              = get(0,'FactoryUicontrolForegroundColor');
 
 
% adjust button height and width
btnMargin=1.4;
ExtControl=uicontrol(InputFig   ,BtnInfo     , ...
  'String'   ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))        , ...
  'Visible'  ,'off'         ...
  );
 
% BtnYOffset  = DefOffset;
BtnExtent = get(ExtControl,'Extent');
BtnWidth  = max(DefBtnWidth,BtnExtent(3)+8);
BtnHeight = max(DefBtnHeight,BtnExtent(4)*btnMargin);
delete(ExtControl);
 
% Determine # of lines for all Prompts
TxtWidth=FigWidth-2*DefOffset;
ExtControl=uicontrol(InputFig   ,StInfo     , ...
  'String'   ,''         , ...
  'Position' ,[ DefOffset DefOffset 0.96*TxtWidth BtnHeight ] , ...
  'Visible'  ,'off'        ...
  );
 
WrapQuest=cell(NumQuest,1);
QuestPos=zeros(NumQuest,4);
 
for ExtLp=1:NumQuest
  if size(NumLines,2)==2
    [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
      textwrap(ExtControl,Prompt(ExtLp),NumLines(ExtLp,2));
  else
    [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
      textwrap(ExtControl,Prompt(ExtLp),80);
  end
end % for ExtLp
 
delete(ExtControl);
QuestWidth =QuestPos(:,3);
QuestHeight=QuestPos(:,4);
if ismac % Change Edit box height to avoid clipping on mac.
    editBoxHeightScalingFactor = 1.4;
else 
    editBoxHeightScalingFactor = 1;
end
TxtHeight=QuestHeight(1)/size(WrapQuest{1,1},1) * editBoxHeightScalingFactor;
EditHeight=TxtHeight*NumLines(:,1);
EditHeight(NumLines(:,1)==1)=EditHeight(NumLines(:,1)==1)+4;
 
FigHeight=(NumQuest+2)*DefOffset    + ...
  BtnHeight+sum(EditHeight) + ...
  sum(QuestHeight);
 
TxtXOffset=DefOffset;
 
QuestYOffset=zeros(NumQuest,1);
EditYOffset=zeros(NumQuest,1);
QuestYOffset(1)=FigHeight-DefOffset-QuestHeight(1);
EditYOffset(1)=QuestYOffset(1)-EditHeight(1);
 
for YOffLp=2:NumQuest,
  QuestYOffset(YOffLp)=EditYOffset(YOffLp-1)-QuestHeight(YOffLp)-DefOffset;
  EditYOffset(YOffLp)=QuestYOffset(YOffLp)-EditHeight(YOffLp);
end % for YOffLp
 
QuestHandle=[];
EditHandle=[];
 
AxesHandle=axes('Parent',InputFig,'Position',[0 0 1 1],'Visible','off');
 
inputWidthSpecified = false;
 
for lp=1:NumQuest,
  if ~ischar(DefAns{lp}),
    delete(InputFig);
    error(message('MATLAB:inputdlg:InvalidInput'));
  end
 
  
  EditHandle(lp)=uicontrol(InputFig    , ...
    EdInfo      , ...
    'Max'        ,NumLines(lp,1)       , ...
    'Position'   ,[ TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp)], ...
    'String'     ,DefAns{lp}           , ...
    'Tag'        ,'Edit'                 ...
    );
 
  QuestHandle(lp)=text('Parent'     ,AxesHandle, ...
    TextInfo     , ...
    'Position'   ,[ TxtXOffset QuestYOffset(lp)], ...
    'String'     ,WrapQuest{lp}                 , ...
    'Interpreter',Interpreter                   , ...
    'Tag'        ,'Quest'                         ...
    );
 
  MinWidth = max(QuestWidth(:));
  if (size(NumLines,2) == 2)
    % input field width has been specified.
    inputWidthSpecified = true;
    EditWidth = setcolumnwidth(EditHandle(lp), NumLines(lp,1), NumLines(lp,2));
    MinWidth = max(MinWidth, EditWidth);
  end
  FigWidth=max(FigWidth, MinWidth+2*DefOffset);
 
end % for lp
 
% fig width may have changed, update the edit fields if they dont have user specified widths.
if ~inputWidthSpecified
  TxtWidth=FigWidth-2*DefOffset;
  for lp=1:NumQuest
    set(EditHandle(lp), 'Position', [TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp)]);
  end
end
 
FigPos=get(InputFig,'Position');
 
FigWidth=max(FigWidth,2*(BtnWidth+DefOffset)+DefOffset);
FigPos(1)=0;
FigPos(2)=0;
FigPos(3)=FigWidth;
FigPos(4)=FigHeight;
 
set(InputFig,'Position',getnicedialoglocation(FigPos,get(InputFig,'Units')));
 
OKHandle=uicontrol(InputFig     ,              ...
  BtnInfo      , ...
  'Position'   ,[ FigWidth-2*BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
  'KeyPressFcn',@doControlKeyPress , ...
  'String'     ,getString(message('MATLAB:uistring:popupdialogs:OK'))        , ...
  'Callback'   ,@doCallback , ...
  'Tag'        ,'OK'        , ...
  'UserData'   ,'OK'          ...
  );
 
setdefaultbutton(InputFig, OKHandle);
 
CancelHandle=uicontrol(InputFig     ,              ...
  BtnInfo      , ...
  'Position'   ,[ FigWidth-BtnWidth-DefOffset DefOffset BtnWidth BtnHeight ]           , ...
  'KeyPressFcn',@doControlKeyPress            , ...
  'String'     ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))    , ...
  'Callback'   ,@doCallback , ...
  'Tag'        ,'Cancel'    , ...
  'UserData'   ,'Cancel'       ...
  ); %%#ok
 
handles = guihandles(InputFig);
handles.MinFigWidth = FigWidth;
handles.FigHeight   = FigHeight;
handles.TextMargin  = 2*DefOffset;
guidata(InputFig,handles);
set(InputFig,'ResizeFcn', {@doResize, inputWidthSpecified});
 
% make sure we are on screen
movegui(InputFig)
 
% if there is a figure out there and it's modal, we need to be modal too
if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
  set(InputFig,'WindowStyle','modal');
end
 
set(InputFig,'Visible','on');
drawnow;
 
if ~isempty(EditHandle)
  uicontrol(EditHandle(1));
end
 
if ishghandle(InputFig)
  % Go into uiwait if the figure handle is still valid.
  % This is mostly the case during regular use.
  uiwait(InputFig);
end
 
% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(InputFig)
  Answer={};
  if strcmp(get(InputFig,'UserData'),'OK'),
    Answer=cell(NumQuest,1);
    for lp=1:NumQuest,
      Answer(lp)=get(EditHandle(lp),{'String'});
    end
  end
  delete(InputFig);
else
  Answer={};
end
end
%%%--------------

function doFigureKeyPress(obj, evd) %#ok
switch(evd.Key)
  case {'return','space'}
    set(gcbf,'UserData','OK');
    uiresume(gcbf);
  case {'escape'}
    delete(gcbf);
end
end
%%%--------------
 
function doControlKeyPress(obj, evd) %%#ok 
switch(evd.Key)
  case {'return'}
    if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
    else
      delete(gcbf)
    end
  case 'escape'
    delete(gcbf)
end
end
%%%--------------
 
function doCallback(obj, evd) %#ok
if ~strcmp(get(obj,'UserData'),'Cancel')
  set(gcbf,'UserData','OK');
  uiresume(gcbf);
else
  delete(gcbf)
end
end
%%%--------------
 
function doResize(FigHandle, evd, multicolumn) %#ok
% TBD: Check difference in behavior w/ R13. May need to implement
% additional resize behavior/clean up.
 
Data=guidata(FigHandle);
 
resetPos = false;
 
FigPos = get(FigHandle,'Position');
FigWidth = FigPos(3);
FigHeight = FigPos(4);
 
if FigWidth < Data.MinFigWidth
  FigWidth  = Data.MinFigWidth;
  FigPos(3) = Data.MinFigWidth;
  resetPos = true;
end
 
% make sure edit fields use all available space if
% number of columns is not specified in dialog creation.
if ~multicolumn
  for lp = 1:length(Data.Edit)
    EditPos = get(Data.Edit(lp),'Position');
    EditPos(3) = FigWidth - Data.TextMargin;
    set(Data.Edit(lp),'Position',EditPos);
  end
end
 
if FigHeight ~= Data.FigHeight
  FigPos(4) = Data.FigHeight;
  resetPos = true;
end
 
if resetPos
  set(FigHandle,'Position',FigPos);
end
end
%%%--------------
 
% set pixel width given the number of columns
function EditWidth = setcolumnwidth(object, rows, cols)
% Save current Units and String.
old_units = get(object, 'Units');
old_string = get(object, 'String');
old_position = get(object, 'Position');
 
set(object, 'Units', 'pixels')
set(object, 'String', char(ones(1,cols)*'x'));
 
new_extent = get(object,'Extent');
if (rows > 1)
  % For multiple rows, allow space for the scrollbar
  new_extent = new_extent + 19; % Width of the scrollbar
end
new_position = old_position;
new_position(3) = new_extent(3) + 1;
set(object, 'Position', new_position);
 
% reset string and units
set(object, 'String', old_string, 'Units', old_units);
EditWidth = new_extent(3);
end
%%%--------------
 
function figure_size = getnicedialoglocation(figure_size, figure_units)
parentHandle = gcbf;
propName = 'Position';
if isempty(parentHandle),   parentHandle = 0;   propName = 'ScreenSize';    end,
old_u = get(parentHandle,'Units');
set(parentHandle,'Units',figure_units);
container_size=get(parentHandle,propName);
set(parentHandle,'Units',old_u);
figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));
end
%%%--------------
 
function setdefaultbutton(figHandle, btnHandle)
if nargin<1, error('MATLAB:setdefaultbutton:InvalidNumberOfArguments','Too few arguments for setdefaultbutton');  end,
if nargin>2, error('MATLAB:setdefaultbutton:InvalidNumberOfArguments','Too many arguments for setdefaultbutton'); end,
if (usejava('awt') == 1)
    useJavaDefaultButton(figHandle, btnHandle);
else
    useHGDefaultButton(figHandle, btnHandle);
end
end
%%%--------------
 
function useJavaDefaultButton(figH, btnH)
fh = handle(figH);
fh.setDefaultButton(btnH);
end
%%%--------------
 
function useHGDefaultButton(figHandle, btnHandle) %#ok<INUSL>
btnPos = getpixelposition(btnHandle);
leftOffset   = btnPos(1) - 1;
bottomOffset = btnPos(2) - 2;
widthOffset  = btnPos(3) + 3;
heightOffset = btnPos(4) + 3;
h1 = uipanel(get(btnHandle, 'Parent'), 'HighlightColor', 'black', ...
'BorderType', 'etchedout','units', 'pixels', ...
'Position', [leftOffset bottomOffset widthOffset heightOffset]);
uistack(h1, 'bottom');
end
%%% end of the function group <maler_indlg>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%     T H E   E N D     %%%%%%%%%% 15.06.2018 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%