function maler_userbutt_callbacks(varargin)
%%% callbacks for the user-defined buttons
 
if ( numel(varargin) < 1 ),
    h = errordlg('ERROR! CANNOT GET BUTTON HANDLE!');
    h = f_resize_dlgbox(h,14);   uiwait(h),  return,
end
if ( numel(varargin) > 0 ),    hObject   = varargin{1};    end,
if ( numel(varargin) > 1 ),    eventdata = varargin{2};    end,
 
set(hObject,'State','off');
ButtonTag = get(hObject,'tag'); 
 
switch ButtonTag   %%#ok<*CTCH> 
 
    case 'action01'
        %%% specify callback action for button 1 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action02'
        %%% specify callback action for button 2 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action03'
        %%% specify callback action for button 3 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action04'
        %%% specify callback action for button 4 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action05'
        % specify callback action for button 5 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action06'
        % specify callback action for button 6 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    case 'action07'
        % specify callback action for button 7 here
        h = msgbox({' ';'specify callback action in';' ';['<auxiliary/',mfilename,'.m>']});
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
 
    otherwise
        h = warndlg('CALLBACK ACTION IS UNDEFINED');
        h = f_resize_dlgbox(h,14);   uiwait(h),  return,
end
end     % end of the 'maler_userbutt_callbacks' function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function h = f_resize_dlgbox(h,fsize)
%%% resize dialog boxes
    fsize = max(6,min(fsize,20));         % limiter
    htext = findobj(h, 'Type', 'Text');   % find text control in dialog box
    set(htext,'FontSize',fsize);
    set(h,'Resize','on');           pos = get(h,'Position');
    try
        deffsize = get(0,'factoryUicontrolFontSize');
    catch   %#ok
        deffsize = 8;
    end
    set(h,'Position',[pos(1)-(pos(3)*(fsize/deffsize-1)/2), pos(2), pos(3)*fsize/deffsize, pos(4)]);
end    % end of the function <f_resize_dlgbox>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

