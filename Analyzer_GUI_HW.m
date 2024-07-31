function varargout = Analyzer_GUI_HW(varargin) %Call Analyzer_GUI_HW
% Analyzer_GUI_HW MATLAB code for Analyzer_GUI_HW.fig

% Last Modified by GUIDE v2.5 14-Sep-2018 22:54:36
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Analyzer_GUI_HW_OpeningFcn, ...
    'gui_OutputFcn',  @Analyzer_GUI_HW_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT
% --- Executes just before Analyzer_GUI_HW is made visible.
function Analyzer_GUI_HW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Analyzer_GUI_HW (see VARARGIN)
% Choose default command line output for Analyzer_GUI_HW
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Parameter initialization



% UIWAIT makes Analyzer_GUI_HW wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = Analyzer_GUI_HW_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Pathbutton.
function Pathbutton_Callback(hObject, eventdata, handles)
warning('off','all');

path=uigetdir('D:\','Select directory to analyze');
set(handles.text4,'String',path);
guidata(hObject,handles);

function thres_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function thres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_I_peak_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function min_I_peak_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ini_f_Callback(hObject, eventdata, handles)

function ini_f_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fin_f_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function fin_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fin_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes when selected object is changed in LR_selec.
function LR_selec_SelectionChangedFcn(hObject, eventdata, handles)

%% --- Executes when selected object is changed in Visual_Ana.
function Visual_Ana_SelectionChangedFcn(hObject, eventdata, handles)


%% --- Executes on button press in Analbutton.
function Analbutton_Callback(hObject, eventdata, handles)
%% Set criteria
ch_selector=get(get(handles.LR_selec,'SelectedObject'),'Tag');
Visual_opt=get(get(handles.Visual_Ana,'SelectedObject'),'Tag');
color=get(handles.colorindicator,'Backgroundcolor');

set_thres=get(handles.thres,'String');
set_min_I_peak=get(handles.min_I_peak,'String');
set_ini_f=get(handles.ini_f,'String');
set_fin_f=get(handles.fin_f,'String');
mode=get(handles.check_inten,'Value');
mode_I=get(handles.Surf_Int,'Value');
mode_TL=get(handles.check_timelapse,'Value');
Img_min=str2num(get(handles.Img_min,'String'));
Img_max=str2num(get(handles.Img_max,'String'));

set_window=get(handles.TL_window,'String');
set_interval=get(handles.TL_interval,'String');
window_TL=str2num(set_window); %Number of averaging frames
interval_TL=str2num(set_interval); %Interval

if mode==1
    set_dist_cri=get(handles.dist_info,'String');
    dist_cri=str2num(set_dist_cri);
    info_cri=set_dist_cri;
else
    dist_cri=0;
    info_cri='off';
end

threshold=str2num(set_thres); % Image planation
min_I_peak=str2num(set_min_I_peak); % Cutoff value
ini_f=str2num(set_ini_f) ;% Frames to be analyzed
fin_f=str2num(set_fin_f);

%% Obtain Path info
path=get(handles.text4,'String');
dir_info=dir(path);
dir_info(~[dir_info.isdir]) = [];  %remove non-directories
tf = ismember({dir_info.name}, {'.', '..'});
dir_info(tf) = [];

%% query whether subfolder was already analyzed
fol_list=false(length(dir_info),1); % default, all false
for k=1:length(dir_info)
    sub_path=[path '\' dir_info(k).name];
    subdir_info=dir(sub_path);
    tf=ismember({subdir_info.name},{'result.dat'});
    if any(tf) || length(subdir_info)==2
        fol_list(k)=true;
    end
end

%% Tif file searching
dir_info(fol_list)=[];
tic
for j=1:length(dir_info)  % sub-folder searching
    sub_path=[path '\' dir_info(j).name];
   
    cd(sub_path);
    disp(sub_path)
    subdir_info=dir([sub_path '\' '*.tif']);
    nmol=zeros(1,length(subdir_info));
    Label_fname=cell(1,1);
   
    if mode_I==1
        Surf_Int=zeros(1,length(subdir_info));
    end
    
    for h=1:length(subdir_info) % file loading

        temp_Img=loadTif16([sub_path '\' subdir_info(h).name]);
        temp_Img=imrotate(temp_Img,-90); 
        %Crop image / Select frames / channel
        switch ch_selector
            case 'Top'
                Tif_raw(:,:)=uint16(mean(temp_Img(:,257:512,ini_f:fin_f),3));
            case 'Bottom'
                Tif_raw(:,:)=uint16(mean(temp_Img(:,1:256,ini_f:fin_f),3));
        end
        [num_mol, x_pos, y_pos]=peakfinder_HW(Tif_raw, threshold, min_I_peak, mode, dist_cri); % nmol=number of selected molecules / x_pos=x-position array & y_pos=y-position array for selected peak
        nmol(h)=num_mol;
        pos=[x_pos, y_pos]; %position information for selected peaks
        FigGenerator_HW(Tif_raw, subdir_info(h).name(1:end-4), pos,Img_min,Img_max,Visual_opt,color);
        %fname=[subdir_info(h).name(1:end-4) '.dat'];
        %save(fname,'pos','-ascii');
        
        if mode_I==1
            Surf_Int(h)=sum(sum(Tif_raw));
        end
        
       %% RT counting mode
        if mode_TL==1
            loop_max=fix((size(temp_Img,3)-window_TL)/interval_TL); % Time Lapse Iteration number
            nmol_TL_temp=zeros(loop_max,1);
            Interval_TL=cell(loop_max+1,1);
            
            Tif_raw_TL=zeros(512,256,size(temp_Img,3));
            switch ch_selector
                case 'Top'
                    Tif_raw_TL(:,:,:)=temp_Img(:,257:512,:);
                case 'Bottom'
                    Tif_raw_TL(:,:,:)=temp_Img(:,1:256,:);
            end
            for j=0:loop_max
                ith_TL=1+j*(interval_TL);
                jth_TL=j*(interval_TL)+window_TL;
                Avg_Tif_TL=zeros(512,256);
                Avg_Tif_TL(:,:)=mean(Tif_raw_TL(:,1:256,ith_TL:jth_TL),3);
                [nmol_TL_temp(j+1,1), x_pos_TL, y_pos_TL]=peakfinder_HW(Avg_Tif_TL, threshold, min_I_peak, mode, dist_cri); % nmol=number of selected molecules / x_pos=x-position array & y_pos=y-position array for selected peak
                Interval_TL{j+1,1}=[num2str(ith_TL) '-' num2str(jth_TL)];
            end
            nmol_TL(:,h)=nmol_TL_temp(:,1);
        end
        Label_fname{1,h}=subdir_info(h).name(1:end-4);
    end
    
    Label_fname{1,size(Label_fname,2)+1}=['Mean'];
    Label_fname{1,size(Label_fname,2)+1}=['Std'];

    nmol(size(nmol,2)+1)=mean(nmol(1:end)); % mean of num_mol in the same folder
    nmol(size(nmol,2)+1)=std(nmol(1:end-1)); % std of num_mol in the same folder
    result_table=table(Label_fname', nmol', 'VariableNames',{'Name' 'Value'});
    writetable(result_table,'result.dat','Delimiter','\t');    
    
    if mode_I==1     
        Surf_Int(size(Surf_Int,2)+1)=mean(Surf_Int(1:end)); % mean of num_mol in the same folder
        Surf_Int(size(Surf_Int,2)+1)=std(Surf_Int(1:end-1)); % std of num_mol in the same folder
        
        result_table_Int=table(Label_fname', Surf_Int', 'VariableNames',{'Name' 'Value'});
        writetable(result_table_Int,'result_Int.dat','Delimiter','\t');
        %save(['result_Int.dat'],'Surf_Int','-ascii');
    end
    
    if mode_TL==1
        nmol_TL(:,size(nmol_TL,2)+1)=mean(nmol_TL,2);
        nmol_TL(:,size(nmol_TL,2)+1)=std(nmol_TL(:,1:(size(nmol_TL,2)-1)),0,2);
        result_TL=array2table(nmol_TL, 'RowNames', Interval_TL);
        result_TL.Properties.VariableNames=Label_fname;
        writetable(result_TL,'result_TL.dat','Delimiter','\t','WriteRowNames',true);
        clear nmol_TL;
    end
        
    %% Creat Analysis information report
    Name_info={'Background';'Min_peak_criteria';'Analyzed frames';'Analyzed ch';'Overlap criteria'};
    Value_info={set_thres;set_min_I_peak;[set_ini_f ' to ' set_fin_f];ch_selector;info_cri};
    T_info=table(Value_info,'RowNames', Name_info);
    writetable(T_info,'Analysis_info.txt','Delimiter','\t','WriteRowNames',true);    
end
cd(path);
toc
disp('Analyzed')
% hObject    handle to Analbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Del_bt.
function Del_bt_Callback(hObject, eventdata, handles)
bt_result=questdlg('Are you sure?','Removing all result.dat in your folders','Yes','No','No');

if bt_result=='Yes'
    path=get(handles.text4,'String');
    dir_info=dir(path);
    dir_info(~[dir_info.isdir]) = [];  %remove non-directories
    tf = ismember( {dir_info.name}, {'.', '..'});
    dir_info(tf) = [];
    
    %% query all subfolders
    for k=1:length(dir_info)
        sub_path=[path '\' dir_info(k).name];
        sub_path1=[sub_path '\' 'result.dat'];
        sub_path2=[sub_path '\' 'result_Int.dat'];
        delete(sub_path1);
        delete(sub_path2);
    end
    
    disp('Clear')
end

% hObject    handle to Del_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Sum_rep.
function Sum_rep_Callback(hObject, eventdata, handles)
path=get(handles.text4,'String');
dir_info=dir(path);
dir_info(~[dir_info.isdir]) = [];  %remove non-directories
tf = ismember({dir_info.name},{'.', '..'});
dir_info(tf) = [];
mode_I=get(handles.Surf_Int,'Value');
%% query all subfolders
for k=1:length(dir_info)
    sub_path=[path '\' dir_info(k).name];
    subdir_info=dir(sub_path);
    tf=ismember({subdir_info.name},{'result.dat'});
    fol_name=sub_path(length(path)+2:end);
       
    if any(tf)
        %temp_dat=importdata([sub_path '\' 'result.dat']);
        temp_dat=readtable([sub_path '\' 'result.dat']);
        load_value=temp_dat.Value';
        temp_name{k,1}=fol_name;
        temp_stat(k,1:2)=num2cell(load_value(end-1:end));
        temp_raw(k,1:length(temp_dat.Value)-2)=num2cell(load_value(1:end-2)');
    else
        %temp_stat(k,1:2)={fol_name, 'result.dat file is not existed',''};
    end
    
    if mode_I==1
        temp_dat2=readtable([sub_path '\' 'result_Int.dat']);
        load_value2=temp_dat2.Value';
        temp_name2{k,1}=fol_name;
        temp_stat2(k,1:2)=num2cell(load_value2(end-1:end));
        temp_raw2(k,1:length(temp_dat2.Value)-2)=num2cell(load_value2(1:end-2)');
    end   
end

fname=get(handles.fname_sum,'String');
stat=[temp_name, temp_stat];
disp(temp_name');
T1=cell2table(stat,'VariableNames',{'folder', 'mean', 'std'});
T2=cell2table([temp_name'; temp_raw']);

writetable(T1,[path '\' fname '.dat'],'Delimiter','\t');
writetable(T2,[path '\' fname 'Raw.dat'],'Delimiter','\t');

if mode_I==1
    stat2=[temp_name2, temp_stat2];
    disp(temp_name2');
    T3=cell2table(stat2,'VariableNames',{'folder', 'mean', 'std'});
    T4=cell2table([temp_name2'; temp_raw2']);
    
    writetable(T3,[path '\' fname '_Int.dat'],'Delimiter','\t');
    writetable(T4,[path '\' fname '_Int_Raw.dat'],'Delimiter','\t');
end



disp('Summarized')

% hObject    handle to Sum_rep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on Top and none of its controls.
% --- Executes on button press in check_inten.
function check_inten_Callback(hObject, eventdata, handles)
status=get(handles.check_inten,'Value');
if status == 1
    %Allow input for editMaxResults
    set(handles.dist_info, 'enable', 'on')
else
    %Block input for editMaxresults
    set(handles.dist_info, 'enable', 'off')
end


% hObject    handle to check_inten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_inten

function dist_info_Callback(hObject, eventdata, handles)
% hObject    handle to dist_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_info as text
%        str2double(get(hObject,'String')) returns contents of dist_info as a double

% --- Executes during object creation, after setting all properties.
function dist_info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in colorindicator.
function colorindicator_Callback(hObject, eventdata, handles)
color=uisetcolor([0 1 0]);
set(handles.colorindicator,'BackgroundColor',color);
% hObject    handle to colorindicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Surf_Int.
function Surf_Int_Callback(hObject, eventdata, handles)
status=get(handles.Surf_Int,'Value');
% hObject    handle to Surf_Int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Surf_Int


% --- Executes on button press in check_timelapse.
function check_timelapse_Callback(hObject, eventdata, handles)
status=get(handles.check_timelapse,'Value');
if status == 1
    %Allow input for editMaxResults
    set(handles.TL_window, 'enable', 'on')
    set(handles.TL_interval, 'enable', 'on')
else
    %Block input for editMaxresults
    set(handles.TL_window, 'enable', 'off')
    set(handles.TL_interval, 'enable', 'off')
end

% hObject    handle to check_timelapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_timelapse



function TL_window_Callback(hObject, eventdata, handles)
% hObject    handle to TL_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TL_window as text
%        str2double(get(hObject,'String')) returns contents of TL_window as a double


% --- Executes during object creation, after setting all properties.
function TL_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TL_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TL_interval_Callback(hObject, eventdata, handles)
% hObject    handle to TL_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TL_interval as text
%        str2double(get(hObject,'String')) returns contents of TL_interval as a double


% --- Executes during object creation, after setting all properties.
function TL_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TL_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on Bottom and none of its controls.
function Bottom_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Bottom (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on Top and none of its controls.
function Top_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Top (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function Img_min_Callback(hObject, eventdata, handles)
% hObject    handle to Img_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Img_min as text
%        str2double(get(hObject,'String')) returns contents of Img_min as a double


% --- Executes during object creation, after setting all properties.
function Img_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Img_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Img_max_Callback(hObject, eventdata, handles)
% hObject    handle to Img_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Img_max as text
%        str2double(get(hObject,'String')) returns contents of Img_max as a double


% --- Executes during object creation, after setting all properties.
function Img_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Img_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fname_sum_Callback(hObject, eventdata, handles)
% hObject    handle to fname_sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fname_sum as text
%        str2double(get(hObject,'String')) returns contents of fname_sum as a double


% --- Executes during object creation, after setting all properties.
function fname_sum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fname_sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
