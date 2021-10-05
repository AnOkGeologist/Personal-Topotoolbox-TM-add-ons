% THIS VERSION OF THE CODE ADDS POINTS IF DIFFERENT SCANS HAVE DIFFERENT x
% (length) DIMENSIONS, WITHOUT DISTORTIONS OF THE SCANS. The points added
% are calculated from the scan, so ARE NOT REAL!
% The csv2DEM conversion portion of this was written by Dr. Riccardo
% Reitano 2021 to convert DEMs exported out of Paraview in csv format to DEMs readable by Topotoolbox (Schwanghart and Strak). 
clear; close all; clc
tic
selpath = uigetdir(path);
cd(selpath)  % We're moving into the folder where the files are
ls
%% Load DEM
prompt          = {'File prefix?', 'Digits from beg. of file?', 'Digits from end of file?'};
dlgtitle        = 'Input';
dims            = [1 35];
definput        = {'Test', '4', '4'};
Answer          = inputdlg(prompt,dlgtitle,dims,definput);

file_name = Answer{1};  % Self-explained
a = length(file_name);
b = str2num(Answer{2});  % Position of the digits indentifing the number of the files (e.g. Ero5_07.vtk, we are looking for "07") (from)
c = str2num(Answer{3});  % Position of the digits indentifing the number of the files (to)
cutx = 0; % How much cut from the laser borders (x), mm. MUST BE POSITIVE!
cuty = 0; % How much cut from the laser borders (y), mm. MUST BE POSITIVE!
pathname = pwd;  % Where we are

%% List directory files (including itself ...)
f_num=0;
d=dir(pathname); % a structure. All the elements in the chosen directory
pre_num=numel(d); % number of entries in directory
for fli=1:pre_num
       % only files
       if(d(fli).isdir==0)
           filn=d(fli).name;   % name field of listing structure
           fstrn=length(filn);   % Number of digits of the files name
           pref=filn(1:a);  % "prefix", just checking whether starting with the input "a"
           
           if(strcmp(pref,file_name))   % If the file name match with my input, then I'm loading it
               f_num=f_num+1;
               figfile{f_num}=[pathname,filesep,filn];  % Adding the file to the path in a new cell
               figstep(f_num)= str2num(filn(fstrn-b:fstrn-c));
           end
       end
end

%% I'm putting every scan in a cell of a cell-array
numfiles = f_num;
Q = cell(1,numfiles);  
for i = 1:numel(figfile)
    A = importdata(figfile{i});
    B = A.data;
%     B = B(~isnan(B));  
    C = reshape(B,numel(B)/3,3);  
    Q{i} = C;
end

%% Processing
for i = 1:numel(Q)
    idx = any((Q{i}(:,1) >= (max(Q{i}(:,1))-cutx)) | (Q{i}(:,2) <= (min(Q{i}(:,2))+cuty)) | (Q{i}(:,2) >= (max(Q{i}(:,2))-cuty)),2);
    Q{i}(idx,:) = [];
    resx = max(Q{i}(:,1))-min(Q{i}(:,1)); % resolution for linspace for griddata
    resy = max(Q{i}(:,2))-min(Q{i}(:,2));
    x = (Q{i}(:,1)+abs(min(Q{i}(:,1))));%*1e2;
    y = (Q{i}(:,2)+abs(min(Q{i}(:,2))));%*1e2;
    z = (Q{i}(:,3)+abs(min(Q{i}(:,3))));%*1e2;
    xq = linspace(min(x),max(x),resx);  % I'm choosing this min and max with resx (or resy),
    yq = linspace(min(y),max(y),resy);  % so that spacing between points should be 1 mm.
    [Xq,Yq] = meshgrid(xq,yq);   % I'm actually creating the grid
    % Now I need to interpolate the z-values coming from the scan, into the
    % new regular grid
    Z = griddata(x,y,z,Xq,Yq);
    Z = (Z.*1e2);
    Z = round(Z);
    Z(isnan(Z)) = -9999;
    Z = fliplr(Z);
    cell_size = round((Xq(1,end)-Xq(1,end-1))*1e2);

%% Rotate
    % Set up rotation
    crnr = -700000.00000000 %set corner
    xr = crnr:cell_size:crnr+length(xq)*cell_size-1; % Xs for GRIDobj
    yr = crnr:cell_size:crnr+length(yq)*cell_size-1; % Ys for GRIDobj
    DEMr = GRIDobj(xr,yr, Z) % create GRIDobj
    imageschs(DEMr) % show image
    i
    promptMessage = sprintf('Continue with rotation')
    titleBarCaption = 'Rotation'
    button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes')
    
    if strcmpi(button,'Yes')==1;
        
    display('select left and right ends or "rotation line" e.g. along strike')
    [Xr Yr] = ginput(2) % select strike
    % show selected points
    imageschs(DEMr); hold on;
    plot(Xr(1),Yr(1),'dk', 'MarkerFaceColor', 'r');
    plot(Xr(2),Yr(2),'dk', 'MarkerFaceColor', 'c')
    display('Press any key to continue and show rotation line')
    pause()
    close
    % fit selected points
    p = polyfit(Xr, Yr, 1)
    slope = p(1);
    f = polyval(p,Xr);
    imageschs(DEMr); hold on; plot(Xr(1),Yr(1),'dk', 'MarkerFaceColor', 'r');
    plot(Xr(2),Yr(2),'dk', 'MarkerFaceColor', 'c');
    plot(Xr,f, '-r')
    % Calculate angle from horizontal
    Angle = -atand(slope);
    display('Press any key to continue and rotate data')
    pause()
    close
    
    % Rotate data
    
    XY = [Xq(:) Yq(:)];                                     % Create Matrix Of Vectors
    R=[cosd(Angle) -sind(Angle); sind(Angle) cosd(Angle)]; %CREATE THE MATRIX
    rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
    
    % Reshape to query points
    Xqr = reshape(rotXY(:,1), size(Xq,1), []);
    Yqr = reshape(rotXY(:,2), size(Yq,1), []);
    Z = griddata(x,y,z, Xqr, Yqr);
    Z = (Z.*1e2);
    Z = round(Z);
    Z(isnan(Z)) = -9999;
    Z = fliplr(Z);
    
    crnr = -700000.00000000 %set corner
    xr = crnr:cell_size:crnr+length(xq)*cell_size-1; % Xs for GRIDobj
    yr = crnr:cell_size:crnr+length(yq)*cell_size-1; % Ys for GRIDobj
    DEMr = GRIDobj(xr,yr, Z) % create GRIDobj
    imageschs(DEMr) % show image 
    display('Press any key to continue')
    pause()
    
    end
    
    %% Save
    

    name = sprintf('DEM_%02i.txt',i);
    fid2 = fopen(name,'wt');
    header_nc = sprintf('ncols         %d',size(Z,2));
    header_nr = sprintf('nrows         %d',size(Z,1));
    header_xll = sprintf('xllcorner     -700000.00000000');
    header_yll = sprintf('yllcorner     -700000.00000000');
    cell_size = round((Xq(1,end)-Xq(1,end-1))*1e2);
    header_cs = sprintf('cellsize      %i.0000000000',cell_size);
    header_NaN = sprintf('NODATA_value  -9999');
    fprintf(fid2,[header_nc '\n' header_nr '\n' header_xll '\n' header_yll '\n' ...
        header_cs '\n' header_NaN '\n']);
    for j = 1:size(Z,1)
        fprintf(fid2,'%.1d ',Z(j,:));
        fprintf(fid2,'\n');
    end
    fclose(fid2);
    
end
toc

