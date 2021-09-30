function [EleClip] = Eleclip(DEM, levels)
% Clips a DEM to a specified elevation
% Input: Original DEM
% Output: clipped DEM
% instructions: Set upper limit based on contour map
figure('units','normalized','outerposition',[0 0 1 1])
imageschs(DEM);
hold on;
contour(DEM, levels, 'color', 'k', 'showtext', 'on')
prompt          = {'Define Upper limit'};
dlgtitle        = 'Input';
dims            = [1 35];
definput        = {'100'};
Answer_range          = inputdlg(prompt,dlgtitle,dims,definput);
UpperLimit = str2num(Answer_range{1});
close
%% Eleclip
L = DEM > UpperLimit; % mask variable
EleClip = clip(DEM, L) % clip function
% Asks if you need to clip the map area to a smaller region
promptMessage = sprintf('Create area mask?') % option to clip to a smaller area
titleBarCaption = 'Mask?'
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes')
if strcmpi(button,'Yes')==1;
    MASK = createmask(EleClip, 'usehillshade')
    % Clip to mask
    EleClip = clip(EleClip, MASK);
    % Crop to area
    EleClip = crop(EleClip, MASK);
    % Display
end

imageschs(EleClip)
end
