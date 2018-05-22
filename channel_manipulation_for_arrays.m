

%% Define the layout of the electrodes on the arrays/probes
% as provided from the manufacturers

layout_PFC =  [NaN   2   1   3   4   6   8  10 NaN NaN;
                65  66  33  34   7  09  11  12  16  18;
                67  68  35  36   5  17  13  23  20  22;
                69  70  37  38  48  15  19  25  27  24;
                71  72  39  40  42  50  54  21  29  26;
                73  74  41  43  44  46  52  62  31  28;
                75  76  45  47  51  56  58  60  64  30;
                77  78  82  49  53  55  57  59  61  32;
                79  80  84  86  87  89  91  94  63  95;
               NaN  81  83  85  88  90  92  93  96  14]'; % 10x10 array : 96 Channels


layout_TEO =  [NaN   2   1   3   4   6   8  10  14  18;
                65  66  33  34   7   9  11  12  16 NaN;
                67  68  35  36   5  17  13  23  20 NaN;
                69  70  37  38  48  15  19  25  27  24;
                71  72  39  40  42  50  54  21  29  26;
                73  74  41  43  44  46  52  62  31  28;
                75  76  45  47  51  56  58  60  64  30;
                77  78  82  49  53  55  57  59  61  32;
                79  80  84  86  87  89  91  94  63  95;
               NaN  81  83  85  88  90  92  93  96  22]'; % 10x10 array : 96 Channels


layout_hippocampus = [[1:16]' [17:32]']; % 2x16 probe : 32 channels



% Insert the distance of the electrodes on the arrays
distance_TEO = 400*10^(-6); % 400um
distance_PFC = 400*10^(-6);
distance_PFC


%% The arrays are 2d objects. In order to place on them on the cortical surface
%  a plane needs to be estimated that is parallel to the cortex at the
%  impantation site. The easiest way to estimate the plane equation, is for
%  users to manually select 3 points near the area where the implant is
%  located, and then insert the values below (right click on an open
%  cortical surface -> Get Coordinates and click on three points near the
%  implant. Insert the MRI coordinates, in mm)

TEO_3_points         = [156.8 119.6 112.1;  % 1st point [X Y Z]
                        156.9 118.8 111.4;  % 2nd point [X Y Z]
                        156.6 117.3 111.9]; % 3rd point [X Y Z]
            
PFC_3_points         = [142.6 117.8 141.7;  % 1st point [X Y Z]
                        141.4 117.6 142.8;  % 2nd point [X Y Z]
                        141.4 116.0 142.0]; % 3rd point [X Y Z]
            
hippocampus_3_points = [127.7 123.5 115.4;  % 1st point [X Y Z]
                        127.8 122.5 113.9;  % 2nd point [X Y Z]
                        127.6 122.5 115.8]; % 3rd point [X Y Z]
                        
                   
%% Estimate the plane equations for each array/probe, based on the 3-points selected

% TEO
TEO_p12 = TEO_3_points(2,:)-TEO_3_points(1,:); % Vector from point 1 to point 2
TEO_p13 = TEO_3_points(3,:)-TEO_3_points(1,:);

normal_Vector_to_plane = cross(TEO_p12,TEO_p13); % The cross product wiil be normal to the plane


equation = normal_Vector_to_plane(1)*(x-TEO_3_points(1,1))+
           normal_Vector_to_plane(2)*(y-TEO_3_points(1,2))+
           normal_Vector_to_plane(3)*(z-TEO_3_points(1,3))  = 0

z = (normal_Vector_to_plane(1)*(TEO_3_points(1,1)-x)+...
     normal_Vector_to_plane(2)*(TEO_3_points(1,2)-y))/normal_Vector_to_plane(3) + TEO_3_points(1,3)












%% Load the channelMat from the recording (You can just right click on the Channel File -> File -> export to Matlab)
channels = load('E:/brainstorm_db/Playground/data/Monkey_No_Headstage/@rawBarcode_f096/channel.mat');
channels2 = channels;

for i = 1:224
    channels2.Channel(i).Type = 'SEEG';

    if i<97
        channels2.Channel(i).Group = 'PFC';
%         channels2.Channel(i).Type = 'ECOG';

    elseif i<193
        channels2.Channel(i).Group = 'TEO';
%         channels2.Channel(i).Type = 'ECOG';

    else
        channels2.Channel(i).Group = 'Hippocampus';
%         channels2.Channel(i).Type = 'SEEG';
    end
end



%PFC
for iElectrode = 1:96
    
    starting_point_PFC = [0.0372;-0.0112;0.0283];
    
    [a,b] = find(iElectrode == layout_PFC); 
    if ~isempty(a)
        
        % This is the equation for the surface where the array gets implanted
        z = -0.581876* ((a-1)*400*10^(-6)) + 0.575517* ((b-1)*400*10^(-6));%+ 0.0562652;
        
        channels2.Channel(iElectrode).Loc = starting_point_PFC + [(a-1)*400*10^(-6);(b-1)*400*10^(-6);z];
    else
        channels2.Channel(iElectrode).Loc = [0;0;0];
    end
end

%TEO
for iElectrode = 1:96
    
    starting_point_TEO = [0.0085;-0.026;0.002];
    
    [a,b] = find(iElectrode == layout_TEO); 
    
    a = a;
    b = b/2.1; % This is a hack to display it well. For some reason it gets fucked up from the equation.
    if ~isempty(a)
         % This is the equation for the surface where the array gets implanted
%          z = 0.0147059 x - 1.79412 y - 0.0447721
        z = 0.0147059* ((a-1)*400*10^(-6)) - 1.79412* ((b-1)*400*10^(-6)) -0.0447721*400*10^(-6);
        
        channels2.Channel(96+iElectrode).Loc =  starting_point_TEO+ [(a-1)*400*10^(-6);z;(b-1)*400*10^(-6)];
    else
        channels2.Channel(96+iElectrode).Loc = [0;0;0];
    end
end

%Hippocampus
for iElectrode = 1:32
    
    starting_point_Hippocampus = [0.0168;-0.0014; 0.0031];

    
    [a,b] = find(iElectrode == layout_hippocampus); 
    if ~isempty(a)
        
        z = 0;
        
        channels2.Channel(192+iElectrode).Loc = starting_point_Hippocampus + [z;(b-1)*150*10^(-6);(a-1)*150*10^(-6)];
    else
        channels2.Channel(192+iElectrode).Loc = [0;0;0];
    end
end


% length of hippocampal linear array: 4.65mm






figure(1);
hold on
for iElectrode = 1:96
    plot3(channels2.Channel(iElectrode).Loc(1), channels2.Channel(iElectrode).Loc(2), channels2.Channel(iElectrode).Loc(3),'.')
end


figure(2);
hold on
for iElectrode = 97:192
    plot3(channels2.Channel(iElectrode).Loc(1), channels2.Channel(iElectrode).Loc(2), channels2.Channel(iElectrode).Loc(3),'.')
end


figure(3);
hold on
for iElectrode = 193:224
    plot3(channels2.Channel(iElectrode).Loc(1), channels2.Channel(iElectrode).Loc(2), channels2.Channel(iElectrode).Loc(3),'.')
end























