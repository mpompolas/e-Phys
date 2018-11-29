%% Select the 



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


layout_hippocampus = [1:32]; % 1:32 probe : 32 channels



% Insert the distance of the electrodes on the arrays (this assumes all electrodes have the same distance on the arrays)
distance_TEO         = 400*10^(-6); % 400um
distance_PFC         = 400*10^(-6);
distance_hippocampus = 150*10^(-6);


%% The arrays are 2d objects. In order to place on them on the cortical surface
%  a plane needs to be estimated that is parallel to the cortex at the
%  impantation site. The easiest way to estimate the plane equation, is for
%  users to manually select 3 points near the area where the implant is
%  located, and then insert the values below (right click on an open
%  cortical surface -> Get Coordinates and click on three points near the
%  implant. Insert the SCS coordinates in mm. It is divided by 1000 to convert to m)

TEO_3_points         = [ 6.5 -29.0 14.3;        % 1st point [X Y Z]
                         7.4 -28.8 15.1;        % 2nd point [X Y Z]
                         7.7 -28.9 14.2]./1000; % 3rd point [X Y Z]
            
PFC_3_points         = [39.5 -10.6 27.4;        % 1st point [X Y Z]
                        39.0  -9.5 28.2;        % 2nd point [X Y Z]
                        39.2 -11.8 26.6]./1000; % 3rd point [X Y Z]
            
hippocampus_3_points = [15.1  -0.9  6.5;        % 1st point [X Y Z]
                        17.2  -0.5  6.2;        % 2nd point [X Y Z]
                        15.7  -0.7  5.7]./1000; % 3rd point [X Y Z]
                        
                    
%% Load the channelMat from the recording (You can just right click on the Channel File -> File -> export to Matlab)
channels = load('E:/brainstorm_db/Playground/data/Monkey_No_Headstage/@rawBarcode_f096/channel.mat');
channels2 = channels; % Everything will be saved on a channels2 struct

for i = 1:224
    channels2.Channel(i).Type = 'SEEG';

    if i<97
        channels2.Channel(i).Group = 'TEO';
%         channels2.Channel(i).Type = 'ECOG';

    elseif i<193
        channels2.Channel(i).Group = 'PFC';
%         channels2.Channel(i).Type = 'ECOG';

    else
        channels2.Channel(i).Group = 'Hippocampus';
%         channels2.Channel(i).Type = 'SEEG';
    end
end






%% Do the computation


% Plane equation: Ax+By+Cz+D =0;
%                  z = -(Ax+By-D)/C

                 
%PFC
% Compute the equation of the surface
PFC_p12 = PFC_3_points(2,:)-PFC_3_points(1,:); % Vector from point 1 to point 2
PFC_p13 = PFC_3_points(3,:)-PFC_3_points(1,:);
normal_Vector_to_plane = cross(PFC_p12,PFC_p13); % The cross product wiil be normal to the plane
D = -(normal_Vector_to_plane(1)*PFC_3_points(1,1)+normal_Vector_to_plane(2)*PFC_3_points(1,2)+normal_Vector_to_plane(3)*PFC_3_points(1,3));

% Get two vectors on the  plane
% v1 = [B/A,-1,0];
% v2 = [C/A,0,-1];
V = [normal_Vector_to_plane(2)/normal_Vector_to_plane(1) normal_Vector_to_plane(3)/normal_Vector_to_plane(1);-1 0;0 -1];

% Perform Gram-Schmidt Process (Get orthonormal vectors from two vectors on
% a plane) - Columns of U are the orthonormal basis
n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:k
  U(:,i) = V(:,i);
  for j = 1:i-1
    U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
  end
  U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end

for iElectrode = 1:96
    starting_point_PFC = PFC_3_points(1,:)';
    [a,b] = find(iElectrode == layout_PFC); 
   
    if ~isempty(a)        
%         z = -(normal_Vector_to_plane(1)*(PFC_3_points(1,1)+(a-1)*distance_PFC)+...
%               normal_Vector_to_plane(2)*(PFC_3_points(1,2)+(b-1)*distance_PFC)+...
%               D                                                                   )./normal_Vector_to_plane(3);
%           
%         channels2.Channel(iElectrode).Loc = starting_point_PFC + [normal_Vector_to_plane(1)/normal_Vector_to_plane(3)*(a-1)*distance_PFC;normal_Vector_to_plane(2)/normal_Vector_to_plane(3)*(b-1)*distance_PFC;z];
        channels2.Channel(iElectrode).Loc = starting_point_PFC + (a-1)*distance_PFC*U(:,1) + (b-1)*distance_PFC*U(:,2);
    else
        channels2.Channel(iElectrode).Loc = [0;0;0];
    end
end


%TEO
% Compute the equation of the surface
TEO_p12 = TEO_3_points(2,:)-TEO_3_points(1,:); % Vector from point 1 to point 2
TEO_p13 = TEO_3_points(3,:)-TEO_3_points(1,:);
normal_Vector_to_plane = cross(TEO_p12,TEO_p13); % The cross product wiil be normal to the plane
D = -(normal_Vector_to_plane(1)*TEO_3_points(1,1)+normal_Vector_to_plane(2)*TEO_3_points(1,2)+normal_Vector_to_plane(3)*TEO_3_points(1,3));

% Get two vectors on the  plane
% v1 = [B/A,-1,0];
% v2 = [C/A,0,-1];
V = [normal_Vector_to_plane(2)/normal_Vector_to_plane(1) normal_Vector_to_plane(3)/normal_Vector_to_plane(1);-1 0;0 -1];

% Perform Gram-Schmidt Process (Get orthonormal vectors from two vectors on
% a plane) - Columns of U are the orthonormal basis
n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:k
  U(:,i) = V(:,i);
  for j = 1:i-1
    U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
  end
  U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end

for iElectrode = 1:96
    starting_point_TEO = TEO_3_points(1,:)';
    [a,b] = find(iElectrode == layout_TEO); 
    
    if ~isempty(a)
        
%         z = -(normal_Vector_to_plane(1)*(TEO_3_points(1,1)+(a-1)*distance_TEO)+...
%               normal_Vector_to_plane(2)*(TEO_3_points(1,2)+(b-1)*distance_TEO)+...
%               D                                                                   )/normal_Vector_to_plane(3);
          
%         channels2.Channel(96+iElectrode).Loc = starting_point_TEO + [normal_Vector_to_plane(1)*(a-1)*distance_TEO;normal_Vector_to_plane(2)*(b-1)*distance_TEO;z];
          channels2.Channel(96+iElectrode).Loc = starting_point_TEO + (a-1)*distance_TEO*U(:,1) + (b-1)*distance_TEO*U(:,2);

    else
        channels2.Channel(96+iElectrode).Loc = [0;0;0];
    end
end



%Hippocampus
% Compute the equation of the surface
hippocampus_p12 = hippocampus_3_points(2,:)-hippocampus_3_points(1,:); % Vector from point 1 to point 2
hippocampus_p13 = hippocampus_3_points(3,:)-hippocampus_3_points(1,:);
normal_Vector_to_plane = cross(hippocampus_p12,hippocampus_p13); % The cross product wiil be normal to the plane
D = -(normal_Vector_to_plane(1)*hippocampus_3_points(1,1)+normal_Vector_to_plane(2)*hippocampus_3_points(1,2)+normal_Vector_to_plane(3)*hippocampus_3_points(1,3));

% Get two vectors on the  plane
% v1 = [B/A,-1,0];
% v2 = [C/A,0,-1];
V = [normal_Vector_to_plane(2)/normal_Vector_to_plane(1) normal_Vector_to_plane(3)/normal_Vector_to_plane(1);-1 0;0 -1];

% Perform Gram-Schmidt Process (Get orthonormal vectors from two vectors on
% a plane) - Columns of U are the orthonormal basis
n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:k
  U(:,i) = V(:,i);
  for j = 1:i-1
    U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
  end
  U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end


for iElectrode = 1:32
    starting_point_hippocampus = hippocampus_3_points(1,:)';
    [a,b] = find(iElectrode == layout_hippocampus); 
    
    if ~isempty(a)
        
        z = -(normal_Vector_to_plane(1)*(hippocampus_3_points(1,1)+(a-1)*distance_hippocampus)+...
              normal_Vector_to_plane(2)*(hippocampus_3_points(1,2)+(b-1)*distance_hippocampus)+...
              D                                                                                )/normal_Vector_to_plane(3);
          
%         channels2.Channel(192+iElectrode).Loc = starting_point_hippocampus + [normal_Vector_to_plane(1)*(a-1)*distance_hippocampus;normal_Vector_to_plane(2)*(b-1)*distance_hippocampus;z];
          channels2.Channel(192+iElectrode).Loc = starting_point_hippocampus + (a-1)*distance_hippocampus*U(:,1) + (b-1)*distance_hippocampus*U(:,2);
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
hold off


figure(2);
hold on
for iElectrode = 97:192
    plot3(channels2.Channel(iElectrode).Loc(1), channels2.Channel(iElectrode).Loc(2), channels2.Channel(iElectrode).Loc(3),'.')
end
hold off


figure(3);
hold on
for iElectrode = 193:224
    plot3(channels2.Channel(iElectrode).Loc(1), channels2.Channel(iElectrode).Loc(2), channels2.Channel(iElectrode).Loc(3),'.')
end
hold off


clear channels




%% Check if I can get back the indices of the array to create rows and columns


Channels = channels2.Channel;


channelsMontage = zeros(length(Channels),1); % This holds the code of the montage each channel holds 
channelsCoords  = zeros(length(Channels),3); % THE 3D COORDINATES
Montages = unique({Channels.Group});

for iChannel = 1:length(Channels)
    for iMontage = 1:length(Montages)
        if strcmp(Channels(iChannel).Group, Montages{iMontage})
            channelsMontage(iChannel)    = iMontage;
            channelsCoords(iChannel,1:3) = Channels(iChannel).Loc;
        end
    end
end


for iMontage = 1:length(Montages)
    clear single_array_coords
    single_array_coords = channelsCoords(channelsMontage==iMontage,:);

    
    distances = zeros(size(single_array_coords,1),size(single_array_coords,1));
    for iChannel = 1:size(single_array_coords,1)
        for jChannel = 1:size(single_array_coords,1)
            distances(iChannel,jChannel) = sqrt((single_array_coords(iChannel,1)-single_array_coords(jChannel,1))^2 + ...
                                                (single_array_coords(iChannel,2)-single_array_coords(jChannel,2))^2 + ...
                                                (single_array_coords(iChannel,3)-single_array_coords(jChannel,3))^2 ...
                                               );
            if iChannel==jChannel || iChannel>jChannel
                distances(iChannel, jChannel) = NaN;
            end
        end
    end

    
    if isMontageLinear(iMontage)
        
        quantum_distance = min(min(distances));
    
        % Find the longest distance between eletctrodes, and assign the first of
        % them as the origin
        [temp,~] = max(distances');
        [temp,iStart] = max(temp);
        [temp,~] = max(distances);
        [temp,iEnd] = max(temp);
        iStart
        iEnd
        
        indices_X = round(distances(iStart,:)./quantum_distance)+1;
        indices_X(1) = 1;
        indices_Y = ones(1,length(indices_X));
        
    else
        quantum_distance = min(min(distances));

        two_perpendicular_vectors = zeros(2,3);
        two_perpendicular_vectors(1,:) = single_array_coords(2,:) - single_array_coords(1,:); % Assign a random first vector
        justStop = 0;
        iVector = 1;
        % Find two vectors that are perpendicular to each other
        for iChannel = 1:size(single_array_coords,1)
            for jChannel = 1:size(single_array_coords,1)
                other_vector  = single_array_coords(jChannel,:) - single_array_coords(iChannel,:);
                other_vector(abs(other_vector)<10^(-8))=0; % Get rid of precision errors
                % Run until you find a second 
                if sum(cross(two_perpendicular_vectors(1,:),other_vector'))~=0 && iChannel~=jChannel && (distances(iChannel,jChannel) - quantum_distance <10^(-10))
                    two_perpendicular_vectors(iVector,:) = single_array_coords(jChannel,:) - single_array_coords(iChannel,:);
                    iChannel
                    jChannel
                    iVector = iVector +1;
                    justStop = 1;
                    break
                end
            end
            if justStop && iVector == 3
                break
            end
        end
     
        
        % Find the longest distance between eletctrodes, and assign the first of
        % them as the origin
        [temp,~] = max(distances');
        [temp,iStart] = max(temp);
        [temp,~] = max(distances);
        [temp,iEnd] = max(temp);
        iStart
        iEnd
        
        
        % Create a grid based on the quantumCoords starting from iStart
        
        grid = zeros(size(single_array_coords,1)*2+1,size(single_array_coords,1)*2+1,3); % 96 channel array, I create 193x193x3 in case the start is not at the start
        NEW_2D_COORDS = zeros(size(single_array_coords,1),2);
        
        
        for xGrid = -size(single_array_coords,1):size(single_array_coords,1)
            for yGrid = -size(single_array_coords,1):size(single_array_coords,1)
                
                grid(size(single_array_coords,1)+xGrid +1,size(single_array_coords,1)+yGrid+1,:) = single_array_coords(iStart,:)+ ...
                                                                             xGrid * two_perpendicular_vectors(1,:) + yGrid * two_perpendicular_vectors(2,:);
                
                for iChannel = 1:size(single_array_coords,1)
                    point_from_grid_entry_distance = single_array_coords(iStart,:)+xGrid * two_perpendicular_vectors(1,:) + yGrid * two_perpendicular_vectors(2,:) - single_array_coords(iChannel,:);                                           
                    point_from_grid_entry_distance(abs(point_from_grid_entry_distance)<10^(-10))=0;
                    if sum(point_from_grid_entry_distance)==0
                        NEW_2D_COORDS(iChannel,:) = [xGrid yGrid];
                        break
                    end
                end
                                                                         
                                                                         
            end
        end
        
        shift = min(NEW_2D_COORDS); % Check if the origin was not at the edge
        if shift<1
            NEW_2D_COORDS = NEW_2D_COORDS + abs(min(NEW_2D_COORDS))+1;
        end
    end
end
















