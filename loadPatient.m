%{
Calibrates the RDFDM forward evaluation in 3D for GBM growth and treatment response
dN/dt = D*Delta(N) + kp*N*(1-N/theta)
kp = global or local depending on model_flags(2)
    *note - local kp requires parellel processing or becomes a lengthy calibration

Inputs:
    - loc: location for MRI and SPECT data
    - reduce: if on (reduce = 1) the simulation domain is shrunk to contain only the measured tumor we are calibrating + padding
    - target_dim: dimensions for dx, dy, and dz
Outputs:
    - tumor: information needed to run simulations other than parameters

%}
function tumor = loadPatient(loc, reduce, target_dim)
    adc_w = 2.5e3;
    theta = 818503; %Based on MRI voxel volume and size of GBM cells
    per = 0.10; %Simulation threshold, 10% carrying capacity
    
    mri_cnt = 0; mri_str = {}; t = [];
    tumor_check = 0; brain_check = 0;
    files = dir(loc);
    for i = 1:numel(files)
        if(contains(files(i).name, 'MRI')) %Check for MRI images
            %Load into temp file
            temp = load([loc, files(i).name]);
            mri_str{end+1} = ['day',num2str(temp.day)];
            
            t = [t, temp.day];
            
            eval(['day',num2str(temp.day),'_MRI = temp;']);
            
            mri_cnt = mri_cnt+1;
        elseif(contains(files(i).name, 'bulkTumor'))
            %Load masks, times should match days in mri_str
            TumorMasks = load([loc,files(i).name]);
            
            tumor_check = 1;
        elseif(contains(files(i).name, 'BrainMask'))
            %Load masks
            BrainMask = load([loc,files(i).name]);
            
            brain_check = 1;
        end
    end
    
    if(mri_cnt==0 || tumor_check==0 || brain_check==0)
        error('Required data not found in pathway, check data and try again');
    end
    
    %Cell Count Conversion
    %Pull out ADC maps and get minimum value in tumor
    adc_min = adc_w; %starts at max before search
    for i = 1:mri_cnt
        str = mri_str{i};
        
        %Get tumor mask
        eval([str,'_mask = logical(TumorMasks.',str,'_tumorMask);']);
        
        %Get ADC and minimum tumor ADC
        eval([str,'_adc = ',str,'_MRI.ADC;']);
        eval(['adc_min = min([min(',str,'_adc(find(',str,'_mask))) , adc_min], [], ''all'');']);
    end
    
    %Get brain slice and size
    brain = BrainMask.brain;
    [sy,sx,sz] = size(brain);
    
    %Set measured cell counts at each TP, N
    N = zeros(sy,sx,sz,mri_cnt);
    for i = 1:mri_cnt
        str = mri_str{i};
        eval([str,'_mask(~find(brain))=0;']); %zero mask outside the brain
        eval(['temp = ',str,'_mask.*(theta*(adc_w - ',str,'_adc)./(adc_w - adc_min));']); %convert to cells
        N(:,:,:,i) = temp;
    end
    N(N<0) = 0;
    
    dx = target_dim(1); dy = target_dim(2); dz = target_dim(3);
    %Interpolate measured data to desired spacing
    [X,Y,Z] = meshgrid(1:sx,1:sy,1:sz);
    [X_new,Y_new,Z_new] = meshgrid(linspace(1,sx,sx/dx), linspace(1,sy,sy/dy), linspace(1,sz,sz/dz));
    
    N_new = zeros([size(X_new),mri_cnt]);
    for i = 1:mri_cnt
        N_new(:,:,:,i) = dx*dy*dz*interp3(X,Y,Z,N(:,:,:,i),X_new,Y_new,Z_new); %Assumes initial image has resolution 1x1x1 mm
    end
    brain = interp3(X,Y,Z,brain,X_new,Y_new,Z_new); brain(brain>0) = 1;
    N = N_new;
    
    [sy,sx,sz,~] = size(N);
    theta = theta*dx*dy*dz;

    %Update thresh after resizing
    thresh = theta*per;
    
    %Apply image smoothing, leave theta at max of unfiltered image
    temp_idx = [];
    for i = 1:mri_cnt
            N(:,:,:,i) = imgaussfilt3(N(:,:,:,i));
            temp_idx = union(temp_idx, find(N(:,:,:,i)));
    end
    N(N<thresh) = 0;
    
    %Ensure boundary includes whole tumor region
    brain(temp_idx) = 1;
    
    if(reduce == 1)
        padding = 5;
        %Find edges
        temp = sum(N,4); temp(temp~=0) = 1;
        temp_xy = sum(temp,3);
        edges = bwboundaries(temp_xy);
        y_low = min(edges{1}(:,1)) - padding;
        if(y_low<1)
            y_low = 1;
        end
        y_up  = max(edges{1}(:,1)) + padding;
        if(y_up>sy)
            y_up = sy;
        end
        
        x_low = min(edges{1}(:,2)) - padding;
        if(x_low<1)
            x_low = 1;
        end
        x_up = max(edges{1}(:,2)) + padding;
        if(x_up>sx)
            x_up = sx;
        end
        temp_xz = squeeze(sum(temp,1));
        edges = bwboundaries(temp_xz);
        z_low = min(edges{1}(:,2)) - padding;
        if(z_low<1)
            z_low = 1;
        end
        z_up = max(edges{1}(:,2)) + padding;
        if(z_up>sz)
            z_up = sz;
        end
        
        %Reduce cell map
        N = N(y_low:y_up, x_low:x_up, z_low:z_up, :);
        
        %Reduce brain
        brain = brain(y_low:y_up, x_low:x_up, z_low:z_up);
        
    end
    clear day0_mask day28_mask day56_mask day112_mask 
    clear temp temp_xy temp_xz edges
    clear x_low x_up y_low y_up z_low z_up
    
    [sy,sx,sz,~] = size(N);
    
    %Build boundaries
    bcs = BuildBoundaries(brain);
    
    tumor.N = N;
    tumor.bcs = bcs;
    tumor.dx = dx;
    tumor.dy = dy;
    tumor.dz = dz;
    tumor.t = t;
    tumor.theta = theta;
end

