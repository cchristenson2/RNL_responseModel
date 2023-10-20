%{
Calibrates the RDFDM forward evaluation in 3D for GBM growth and treatment response
dN/dt = D*Delta(N) + kp*N*(1-N/theta)
kp = global or local depending on model_flags(2)
    *note - local kp requires parellel processing or becomes a lengthy calibration

Author: Chase Christenson - 10/18/2023

Reference:
Hormuth II DA, Eldridge SL, Weis JA, et al. 
Mechanically coupled reaction-diffusion model to predict Glioma growth: methodological details. 
Methods Mol Biol. 2018;1711:225–241.


Inputs:
    - loc: location for MRI and SPECT data
    - iterations: maximum number of fitting iterations
    - model_flags: model specifiers [cal_flag, kp_flag]
        % cal_flag: 1 = calibrate to all data available; 2 = calibrate to all-1 and predict for extra
        % kp_flag:  1 = global; 2 = local(non-zero inside cal TP bounds + zero outside (adaptive proliferation));
    - reduce: if on (reduce = 1) the simulation domain is shrunk to contain only the measured tumor we are calibrating + padding
Outputs:
    - Params_cal: calibrated diffusivity, proliferation
    - Outputs_cal: Calibrated and simulated cell maps

%}

function [Params_cal, Outputs_cal] = Calibrate_RDFDM_3D_example(loc, max_it, model_flags, reduce)
    %LM Parameters
    alpha = 1e0;
    alpha_min = 1e-10;
    alpha_max = 1e10;
    pass = 7;
    fail = 10;
    
    e_change    = 1e-7;   %minimum change in error for stopping, convergence reached
    e_tol       = 1e-5;   %tolerance value, error minimized in residuals
    
    normalize = 0; %Normalize scale of parameters for calibration (1 = on)
    %%% Can help with convergence for local kp %%%
    
    %Jacobian perturbations
    delta = 1.001;
    j_freq = 1; j_change = j_freq; %How frequently to update jacobian

    %Simulation Parameters
    dx = 2; %mm
    dy = 2; %mm
    dz = 2; %mm
    dt = 0.25; %days
    
    %Cell/Voxel parameters
    theta = 818503; %Based on MRI voxel volume and size of GBM cells
    per = 0.10; %Simulation threshold, 10% carrying capacity
    adc_w = 2.5*(10^3); %ADC of water

    %% Functional Block 1 - Load data and apply boundaries for model
    disp('Loading data...');
    
    %Searches file in the pathway (searches since different patients will have scans on different days)
    %All data is registered prior to loading
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
    disp('Cell count conversion...');
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

    clear day0_MRI day28_MRI day56_MRI day112_MRI
    clear day0_adc day28_adc day56_adc day112_adc
    
    clear BrainMask
    clear TumorMasks
    clear N_new day0_mask day28_mask day56_mask day112_mask
    
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
    %% Function Block 2 - Initialize parameters according to Model_flags
    num = numel(N(:,:,:,1)); %total voxels in domain after resizing
    
    %Number of calibrations
    if(model_flags(1) == 1)
        ntp_cal = mri_cnt-1;
        ntp_pred = 0;
        
        tspan = t(2:end);
        
    elseif(model_flags(1)==2)
        ntp_cal = mri_cnt-2;
        ntp_pred = 1;
        
        tspan = t(2:end-1);
        tspan_pred = t(end);
    end
    N_true = N(:,:,:,1:ntp_cal+1);
    
    %Initial Parameter Guesses
    d_g  = 0.01; num_d = 1;
    if(model_flags(2)==1)
        kp_g = 0.01; %Global proliferation
        num_kp = 1;
    elseif(model_flags(2)==2)
        kp_g  = zeros(sy,sx,sz); %Local proliferation
        
        %Find extent of calibrated tumor
        temp = sum(N_true,4);
        idx_fitKp = find(temp>0);
        
        kp_g(idx_fitKp) = 0.1;
        
        num_kp = numel(idx_fitKp);
    end
    
    %Set number of fit parameters
    p = num_d + num_kp;
    
    %Parameter tracking variables
    kp_track    = kp_g(:);
    d_track     = d_g;
    alpha_track = alpha;
    
    
    %% Functional Block 3 - Parameter optimization
    disp('Initialize model SSE...');
    %Initialize Model Fit with parameter guesses
    % ntp_cal = num simulated time points; ntp_pred = num predictions
    
    N0 = N_true(:,:,:,1);
    N_meas = N_true(:,:,:,2:end);
    N_g = RDFDM_3D(N0,d_g,kp_g,dx,dy,dz,dt,tspan,bcs,theta);
  
    SSE = sum((N_meas - N_g).^2,'all');
    
    %Parpool if local proliferation is on
    if(model_flags(2)==2)
        pp = gcp('nocreate');
        if(isempty(pp))
           parpool 
        end
    end
    
    %Set parameter bounds for diffusion, proliferation 
    d_up   = 0.09999;
    d_low  = 1e-3;
    
    kp_up  = 0.09999;
    kp_low = -.05;
    

    %Normalization applied or not, helps with scaling differences between parameters
    if(normalize==1)
        norm_kp = kp_up;
        norm_d = d_up;
    else
        norm_kp   = 1;
        norm_d    = 1;
    end

    %Apply normalization to bounds
    kp_up = kp_up/norm_kp;
    d_up = d_up/norm_d;

    kp_low = kp_low/norm_kp;
    d_low = d_low/norm_d;

    %Normalize to start
    kp_g = kp_g./norm_kp;
    d_g = d_g./norm_d;
    
    iteration = 1;    
    stuck_check = 0; %Iterations since last successful update 
    disp('Starting calibration...');
    while iteration < max_it && SSE>e_tol %Stop on error minimization or max iterations
        %Jacobian = f((p*del) - f(p)) / (p*del - p)
        if(j_change==j_freq)
            J = zeros(num*ntp_cal, p);
            %Unnormalize for forward evaluations
            kp_g = kp_g.*norm_kp;
            d_g = d_g.*norm_d;

            %Set jacobian for changing Kp parameters
            %Renormalize kp for perturbation
            kp_g = kp_g./norm_kp;
            if(model_flags(2)==1) % Global kp
                %Forward Eval
                kp_test = kp_g*delta; 
                dif = kp_test - kp_g;
                %Unnormalize test Kp
                kp_test = kp_test.*norm_kp;

                model_kp = RDFDM_3D(N0,d_g,kp_test,dx,dy,dz,dt,tspan,bcs,theta);

                J(:,1) = reshape((model_kp - N_g),[],1)./(dif);

            elseif(model_flags(2)==2) %Local proliferations in largest extent
                parfor i = 1:num_kp
                    %Forward Eval
                    kp_test = kp_g; kp_test(idx_fitKp(i)) = kp_test(idx_fitKp(i))*delta;
                    dif = kp_test(idx_fitKp(i)) - kp_g(idx_fitKp(i));
                    %Unnormalize test Kp
                    kp_test = kp_test.*norm_kp;

                    model_kp = RDFDM_3D(N0,d_g,kp_test,dx,dy,dz,dt,tspan,bcs,theta);

                    J(:,i) = reshape((model_kp - N_g),[],1)./(dif);

                end
            end
            %End Kp perturbs, unnormalize curr best fit
            kp_g = kp_g.*norm_kp;

            %Normalize D for perturbs
            d_g = d_g./norm_d;

            %Set jacobian for changing D parameters
            d_test = d_g*delta;
            dif = d_test(1) - d_g(1);
            d_test = d_test.*norm_d;

            model_d = RDFDM_3D(N0,d_test,kp_g,dx,dy,dz,dt,tspan,bcs,theta);

            J(:,num_kp+1) = reshape((model_d - N_g),[],1)./(dif);

            %End d perturbs, unnormalize curr best fit
            d_g = d_g.*norm_d;
            
            %%%              Add extra parameter perturbations here                %%%
            %%% Can replace individual names with beta vector of calibrated params %%%

            j_change = 0; %Do not build again until J_freq updates have occured

            %Renormalize after jacobian complete
            kp_g = kp_g./norm_kp;
            d_g = d_g./norm_d;

            clear model_kp model_d
        end


        %Set residuals into column vector
        residuals = reshape((N_meas - N_g), [], 1);

        [update,~] = bicgstab(J'*J + alpha*diag(diag(J'*J)), J'*residuals, 1e-10,100);

        %Proliferation Update
        if(model_flags(2)==1) %Global calibration
            kp_check = kp_g + update(1);
            if(kp_check < kp_low)
                kp_check = kp_low;
            elseif(kp_check > kp_up)
                kp_check = kp_up;
            end
        elseif(model_flags(2)==2) %Local calibration in tumor extent
            kp_changing = kp_g(idx_fitKp) + update(1:num_kp);
            kp_changing(kp_changing<kp_low) = kp_low;
            kp_changing(kp_changing>kp_up)  = kp_up;

            kp_check = zeros(sy,sx,sz);
            kp_check(idx_fitKp) = kp_changing;
        end

        %Diffusivity Update
        d_check = d_g + update(num_kp + 1);
        if(d_check < d_low)
            d_check = d_low;
        elseif(d_check > d_up)
            d_check = d_up;
        end

        %Unnormalize params for forward evaluation
        kp_check = kp_check.*norm_kp;
        d_check = d_check.*norm_d;

        %Perform RD evaluation with testing parameters
        N_test = RDFDM_3D(N0,d_check,kp_check,dx,dy,dz,dt,tspan,bcs,theta);
        
        temp = numel(find(N_test>theta));
        if(temp~=0)
            disp(['# Larger than theta = ',num2str(temp)]);
        end
        
        temp = numel(find(isnan(N_test)));
        if(temp~=0)
            disp(['# NaNs = ',num2str(temp)]);
        end

        %Set new trial SSE
        SSE_test = sum((N_meas - N_test).^2,'all');

        %Normalize after forward eval
        kp_check = kp_check./norm_kp;
        d_check = d_check./norm_d;


        %Update successful if SSE reduces, change best fit parameters
        if(SSE_test <= SSE)
            %Updates Sim to new best calibration
            N_g = N_test;

            %Updates all parameters to new best fit
            kp_g = kp_check;
            d_g  = d_check;

            alpha = alpha/pass; %Move towards gauss-newton step

            j_change = j_change+1; %Successful update, increment j_change

            %Check if local convergence has been reached
            if(abs(SSE_test-SSE) < e_change)
               disp(['Stopped on iteration: ',num2str(iteration),', change in error less than ',num2str(e_change)]); 
               %Unnormalize for storing
               kp_g = kp_g.*norm_kp;
               d_g = d_g.*norm_d;

               d_track = [d_track, d_g];
               kp_track = [kp_track, kp_g(:)];
               alpha_track = [alpha_track, alpha];

               %Renormalize after
               kp_g = kp_g./norm_kp;
               d_g = d_g./norm_d;

               break;
            end    

            SSE = SSE_test;

        else %Update not successful
            alpha = alpha*fail; %Move towards gradient descent

            if(stuck_check == 0 && alpha > alpha_max) %Move to minimum if checked full alpha range
                alpha = alpha_min;
                stuck_check = stuck_check+1;
                %Stuck so reset jacobian
                j_change = j_freq;
            elseif(stuck_check == 1 && alpha > alpha_max) %Already tried full alpha range, no update found
                d_track = [d_track, d_g];
                kp_track = [kp_track, kp_g(:)];
                alpha_track = [alpha_track, alpha];

                disp(['Converged on iteration: ',num2str(iteration),', Full mu range tested with no change'])
                break;
            end
        end

        if(mod(iteration,50) == 0)
            disp(['Iteration: ', num2str(iteration)]);
        end

        if(SSE<e_tol)
            disp(['Error meets threshold, less than ', num2str(e_tol)]);
        end
        %Add current parameters to tracking variables for diagnsostic plotting
        %Unnormalize
        kp_g = kp_g.*norm_kp;
        d_g = d_g.*norm_d;

        d_track = [d_track, d_g];
        kp_track = [kp_track, kp_g(:)];
        alpha_track = [alpha_track, alpha];

        %Normalize new params or renormalize old
        kp_g = kp_g./norm_kp;
        d_g = d_g./norm_d;
        
        iteration = iteration+1;
    end %End of while Loop
    
    %% Post calibration
    kp_g = kp_g.*norm_kp;
    d_g  = d_g*norm_d;
    
    %Make prediction if needed
    if(ntp_pred~=0)
        N_pred = RDFDM_3D(N0,d_g,kp_g,dx,dy,dz,dt,tspan_pred,bcs,theta);
    else
        N_pred = [];
    end

    Outputs_cal.Sim = N_g;
    Outputs_cal.Pred = N_pred;
    
    %Save measured tumor details for visualization
    Outputs_cal.meas.N = N;
    Outputs_cal.meas.bcs = bcs;
    Outputs_cal.meas.dx = dx;
    Outputs_cal.meas.dy = dy;
    Outputs_cal.meas.dz = dz;


    Params_cal.kp = kp_g;
    Params_cal.d = d_g;

    
    %% Fit figures and save
    figure
    subplot(1,3,1)
    n = numel(d_track);
    plot(1:n, d_track); title('D Parameter Fit'); xlabel('Iteration'); ylabel('Diffusion (mm^2/day)');
    
    subplot(1,3,2)
    plot(1:n, kp_track); title('Kp Parameter Fits'); xlabel('Iteration'); ylabel('Proliferation (1/day)');
    
    subplot(1,3,3)
    semilogy(1:n, alpha_track(:,1:n)); title('Alpha Tracking'); xlabel('Iteration'); ylabel('Value');
end

function ccc = overlapCCC(m,sim,thresh)
    N_idx = find(m>0);
    Sim_idx = find(sim>thresh);
    idx = intersect(N_idx,Sim_idx);
    
    a = m(idx);
    b = sim(idx);
    
    mu_a = mean(a);
    mu_b = mean(b);

    cov = 0;
    varx = 0;
    vary = 0;
    for i = 1:length(idx)
       cov = cov + ((a(i)-mu_a)*(b(i)-mu_b)); 
       varx = varx + (a(i)-mu_a)^2;
       vary = vary + (b(i)-mu_b)^2;
    end
    cov = cov/length(idx);
    varx = varx/length(idx);
    vary = vary/length(idx);

    ccc = 2*cov/(varx+vary+(mu_a - mu_b)^2);
end

function dice = calcDice(m,sim,thresh)
    N_idx = find(m>0);
    Sim_idx = find(sim>thresh);
    idx = intersect(N_idx,Sim_idx);
    dice = 2*numel(idx)/(numel(N_idx)+numel(Sim_idx));
end

function V = getVol(IM,thresh,dx,dy,dz)
    test = zeros(size(IM));
    test(IM>thresh) = 1;
    test = imfill(test,'holes');
    V = numel(find(test))*dx*dy*dz;
end

function C = getCells(IM,thresh)
    C = sum(IM(IM>thresh),'all');
end