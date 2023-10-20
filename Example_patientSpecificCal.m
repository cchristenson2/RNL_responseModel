%{
Example script for the reaction-diffusion calibration

Calibrates a scalar diffusivity (d) and proliferation (kp) for:
    dN/dt = grad * (d * grad(N)) + kp*N*(1-N/theta)

Can be modified to calibrate spatially varying prolfieration rate - kp(x)
or to include treatment.
    Calibration and forward model would need to be adjusted depending on 
    format of the radiation time course
%}
clear
clc

%setup for problem
loc = 'example_patient\'; %in silico tumor with noise added

model = 1; %M0 model is what the data is prepared for
%if model = 2; uses parpool for local proliferation jacobian (not recommended without HPC)

cal_flag = 2; %Calibrates to first time point only, saves last for prediction
%Use 1 for full time course calibration

model_flags = [cal_flag, model];

max_it = 100;

reduce = 1; %Simulates on a smaller domain surrounding the tumor location

%Call calibration
[Params_cal, Outputs_cal] = Calibrate_RDFDM_3D_example(loc, max_it, model_flags, reduce);

%Load true values
Params_true = load([loc,'Example_true_params.mat']);


%Visualize calibration
[sy,sx,sz,~] = size(Outputs_cal.Sim);
slice = round(sz/2);
theta = 818503*Outputs_cal.meas.dx*Outputs_cal.meas.dy*Outputs_cal.meas.dz;

figure
subplot(2,3,1)
imagesc(squeeze(Outputs_cal.meas.N(:,:,slice,2)), [0,theta]); 
axis image; xticks([]); yticks([]); ylabel('Calibration Day 28'); title('Measured');

subplot(2,3,2)
imagesc(squeeze(Outputs_cal.Sim(:,:,slice,1)), [0,theta]);
axis image; xticks([]); yticks([]); title('Simulation');

subplot(2,3,3)
temp_meas = Outputs_cal.meas.N(:,:,slice,2); temp_sim = Outputs_cal.Sim(:,:,slice,1);
idx = intersect(find(temp_meas), find(temp_sim));
scatter(temp_meas(idx), temp_sim(idx), 50, 'filled');
refline(1,0);
xlabel('Measured Cells'); ylabel('Calibrated Cells'); axis square

subplot(2,3,4)
imagesc(squeeze(Outputs_cal.meas.N(:,:,slice,3)), [0,theta]); 
axis image; xticks([]); yticks([]); 
if(cal_flag==1)
    ylabel('Calibration Day 56');
    temp_sim = squeeze(Outputs_cal.Sim(:,:,slice,2));
else
    ylabel('Prediction Day 56');
    temp_sim = squeeze(Outputs_cal.Pred(:,:,slice));
end


subplot(2,3,5)
imagesc(temp_sim, [0,theta]);
axis image; xticks([]); yticks([]);

subplot(2,3,6)
temp_meas = squeeze(Outputs_cal.meas.N(:,:,slice,3));
scatter(temp_meas(idx), temp_sim(idx), 50, 'filled');
refline(1,0);
xlabel('Measured Cells'); 
if(cal_flag==1)
    ylabel('Calibrated Cells'); 
else
    ylabel('Predicted Cells'); 
end
axis square


disp(['Diffusivity % error = ',num2str(100*(Params_cal.d-Params_true.d)/Params_true.d),'%']);
disp(['Proliferation % error = ',num2str(100*mean((Params_cal.kp-Params_true.kp)./Params_true.kp,'all')),'%']);
