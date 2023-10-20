%{
Example script for the reaction-diffusion calibration with cohort sampling

Calibrates a scalar diffusivity (d) and proliferation (kp) for:
    dN/dt = grad * (d * grad(N)) + kp*N*(1-N/theta)

Can be modified to calibrate spatially varying prolfieration rate - kp(x)
or to include treatment.
    Calibration and forward model would need to be adjusted depending on 
    format of the radiation time course
%}
clear
clc
close all

k_bounds = [-0.05, 0.1];
d_bounds = [1e-6, 0.1];
dose_bounds = [0, 200];

num_patients = 10;

%Build example distributions for cohort parameters
%D is normally distributed
d_mu = mean(d_bounds);
d_sigma = (d_bounds(2) - d_mu)/2;
d_pd = makedist('normal', 'mu', d_mu, 'sigma', d_sigma);
d_pd = truncate(d_pd, d_bounds(1), d_bounds(2));

d_array = random(d_pd, num_patients, 1);

%Proliferation follows exponential distribution with dose absorbed
dose_array = dose_bounds(1) + rand(num_patients,1).*diff(dose_bounds);
fun = @(t) (k_bounds(2) - k_bounds(1)).*exp(-t*0.005);
kp_array = (fun(dose_array) + k_bounds(1))+randn(num_patients,1)*1e-2; %Creates noisy kp data


tumor = loadPatient('example_patient\', 1, [2,2,2]); %All tumors start from same N0
N0 = tumor.N(:,:,:,1);
%Set params for fwd runs
theta = tumor.theta; 
dx = tumor.dx; dy = tumor.dy; dz = tumor.dz;
tspan = tumor.t(2:end);
dt = 0.25;
bcs = tumor.bcs;

num_samples = 10;

f1 = figure; %Visualizes exponential fits with true parameter vs fits
f2 = figure; %Visualize cell count plots

y_max = 0;

for n = 1:num_patients
    %Get true parameters to build measured data
    d_true = d_array(n);
    kp_true = kp_array(n);
    dose_true = dose_array(n);
    N_true = RDFDM_3D(N0, d_true, kp_true, dx, dy, dz, dt, tspan, bcs, theta);
    
    %Set training distributions
    d_training = d_array; d_training(n) = [];
    kp_training = kp_array; kp_training(n) = [];
    dose_training = dose_array; dose_training(n) = [];
    
    %Fit training kp to exponential decay to get proliferation samples
    fkp = fit(dose_training, kp_training, 'exp1');
    CI = confint(fkp);
    coeff = coeffvalues(fkp);
    a_STD = sqrt(1) * (max(CI(:,1)) - min(CI(:,1)))/4;
    b_STD = sqrt(1) * (max(CI(:,2)) - min(CI(:,2)))/4;
    a_dist = makedist('Normal','mu',coeff(1),'sigma',a_STD);
    b_dist = makedist('Normal','mu',coeff(2),'sigma',b_STD);
    a_samp = random(a_dist, num_samples, 1);
    b_samp = random(b_dist, num_samples, 1);
    
    kp_samples = zeros(num_samples,1);
    for i = 1:num_samples
        fh = @(x) a_samp(i)*exp(b_samp(i)*x);
        kp_samples(i) = fh(dose_true); %Use true dose to get proliferation samples
    end
    
    %Fit training d to normal to get diffusivity samples
    mu_d = mean(d_training);
    std_d = std(d_training);
    pd = makedist('Normal','mu',mu_d,'sigma',std_d);
    pd = truncate(pd, d_bounds(1), d_bounds(2));
    
    d_samples = random(pd, num_samples, 1);
    
    %Predict for first time point
    Cell_TC_pred1 = zeros(num_samples, tspan(1)/dt + 1);
    for i = 1:num_samples
        [~, temp_TC] = RDFDM_3D(N0, d_samples(i), kp_samples(i), dx, dy, dz, dt, tspan(1), bcs, theta);
        Cell_TC_pred1(i,:) = squeeze(sum(sum(sum(temp_TC,3),2),1));
    end
    
    %Typically would calibrate parameters with current data at this point
    %%%% Baseline to first visit
    %We will just use the true parameters in the weighting scheme
    
    d_samples_2  = 0.5*d_true + 0.5*d_samples;
    kp_samples_2 = 0.5*kp_true + 0.5*kp_samples;
    
    %Predict for second time point, starting at most recent measured time point
    Cell_TC_pred2 = zeros(num_samples, (tspan(2)-tspan(1))/dt + 1);
    for i = 1:num_samples
        [~, temp_TC] = RDFDM_3D(N_true(:,:,:,1), d_samples_2(i), kp_samples_2(i), dx, dy, dz, dt, tspan(2)-tspan(1), bcs, theta);
        Cell_TC_pred2(i,:) = squeeze(sum(sum(sum(temp_TC,3),2),1));
    end
    
    %%% Would predict for third time point and further here if possible %%%
    
    %Visualize fit ditributions for proliferation
    set(0,'currentfigure',f2);
    subplot(2,5,n)
    scatter(dose_training, kp_training, 50, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 'Training');
    hold on
    scatter(dose_true,kp_true, 50, 'filled', 'MarkerFaceColor', 'black', 'DisplayName', 'True Testing');
    p1=plot(fkp);
    set(p1, 'Color','green'); set(p1, 'DisplayName','Fit'); set(p1, 'LineStyle','-');
    CI = predint(fkp,dose_training,0.95,'functional','on');
    p2=plot(sort(dose_training),sort(CI,'descend'));
    set(p2, 'Color','green'); set(p2, 'DisplayName','Fit - 95% CI'); set(p2, 'LineStyle','--');
    errorbar(dose_true,mean(kp_samples),std(kp_samples),'Color','red','LineStyle','none','Marker','o','MarkerFaceColor','red','DisplayName','Samples');
    legend off;
    if(n==10)
        legend;
    end
    xlabel('Dose (Gy)'); ylabel('Proliferation (1/day)'); title(['Patient ',num2str(n)]);
    
    %Visualize hurricane plots
    set(0,'currentfigure',f1);
    mu_pred1 = mean(Cell_TC_pred1(:,1:end-1),1);
    std_pred1 = std(Cell_TC_pred1(:,1:end-1),1);
    
    mu_pred2 = mean(Cell_TC_pred2,1);
    std_pred2 = std(Cell_TC_pred2,1);
    
    CI_pred1 = 1.96.*std_pred1./sqrt(num_samples);
    CI_pred2 = 1.96.*std_pred2./sqrt(num_samples);
    
    cells_true = [sum(N0,'all'), sum(N_true(:,:,:,1),'all'), sum(N_true(:,:,:,2),'all')];
    
    t1 = 0:dt:tspan(1)-dt;
    t2 = tspan(1):dt:tspan(2);
    
    subplot(2,5,n)
    %First interval
    fill([t1, fliplr(t1)], [mu_pred1 + CI_pred1, fliplr(mu_pred1 - CI_pred1)],[0,0,1], 'facealpha', 0.5);
    hold on
    plot(t1, mu_pred1, '-b');
    %Second interval
    fill([t2, fliplr(t2)], [mu_pred2 + CI_pred2, fliplr(mu_pred2 - CI_pred2)],[0,0,1], 'facealpha', 0.5);
    plot(t2, mu_pred2, '-b');
    %Measured data
    scatter(tumor.t, cells_true, 100,'filled','MarkerFaceColor','black');
    xlabel('Time Days'); ylabel('Cell Count'); title(['Patient ',num2str(n)]);
    
    lim = ylim;
    y_max = max([y_max, lim(2)],[],'all');
    
    disp(['Patient ', num2str(n), ' complete...']);
end
%reset figure limits
set(0,'currentfigure',f1);
for i = 1:num_patients
    subplot(2,5,i)
    ylim([0,y_max]);
end


