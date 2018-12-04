%% list subject data files:
group_subjects = {...
'Data_S027_03262018_E1.mat'...
'Data_S028_03262018_E1.mat'...
'Data_S029_03272018_E1.mat'...
'Data_S030_03272018_E1.mat'...
'Data_S031_03272018_E1.mat'...
'Data_S032_03272018_E1.mat'...
'Data_S033_03272018_E1.mat'...
'Data_S034_03272018_E1.mat'...
'Data_S035_03272018_E1.mat'...
'Data_S036_03272018_E1.mat'...
'Data_S037_03272018_E1.mat'...
'Data_S038_03282018_E1.mat'...
'Data_S039_03282018_E1.mat'...
'Data_S040_03282018_E1.mat'...
'Data_S041_03282018_E1.mat'...
'Data_S042_03282018_E1.mat'...
'Data_S043_04042018_E1.mat'...
'Data_S044_04042018_E1.mat'...
'Data_S045_04042018_E1.mat'...
'Data_S071_04092018_E1.mat'...
'Data_S073_04102018_E1.mat'...
};

%% define some constants and parameters:
EARLIEST_VALID_PT = -.2;
LATEST_VALID_PT = .85;

%% Analyze each subject:

indiv_error_queue = {};

a = [];
all_timer = tic;
p_bins = nan(8, length(group_subjects), 3);
z_bins = nan(6, length(group_subjects), 2);
view_time = nan(70, length(group_subjects), 3);
min_pt = nan(1, length(group_subjects));
dir_error = nan(70, length(group_subjects), 3);
dir_absolute = nan(70, length(group_subjects), 3);
catch_tr_apt = nan(12, length(group_subjects), 2);
catch_tr_ppt = nan(12, length(group_subjects), 2);
catch_tr_de = nan(12, length(group_subjects), 2);

h1 = figure; h2 = figure;
target_distances = nan(1, length(group_subjects));
for i_sub = 1:length(group_subjects)
    sub_timer = tic;
    try
        load(group_subjects{i_sub})
        data_indiv = analyze_retention_individual_v1(Data, 0);
        
        p_bin0 = data_indiv.p_bins{1};
        p_bin1 = data_indiv.p_bins{2};
        p_bin2 = data_indiv.p_bins{3};
        pt_z_bin0 = data_indiv.z_bins{1};
        pt_z_bin2 = data_indiv.z_bins{2};
        
        p_bins(:, i_sub, 1) = p_bin0;
        p_bins(:, i_sub, 2) = p_bin1;
        p_bins(:, i_sub, 3) = p_bin2;

        z_bins(:, i_sub, 1) = pt_z_bin0;
        z_bins(:, i_sub, 2) = pt_z_bin2;

        min_pt(i_sub) = data_indiv.min_pt;
        
        type0 = data_indiv.types{1};
        type1 = data_indiv.types{2};
        type2 = data_indiv.types{3};
        type3 = data_indiv.types{4};
        type4 = data_indiv.types{5};
        
        Dir_e = data_indiv.Dir_e;
        Dir_a = data_indiv.Dir_a;
        
        Kin_x = data_indiv.kinematics{1};
        Kin_y = data_indiv.kinematics{2};
        
        vt_temp = Data.ViewTime(Data.Type == 0);
        de_temp = Dir_e(Data.Type == 0);
        da_temp = Dir_a(Data.Type == 0);
        view_time(1:length(vt_temp), i_sub, 1) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 1) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 1) = da_temp;

        vt_temp = Data.ViewTime(Data.Type == 1);
        de_temp = Dir_e(Data.Type == 1);
        da_temp = Dir_a(Data.Type == 1);
        view_time(1:length(vt_temp), i_sub, 2) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 2) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 2) = da_temp;

        vt_temp = Data.ViewTime(Data.Type == 2);
        de_temp = Dir_e(Data.Type == 2);
        da_temp = Dir_a(Data.Type == 2);
        view_time(1:length(vt_temp), i_sub, 3) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 3) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 3) = da_temp;
        
        catch_tr_apt(1:length(type3), i_sub, 1) = Data.ViewTime(type3);
        catch_tr_apt(1:length(type4), i_sub, 2) = Data.ViewTime(type4);
        
        catch_tr_ppt(1:length(type3), i_sub, 1) = Data.pPT(type3);
        catch_tr_ppt(1:length(type4), i_sub, 2) = Data.pPT(type4);
        
        catch_tr_de(1:length(type3), i_sub, 1) = Dir_e(type3);
        catch_tr_de(1:length(type4), i_sub, 2) = Dir_e(type4);

        open_winds = findobj('type', 'figure');
        for i_wind = 1:(length(open_winds)-2)
            close(i_wind + 2);
        end
        figure(h1); 
        subplot(ceil(sqrt(length(group_subjects))), ceil(sqrt(length(group_subjects))), i_sub); hold on;
        plot(Data.ViewTime(type0), Dir_e(type0) , '.', 'MarkerSize', 16, 'Color', [172, 59, 59]/255)
        plot(Data.ViewTime(type1), Dir_e(type1), '.', 'MarkerSize', 15, 'Color', [85, 170, 85]/255)
        plot(Data.ViewTime(type2), Dir_e(type2), '.', 'MarkerSize', 14, 'Color', [86/255 85/255 149/255])
        plot(Data.ViewTime(type3), Dir_e(type3), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [225, 244, 162]/255)
        plot(Data.ViewTime(type4), Dir_e(type4), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [71, 113, 134]/255)
        plot([-.2 .8], [45 45], 'k')
        plot([-.2 .8], -[45 45], 'k')
        plot([min_pt(i_sub), min_pt(i_sub)], [-200 200], 'k-');
        axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
        
        figure(h2);
        subplot(ceil(sqrt(length(group_subjects))), ceil(sqrt(length(group_subjects))), i_sub); hold on;
        plot(Kin_x(:,type0), -Kin_y(:,type0), '-', 'LineWidth', 2, 'Color', [172, 59, 59]/255)
        plot(Kin_x(:,type1), -Kin_y(:,type1), '-', 'LineWidth', 1, 'Color', [85, 170, 85]/255)
        plot(Kin_x(:,type2), -Kin_y(:,type2), '-', 'LineWidth', .5, 'Color', [86/255 85/255 149/255])
        axis([-20 20 -20 20])
        
        target_distances(i_sub) = data_indiv.targ_len;
    catch individ_err
        warning(['Subject ', num2str(i_sub), ' failed.'])
        indiv_error_queue{length(indiv_error_queue) + 1} = individ_err;
    end
    a(i_sub) = toc(sub_timer);
end
toc(all_timer)

%% plot group aggregate results:
n_bins = 8;
vt = view_time(:);
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
[n_pt_all, edge_pt_all] = histcounts(vt, hist_bins);
x_ind = edge_pt_all(2:end) - diff(edge_pt_all)/2;

figure; hold on;
errorbar(x_ind + 0.02, nanmean(p_bins(:,:,1), 2), sqrt(nanvar(p_bins(:,:,1),0,2)./sum(~isnan(p_bins(:,:,1)),2)), '.-', 'MarkerSize', 20, 'Color', [172, 59, 59]/255);
errorbar(x_ind - 0.02, nanmean(p_bins(:,:,2), 2), sqrt(nanvar(p_bins(:,:,2),0,2)./sum(~isnan(p_bins(:,:,2)),2)), '.-', 'MarkerSize', 20, 'Color', [85, 170, 85]/255);
errorbar(x_ind, nanmean(p_bins(:,:,3), 2), sqrt(nanvar(p_bins(:,:,3),0,2)./sum(~isnan(p_bins(:,:,3)),2)), '.-', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255]);
set(gca,'fontsize',20)
legend('None', 'Direct', 'Symbolic')

n_bins = 6; hist_bins = linspace(-.005, .605, n_bins+1);
[pt_sort, i_sort] = sort(Data.pPT);
[n_ppt_all, edge_ppt_all] = histcounts(pt_sort, hist_bins);
x_ind = edge_ppt_all(2:end) - diff(edge_ppt_all)/2;
figure; hold on;
errorbar(x_ind + 0.01, nanmean(z_bins(:,:,1), 2), sqrt(nanvar(z_bins(:,:,1),0,2)./length(group_subjects)), '.-', 'MarkerSize', 20, 'Color', [172, 59, 59]/255);
errorbar(x_ind - 0.01, nanmean(z_bins(:,:,2), 2), sqrt(nanvar(z_bins(:,:,2),0,2)./length(group_subjects)), '.-', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255]);
set(gca,'fontsize',20)
legend('None', 'Symbolic')