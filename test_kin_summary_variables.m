% script to do test of kinematic summary variables:
%
% Are the kinematics different across trial conditions and PT?
%
% David Huberdeau, NTB lab, 07/11/2019

data_pv = nan(size(type_all, 2), 3); %subjects X types
data_mt = nan(size(type_all, 2), 3); %subjects X types
data_var = nan(size(type_all, 2), 3); %subjects X types

types = [0,1,2];
for i_sub = 1:size(type_all,2)
    for i_type = 1:3
        type_temp = type_all(:, i_sub) == types(i_type);

        data_pv(i_sub, i_type) = nanmean(kin_PV(type_temp, i_sub));
        data_mt(i_sub, i_type) = nanmean(kin_MT(type_temp, i_sub));
        data_var(i_sub, i_type) = nanmean(kin_VAR(type_temp, i_sub));
    end
end

%% do 1-way anova stat test:
anova1(data_pv);
anova1(data_mt);
anova1(data_var);

%% for 2-way anova (factoring by PT):

data_pv_low = nan(size(type_all, 2), 3); %subjects X types
data_mt_low = nan(size(type_all, 2), 3); %subjects X types
data_var_low = nan(size(type_all, 2), 3); %subjects X types

data_pv_high = nan(size(type_all, 2), 3); %subjects X types
data_mt_high = nan(size(type_all, 2), 3); %subjects X types
data_var_high = nan(size(type_all, 2), 3); %subjects X types

types = [0,1,2];
for i_sub = 1:size(type_all,2)
    for i_type = 1:3
        type_temp = type_all(:, i_sub) == types(i_type) & ...
            pt_all(:, i_sub) < min_pt(i_sub);

        data_pv_low(i_sub, i_type) = nanmean(kin_PV(type_temp, i_sub));
        data_mt_low(i_sub, i_type) = nanmean(kin_MT(type_temp, i_sub));
        data_var_low(i_sub, i_type) = nanmean(kin_VAR(type_temp, i_sub));
        
        type_temp = type_all(:, i_sub) == types(i_type) & ...
            pt_all(:, i_sub) > min_pt(i_sub);

        data_pv_high(i_sub, i_type) = nanmean(kin_PV(type_temp, i_sub));
        data_mt_high(i_sub, i_type) = nanmean(kin_MT(type_temp, i_sub));
        data_var_high(i_sub, i_type) = nanmean(kin_VAR(type_temp, i_sub));
    end
end

%%
anova2([data_pv_low; data_pv_high], 2);
anova2([data_mt_low; data_mt_high], 2);
anova2([data_var_low; data_var_high], 2);