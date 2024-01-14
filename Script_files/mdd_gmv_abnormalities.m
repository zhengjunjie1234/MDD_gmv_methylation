%% loading data
clear
clc

load('../Data_files/hc_gmv.mat')
load('../Data_files/mdd_gmv.mat')
load('../Data_files/HCP360_regions_name.mat')
load('../Data_files/age_sex_edu_data.mat')
load('../Data_files/tiv_data.mat')

disp('loading data finished.....')
%% extract left hemisphere regional gmv
for j=1:180
    for i=1:length(hc_gmv)
        hc_gmv_l(i,j) = hc_gmv{i}(j);
    end
    for i=1:length(mdd_gmv)
        mdd_gmv_l(i,j) = mdd_gmv{i}(j);
    end
end
gmv_all= [hc_gmv_l;mdd_gmv_l];

disp('extract data finished.....')
%% collect covariates 
age_sex_edu_all = [age_sex_edu_hc;age_sex_edu_mdd];
tiv_all = [tiv_hc;tiv_mdd];
cov_all = [age_sex_edu_all,tiv_all];

disp('collect covariates finished.....')
%% glm statistics t-test

n1 = size(hc_gmv_l,1);
n2 = size(mdd_gmv_l,1);
label_all = [];
label_all(1:n1,1)=0;
label_all(n1+1:n1+n2,1)=1;

x=[label_all,cov_all];
for i = 1:180
    y= gmv_all(:,i);
    [bb,dev,stats] = glmfit(x,y);
    t_glm(i) = stats.t(2);
    p_glm(i) = stats.p(2);
    y_res= y - x(:,2:end)*bb(3:end);
    gmv_all_r_glm(:,i) = y_res; 
end
p_fdr = mafdr(p_glm,'BHFDR',true);
regions_mdd_hc = find(p_fdr<0.05);
regions_mdd_hc_name = regions_net(regions_mdd_hc,:)
save('gmv_all_r.mat','gmv_all_r_glm','regions_mdd_hc');

disp('glm statistics t-test finished.....')
%% abnormal gmv t values and p values
t_glm_sig = t_glm(regions_mdd_hc);
p_fdr_sig = p_fdr(regions_mdd_hc);
save('t_values_pls.mat','t_glm');
%% glm medication covarates validation
load('../Data_files/medication_effects.mat')
med_all = [med_hc;med_mdd];

x=[label_all,age_sex_edu_all,med_all];
for i = 1:180
    y= gmv_all(:,i);
    [bb,dev,stats] = glmfit(x,y);
    t_glm_med(i) = stats.t(2);
    p_glm_med(i) = stats.p(2);
end

t_glm_med_sig = t_glm_med(regions_mdd_hc);
p_glm_med_sig = p_glm_med(regions_mdd_hc);

%%


