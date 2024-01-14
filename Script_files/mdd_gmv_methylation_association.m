%%
clear
clc

load('../Data_files/gmv_all_r.mat') 
load('../Data_files/gmv_id.mat') %%%% gmv samples id
load('../Data_files/methy_index.mat')  %%%%% methylation samples id
load('../Data_files/sig_methy_cpg_2.mat')
load('../Data_files/pathway_gene_targetpath.mat')
load('../Data_files/demethy_pls_genes05.mat')  %%% PLS and DMP overlapped genes 
load('../Data_files/cpg_site.mat') %%% PLS and DMP overlapped genes' cpg sites

disp('loading data finished.....')
%% get the matched samples having both GMV and DNA methylation data
all_ids = [hc_id_gmv;mdd_id_gmv];

sig_region = regions_mdd_hc;
%%% hc
[a,id_hc] = ismember(hc_index_methy,all_ids);
hc_gmv = gmv_all_r_glm(id_hc,sig_region);

%%% mdd
[a,id_mdd] = ismember(mdd_index_methy,all_ids);
mdd_gmv = gmv_all_r_glm(id_mdd,sig_region);

disp('matching data finished.....')
%% DNA methylation PCA model 
gmv_all = [mdd_gmv;hc_gmv];
sig_methy_data_all = [sig_methy_cpg_mdd;sig_methy_cpg_hc];

sig_methy_data_all_overlap = sig_methy_cpg_mdd(:,pathway_gene_targetpath==1);
sig_methy_data_all_overlap = zscore(sig_methy_data_all_overlap);


[coeff,score,latent,tsquared,explained,mu] = pca(sig_methy_data_all_overlap);

sum(explained(1:25)) %%%% explain rates;
sum_exp = [];
for i=1:25
    sumexp(i) = sum(explained(1:i));
end

dim=25;
figure
plot(1:dim,sumexp,'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PCA Components','FontSize',14);
ylabel('Percent of Explanation Variance','FontSize',14);
xlim([1,28])
grid on

disp('PCA analysis finished.....')
%% PCR analysis of mdd gmv and mdd methylation
F_step_gmv_methy = [];
P_step_gmv_methy = [];

for j=1:15
    k=10+j;
    [coeff,score,latent,tsquared,explained,mu] = pca(sig_methy_data_all_overlap,'NumComponents',k);
    for i=1:10
        [b,se,pval,finalmodel,stats] = stepwisefit(score,mdd_gmv(:,i),'display','off');
        F_step_gmv_methy(j,i) = stats.fstat;
        P_step_gmv_methy(j,i) = stats.pval;
    end 
end

disp('PCR model analysis finished.....')
%% get PCR compents of FFC , IFG, ACC
R = 1;

[b,se,pval,finalmodel,stats] = stepwisefit(score,mdd_gmv(:,R),'display','off');
pred_y_swlr1 = stats.intercept + score(:,finalmodel)*b(finalmodel);
[r,p] = corr(mdd_gmv(:,R),pred_y_swlr1)
figure
plot(mdd_gmv(:,R),pred_y_swlr1,'.')
lsline
comps_1 = find(finalmodel==1);
comps_1_v = b(finalmodel==1)


%%% IFG
R = 8;
[b,se,pval,finalmodel,stats] = stepwisefit(score,mdd_gmv(:,R),'display','off');
pred_y_swlr1 = stats.intercept + score(:,finalmodel)*b(finalmodel);
[r,p] = corr(mdd_gmv(:,R),pred_y_swlr1)
figure
plot(mdd_gmv(:,R),pred_y_swlr1,'.')
lsline
comps_2 = find(finalmodel==1);
comps_2_v = b(finalmodel==1)

%%% ACC
R = 10;
[b,se,pval,finalmodel,stats] = stepwisefit(score,mdd_gmv(:,R),'display','off');
pred_y_swlr1 = stats.intercept + score(:,finalmodel)*b(finalmodel);
[r,p] = corr(mdd_gmv(:,R),pred_y_swlr1)
figure
plot(mdd_gmv(:,R),pred_y_swlr1,'.')
lsline 
comps_3 = find(finalmodel==1);
comps_3_v = b(finalmodel==1)

disp('IFG, ACC, FFC PCR model analysis finished.....')
%% Leave one out cross valition of PCR model

%%% FFC
R = 1;
x=score(:,comps_1);
y=mdd_gmv(:,R);
y_pred = [];

for i=1:length(y)
    train_data_index = setdiff([1:length(y)],i);
    test_y = y(i);
    test_x = x(i,:);
    train_y = y(train_data_index);
    train_x = x(train_data_index,:);
    [b,bint,r,rint,stats] = regress(train_y,[ones(length(train_x),1),train_x]);
    pred_y = [1,test_x]*b;  
    y_pred(i) = pred_y;
end
[r,p] = corr(y,y_pred')
figure
plot(y,y_pred','.')
lsline


%%% IFG
R = 8;
x=score(:,comps_2);
y=mdd_gmv(:,R);
y_pred = [];

for i=1:length(y)
    train_data_index = setdiff([1:length(y)],i);
    test_y = y(i);
    test_x = x(i,:);
    train_y = y(train_data_index);
    train_x = x(train_data_index,:); 
    [b,bint,r,rint,stats] = regress(train_y,[ones(length(train_x),1),train_x]);
    pred_y = [1,test_x]*b; 
    y_pred(i) = pred_y;
end

[r,p] = corr(y,y_pred')
figure
plot(y,y_pred','.')
lsline


%%% ACC
R = 10;
x=score(:,comps_3);
y=mdd_gmv(:,R);
y_pred = [];


for i=1:length(y)
    train_data_index = setdiff([1:length(y)],i);
    test_y = y(i);
    test_x = x(i,:);
    train_y = y(train_data_index);
    train_x = x(train_data_index,:); 
    [b,bint,r,rint,stats] = regress(train_y,[ones(length(train_x),1),train_x]);
    pred_y = [1,test_x]*b;
    y_pred(i) = pred_y;
end

[r,p] = corr(y, y_pred')
figure
plot(y, y_pred', '.')
lsline


disp('IFG, ACC, FFC PCR LOOCV predict finished.....')
%% get methylation cpg weights in FFC,IFG,ACC  

%%% ffc
[r1,p] = corr(sig_methy_data_all_overlap,score(:,comps_1));
r_comp1= sum(abs(r1)/3,2);
r_comp1_r = sort(r_comp1,'descend');
r1_comp1 = r1;

%%%IFG
[r1,p] = corr(sig_methy_data_all_overlap,score(:,comps_2));
r_comp2 = sum(abs(r1)/2,2);
r_comp2_r = sort(r_comp2,'descend');
r2_comp1 = r1;

%%% ACC
[r1,p] = corr(sig_methy_data_all_overlap,score(:,comps_3));
r_comp3 = sum(abs(r1)/2,2);
r_comp3_r = sort(r_comp3,'descend');
r3_comp1 = r1;

weight_all = [];
weight_all = [r_comp1,r_comp2,r_comp3];
id1 = find(weight_all(:,1)>r_comp1_r(5));
id2 = find(weight_all(:,2)>r_comp2_r(10));
id3 = find(weight_all(:,3)>r_comp3_r(10));

ids = union(id1,id2);
ids = union(ids,id3); 
ids = union(ids,[24,28,120])  %%% add MDD risk genes 'CRHBP','HTR1A','NTRK3';

r1_comp1_ids = r1_comp1(ids,:);
r2_comp1_ids = r2_comp1(ids,:);
r3_comp1_ids = r3_comp1(ids,:);

weight_all_1 = weight_all(ids,:);
%weight_all_2 = weight_all_1./sum(weight_all_1);

genes_pathway = demethy_pls_genes05(pathway_gene_targetpath==1) ;
cpg_pathway = cpg_site(pathway_gene_targetpath==1);

overlap_genes = demethy_pls_genes05(pathway_gene_targetpath==1);
weight_genes = overlap_genes(ids)
weight_cpg = cpg_pathway(ids)

%%%% remove unrelated genes and the genes had multiple cpgs
weight_genes_sub = [3,5,6,7,8,10,15,16,20,21,22,24];
weight_genes(weight_genes_sub)

weight_all_1([4,9,11,13,8,17,23],:)=[];
weight_genes([4,9,11,13,8,17,23]) =[];
weight_cpg([4,9,11,13,8,17,23]) = [];
ids([4,9,11,13,8,17,23])=[];

disp('IFG, ACC, FFC PCR model top related cpg & genes finished.....')
%% cpg methy and gene exp correlation in IFG
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat');
group_express=parcelExpression(sig_region,2:end);   
GENEdata=group_express;     
gene_name = probeInformation.GeneSymbol;

%%% mean values of cpg methylation
sig_methy_cpg_mdd_overlap = sig_methy_cpg_mdd(:,pathway_gene_targetpath==1);
sig_methy_cpg_mdd_overlap_mean = mean(sig_methy_cpg_mdd_overlap)';

%%% IFG cpgs methylation and gene expression
%%% select cpg weights or gene express higher in IFG;
weight_genes_sub = [3:8,10,11,14:17]

[h,id] = ismember(weight_genes(weight_genes_sub),gene_name);
cpg_gene_exp = GENEdata(:,id);
score_methy1 = sig_methy_cpg_mdd_overlap_mean(ids(weight_genes_sub));

[r,p] = corr(cpg_gene_exp([8],:)',score_methy1,'type','Spearman')

genes_ifg_top = weight_genes(weight_genes_sub);
exp_ifg_cpg = cpg_gene_exp(8,:)';
methy_ifg_cpg = score_methy1;

plot(score_methy1,cpg_gene_exp(8,:)','.')
lsline 

disp('IFG related genes DNA methylation and Gene expression finished.....')
%% ifg random permuation test of correlation between cpg methylation and gene expression
[~,weight_genes_sort] = sort(weight_all(:,2),'descend');
weight_genes_sort_genes = overlap_genes(weight_genes_sort);

[h,id] = ismember(weight_genes_sort_genes,gene_name);

x = sig_methy_cpg_mdd_overlap_mean(weight_genes_sort);
y = GENEdata(8,id)';

r = [];
p=[];
for i=1:5000
    rand_id = randsample([1:length(x)],length(weight_genes_sub));
    [r(i),p(i)] = corr(x(rand_id,1),y(rand_id,1),'type','Spearman');
end
hist(r)
r = r';

p_permutation = sum(r<-0.77)/5000

%%
