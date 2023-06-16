%% pls model analysis
clear
clc

load('t_values_pls.mat');
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat');

disp('finished......')
%% preprocessing : remove nan values and match regional gene expression

mri = t_glm';
[m,n]=find(isnan(mri(1:end,:)))
temp1=m;
temp2=find(isnan(parcelExpression(:,2)));
temp= union(temp1,temp2);
temp3=[181:360];
temp = union(temp,temp3);
region_ind=setdiff(parcelExpression(:,1),temp);
group_express=parcelExpression(region_ind,2:end);   
GENEdata=group_express;      
gene_name = probeInformation.GeneSymbol;

y=mri(region_ind,1:end);
MRIdata = zscore(y);
GENEdata = zscore(GENEdata);
disp('data transform finished......')

%% PLS_calculation
Y = MRIdata;
dim = 15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);
%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],MRIdata);
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
%% plot pls
dim = 15;
PVE = PCTVAR(2,:);

PVE_d = sort(PVE,'descend');

figure
plot(1:dim,PVE_d,'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14,'FontWeight','bold')
xlabel('Number of PLS components','FontSize',14,'FontWeight','bold');
ylabel('Percent Variance Explaination','FontSize',14,'FontWeight','bold');
ylim([0,0.3])
grid on

[~,n_max] = max(PVE);
XS_pls = XS;
Y1=Y;
figure
plot(XS_pls(:,n_max),Y1,'ro','MarkerEdgeColor',[140/255,0,0],'MarkerFaceColor',[140/255,0,0],'MarkerSize',4)
set(gca,'Fontsize',14,'FontWeight','bold')
lsline()
[R,p]=corrcoef(XS_pls(:,1),Y1) 
xlabel('PLS weights for PLS component 1','FontSize',14,'FontWeight','bold');
ylabel('MDD-HC GMD t-statistics ','FontSize',14,'FontWeight','bold');
grid on

%%
gene_name = probeInformation.GeneSymbol;
geneindex=1:size(GENEdata,2);
genes = gene_name;
bootnum=10000;
X=GENEdata;
Y=MRIdata;
dim=15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

n_max = 1;

[R1,p1]=corr(XS(:,n_max),MRIdata);
if R1(1,1)<0
    stats.W(:,n_max)=-1*stats.W(:,n_max);
    XS(:,n_max)=-1*XS(:,n_max);
end
[PLS1w,x1] = sort(stats.W(:,n_max),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
PLS1_ROIscores=XS(:,n_max);
save(['PLS1_ROIscore.mat'],'PLS1_ROIscores');
csvwrite(['PLS1_ROIscores.csv'],XS(:,n_max));
PLS1_score=XS(:,n_max);

PLS1weights = zeros(10027,bootnum);

res = zeros(bootnum,length(region_ind));

parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,n_max);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW; 
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run   
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);

zvalues = Z1;
pvalue = 2*(1-normcdf(abs(zvalues)));
pfdr = mafdr(pvalue,'BHFDR',true);

%%
fid1 = fopen(['PLS1_geneWeights.csv'],'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f,%f,%f\n', PLS1{i},geneindex1(i), Z1(i),pvalue(i),pfdr(i));
end
fclose(fid1);

disp('pls calculation finished......')

%%









