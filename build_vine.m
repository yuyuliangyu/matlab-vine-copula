clc,clear;
load data.mat;
%%
%数据预处理，将不同时间序列合并，保证每个时间点在每个地区都有发生
vine_data=cell(1,5);
vine_data{1}=vine_sse;
vine_data{2}=vine_hsi;
vine_data{3}=vine_twii;
vine_data{4}=vine_n225;
vine_data{5}=vine_ks11;
dimension=length(vine_data);
for i=1:length(vine_data)
    data=vine_data{i};
    time=table2array(data(:,1));
    time=datenum(time);
    time=time(~any(isnan(time),2),:);%any(,2)指定行中是否有零元素。
    data=table2array(data(:,2));
    data=data(~any(isnan(data),2),:);
    vine_data{i}=fints(time,data, sprintf('vine%d', i));
end
combine_data =vine_data{1};
for i=1:length(vine_data)-1
    combine_data =merge(combine_data,vine_data{i+1},'DateSetMethod','Intersection');
end
%%
%找出合并的序列日期与韩国日期（最短）之间的差集
different_dates = setxor(combine_data.dates, vine_data{5}.dates);
date=datestr(different_dates);
%%
vine_returns=diff(log(combine_data));
vine_returns=fts2mat(vine_returns);

%%
%利用garch提取残差
dimension=size(vine_returns,2);
aicbic_matrix=[];
for i=1:dimension
    [estParams, ~, logL, info] = estimate(garch_model(vine_returns(:,i)), vine_returns(:,i));
    [aic,bic]= aicbic(logL,4,size(vine_returns(:,i),1));%aicbic输入的分别是似然函数，参数数量，数据规模
    aicbic_matrix = [aicbic_matrix;aic bic];
    vine_returns(:,i)=infer(garch_model(vine_returns(:,i)),vine_returns(:,i));
end
%%
%进行经验分布的累次积分变换，注意：如果经过累次积几乎不影响肯德尔tao值，已经实证证明。
dimension=size(vine_returns,2);
vine_return_union=[];
vine_return_index=[];
for i=1:dimension %计算每个收益率的经验分布
    [a b]=ecdf(vine_returns(:,i));
    vine_return_union=[vine_return_union a];
    vine_return_index=[vine_return_index b];
end
for i=1:length(vine_returns) %将每个收益率带入经验分布计算累次积分变换后的值
    for j=1:dimension
        a=max( find( vine_return_index(:,j)<=vine_returns(i,j)) );
        vine_returns(i,j)=vine_return_union(a,j);
    end
end
[a b]=find(vine_returns==1);
vine_returns(a,b)=0.99;%累次积分变换的值必须大于0小于1，这里的代码不允许取到0和1.
[a b]=find(vine_returns==0);
vine_returns(a,b)=0.0000001;
% t = trnd(5, 1000, 5);
%%
%对转换后的序列进行均匀分布检验
unionTest_matrix=[];
for i=1:size(vine_returns,2)
[h,p,s] = chi2gof(vine_returns(:,i), 'cdf', {@unifcdf,min(vine_returns(:,i)),max(vine_returns(:,i))});
unionTest_matrix=[unionTest_matrix;h,p,s.chi2stat]
end
%%
[index_matrix,copulafamily_matrix,parameters_matrix,vine_returns] = vine_construction(vine_returns);
[lik]=lik_vine(index_matrix,copulafamily_matrix,parameters_matrix,vine_returns);
% [lik]=lik_vine(index_matrix,copulafamily_matrix,[-0.676181327155250,0.519937870224427,0.280588434296862,11.5635041812958,1.03855149967934,7.81116046027813,-0.387275442720832,0.543849080217753,12.7596036451366,0.140092015817294],vine_returns);
%%
[parameter_max1,ln_lik_max]=optimal_parameters_vine(index_matrix,copulafamily_matrix,parameters_matrix,vine_returns);
% [parameter_max2,ln_lik_max]=optimal_parameters_vine(index_matrix,copulafamily_matrix,parameter_max1,vine_returns);