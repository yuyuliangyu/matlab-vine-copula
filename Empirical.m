function [outputArg1] = Empirical(varargin)
%EMPIRICAL 此处显示有关此函数的摘要
%   此处显示详细说明
if isnumeric(varargin{1})%判断是否为第一层，因为第二层开始，信息就不仅仅有变量，还有别的信息
    inputArg1=varargin{1};
    dimension=size(inputArg1,2);
    empirical_T=zeros(dimension);
    for k=1:(dimension-1)
        for w=(k+1):dimension
            empirical_T(k,w)=corr(inputArg1(:,k),inputArg1(:,w),'Type','Kendall');
        end
    end
    outputArg1= empirical_T;
else
    outputArg1=cell(1,2);
    a=varargin{1,1};
    inputArg1=a(end,:);
    str=a(2,:);
    theta=a(3,:);
    a=a(1,:);
    dimension=size(a,2);
    empirical_T=zeros(dimension);
    %     b=varargin{2};
    for i=1:length(a)-1
        for j=i+1:length(a)
            b=intersect([a{i}],[a{j}]);
            if ~isempty(b)
                %判断变量之间是否有交集，有交集才能有可能连接~ismember([1:dimension],c);
                u=copula_diff( inputArg1{i}(:, ~ismember(a{i},b)) , inputArg1{i}(:, ismember(a{i},b)), str{i} , theta{i}  );
                v=copula_diff( inputArg1{j}(:, ~ismember(a{j},b)) , inputArg1{j}(:, ismember(a{j},b)), str{j} , theta{j}  );
                outputArg1{2}{i,j}=[u,v];             
                empirical_T(i,j)=corr( u ,......
                v ,'Type','Kendall');
            end
        end
    end
    outputArg1{1}= empirical_T;
end
end

