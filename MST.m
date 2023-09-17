function [outputArg1] = MST(empirical_T)
%MST 此处显示有关此函数的摘要
%   此处显示详细说明
dimension=length(empirical_T);
vine_matrix=zeros(dimension,dimension);
for i=1:dimension-1
    [a b]=find(vine_matrix~=0);
    c=unique([a b]);
    index= ~ismember([1:dimension],c);%找到其余点
    if ~isempty(c)
        max_tao=max( max(max(empirical_T(c,index))) , max(max(empirical_T(index,c))) );
        %在其余点和已连接点的边中找到最大值，交换次序是因为只有上三角矩阵存储了信息，
        %不考虑全面的话会遗漏一些相关系数，比如c=[2 3],index为=[1 4 5]；
    else
        max_tao=max( empirical_T(:) );
    end
    [a b]=find(empirical_T==max_tao);
    vine_matrix(a,b)=max_tao;
end
outputArg1 = vine_matrix;
end

