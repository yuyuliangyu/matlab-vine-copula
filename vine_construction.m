function [index_matrix,copulafamily_matrix,parameters_matrix,vine_returns]=vine_construction(vine_returns)
dimension=size(vine_returns,2);
%%
transform_variance=vine_returns;
information_matrix=cell(1,dimension-1);
for d=dimension:-1:2
    if d==dimension
        empirical_T=Empirical(transform_variance);%计算经验肯德尔tao
        vine_matrix=MST(abs(empirical_T));    %找出最小生成树
        %这里取abs是因为要找到相关性绝对值之和最大的路径
        %copula拟合
        [a b]=find(vine_matrix);%找到变量之间的连接关系
        transform_matrix=[];
%         c=[];
        for i=1:(d-1)
            transform_matrix=[transform_matrix ......
            [[a(i),b(i)] AIC( [ transform_variance(:,a(i)) , transform_variance(:,b(i)) ],'aic')......
            [ transform_variance(:,a(i)) , transform_variance(:,b(i)) ] ]' ];......
%             w=AIC( [ transform_variance(:,a(i)) , transform_variance(:,b(i)) ],'aic');
%             disp(w);
            %计算每对的copula的aicbic取最优的
%             c=[c transform_matrix{end,i}];
        end     
        information_matrix{1,dimension+1-d}=transform_matrix;
%         transform_variance=copula_diff(transform_variance,information_matrix{1,dimension+1-d});%求偏导数
    else
        partial_derivative=Empirical(transform_matrix);%计算经验肯德尔tao
        empirical_T=partial_derivative{1};
        partial_derivative=partial_derivative{2};
        vine_matrix=MST(abs(empirical_T));    %找出最小生成树
        [a b]=find(vine_matrix);%找到变量之间的连接关系
        transform_matrix=[];
%         c=[];
        for i=1:(d-1)
            %             inform=information_matrix{dimension-d};
            transform_matrix=[transform_matrix ......
            [[a(i),b(i)] AIC( partial_derivative{ a(i) , b(i) } ,'aic') ......
            partial_derivative( a(i) , b(i) ) ]'];
%             transform_matrix=[transform_matrix ......
%             [[a(i),b(i)] AIC( [ transform_variance(:,a(i)) , transform_variance(:,b(i)) ],'aic') ......
%             partial_derivative( a(i) , b(i) ) ]'];
%           [ information_matrix{dimension-d}{5,a(i)} information_matrix{dimension-d}{5,b(i)} ]
            %计算每对的copula的aicbic取最优的
%             c=[c transform_matrix{end,i}];
        end
%         transform_variance=partial_derivative;
        
%         transform_variance=copula_diff(transform_variance,transform_matrix);%求偏导数
        information_matrix{1,dimension+1-d}=transform_matrix;
%         transform_variance=copula_diff(transform_variance,transform_matrix);%求偏导数
    end
    %%
    disp(empirical_T)
end
save information_matrix;
[index_matrix,copulafamily_matrix,parameters_matrix] = build_matrix(information_matrix);
[index_matrix,origin_index] = transform_index_matrix(index_matrix);
vine_returns=vine_returns(:,fliplr(origin_index));
end