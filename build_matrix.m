 function [outputArg1,outputArg2,outputArg3] = build_matrix(information_matrix)
% clc,clear;
% load information_matrix;
outputArg1=zeros( size(information_matrix,2)+1 , size(information_matrix,2)+1 );%初始化vine矩阵
%存储变量索引
outputArg2=cell( size(information_matrix,2)+1 , size(information_matrix,2)+1 );%初始化vine矩阵
%存储copula
outputArg3=cell( size(information_matrix,2)+1 , size(information_matrix,2)+1 );%初始化vine矩阵
%存储系数
inputArg1={};
for i=1:size(information_matrix,2)
    inputArg1(end+1)={information_matrix{i}(1,:)};
end
%%
%%
%获取并集
union_set=inputArg1;%初始化
for i=2:size(information_matrix,2)
    for j=1:size(inputArg1{i},2)
        union_set{i}{j}=union( union_set{i-1}{ inputArg1{i}{j}(1) },......
        union_set{i-1}{ inputArg1{i}{j}(2) } );
    end
end
%获取差集
diff_set=inputArg1;
for i=2:size(information_matrix,2)
    for j=1:size(inputArg1{i},2)
        diff_set{i}{j}=[ setdiff( union_set{i-1}{ inputArg1{i}{j}(1) },......
            union_set{i-1}{ inputArg1{i}{j}(2) } ) ......
            setdiff( union_set{i-1}{ inputArg1{i}{j}(2) },......
            union_set{i-1}{ inputArg1{i}{j}(1) } ) ];
    end
end
diff_set_copy=diff_set;
max_variance=[];%用来存储边际变量
for w=1:size(information_matrix,2)-1%迭代一次删除边际的变量
    %去max只是为了方便观察变量名但不能保证序号依次向下
    %实际上可以稍微混乱一点，后续进行max操作
    max_variance(end+1)=max(diff_set{end}{1});
    %最后一层的差集中最大的数
    outputArg1(w,w)=max_variance(end);
    for i=1:size(diff_set,2)
        for j=1:size(diff_set{i},2)
            if any( diff_set{i}{j}==outputArg1(w,w) )
                outputArg1(end-i+1,w)=setdiff(diff_set{i}{j},outputArg1(w,w));
                for e=1:size(diff_set_copy{i},2)
                    if isequal(diff_set_copy{i}{e},diff_set{i}{j})
                        outputArg2(end-i+1,w)=information_matrix{i}(2,e);
                        outputArg3(end-i+1,w)=information_matrix{i}(3,e);
                    end
                end
                diff_set{i}(j)=[];
                break;
            end
        end
    end
    diff_set(end)=[];
end
for e=1:size(diff_set_copy{1},2)
    if isequal(diff_set_copy{1}{e},diff_set{1}{1})
        outputArg2(end,end-1)=information_matrix{1}(2,e);
        outputArg3(end,end-1)=information_matrix{1}(3,e);
    end
end
max_variance(end+1)= max( diff_set{end}{1} );
max_variance(end+1)= min( diff_set{end}{1} );
outputArg1(end-1,end-1)=max_variance(end-1);
outputArg1(end,end-1:end)=max_variance(end);
%     index=inputArg1{size(inputArg1,2)}{1};
%     for i=size(inputArg1,2)-1:-1:1%从最后一层的索引开始往上找，
%         %找到边际的索引，继续往上找，直到找到第一层的索引，即为边际的变量
%         index_1=inputArg1{i}{ 1,index(1) };
%         index_2=inputArg1{i}{ 1,index(2) };
%         A=cell2mat(inputArg1{i}(1,:));
%         A=[sum( A == index_1(1)) sum( A == index_1(2))......
%             ;sum(A == index_2(1)) sum( A == index_2(2) )];
%         index(1)=index_1(A(1,:)==1);
%         index(2)=index_2(A(2,:)==1);
%     end
%     outputArg1(w,w)=index(2);
% %     outputArg1(w+1,w)=index(1);
%     max_variance(end+1)=index(2);
%     del_num=max_variance(end);
%     for i=1:size(inputArg1,2)
%         for j=1:size(inputArg1{i},2)
%             if any( inputArg1{1}{j}==del_num )
%                 inputArg1{i}(j)=[];
%                 del_num=j;
%                 break;
%             end
%         end
%     end
%     inputArg1(end)=[];



