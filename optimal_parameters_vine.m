function [outputArg1,outputArg2] = optimal_parameters_vine(varargin)
%OPTIMAL_PARAMETERS_VINE 此处显示有关此函数的摘要
%   此处显示详细说明
%生成参数变量,赋值到参数矩阵中
index_matrix=varargin{1};
copulafamily_matrix=varargin{2};
parameters_matrix=varargin{3};
variance=varargin{4};
%%
dimension=size(index_matrix,2);
%%
%创建参数变量，并且生成句柄函数
variance_str='';
count=0;
for j = 1:dimension
    for i = j+1:dimension
        count=count+1;
        variance_str=[ variance_str ,',',[ '(p','(',num2str(count) ],'))' ];
        if strcmp(copulafamily_matrix{i,j},'t')
            count=count+1;
            variance_str=[ variance_str ,',',[ '(p','(',num2str(count) ],'))' ];
        end
    end
end
variance_str=[ '[',variance_str(2:end) ,']' ];
variance_str=[ 'ln_lik=@(p)(-lik_vine( index_matrix ,copulafamily_matrix ,'......
    , variance_str , ', variance ));'];
%%
%控制参数范围
A=zeros(3*dimension,count);
B=zeros(3*dimension,1);
lowbound=-1.*inf(1,count);
upbound=inf(1,count);
initial_value=zeros(1,count);
count=0;
for j = 1:dimension
    for i = j+1:dimension
        count=count+1;
        if strcmp(copulafamily_matrix{i,j},'Gumbel')
            A(3*count-2,count)=-1;
            B(3*count-2,1)=-1;
            initial_value(count)=parameters_matrix{i,j};
            lowbound(count)=1;
        elseif strcmp(copulafamily_matrix{i,j},'Clayton')
            A(3*count-2,count)=-1;
            B(3*count-2,1)=-0.001;
            initial_value(count)=parameters_matrix{i,j};
            lowbound(count)=0.001;
        elseif strcmp(copulafamily_matrix{i,j},'Frank')
            A(3*count-2,count)=-1;
            B(3*count-2,1)=0.001;
            A(3*count-1,count)=1;
            B(3*count-1,1)=-0.001;
            initial_value(count)=parameters_matrix{i,j};
        elseif strcmp(copulafamily_matrix{i,j},'Gaussian')
            A(3*count-2,count)=1;
            B(3*count-2,1)=1;
            A(3*count-1,count)=-1;
            B(3*count-1,1)=1;
            initial_value(count)=parameters_matrix{i,j};
            lowbound(count)=-1;
            upbound(count)=1;
        elseif strcmp(copulafamily_matrix{i,j},'t')
            A(3*count-2,count)=1;
            B(3*count-2,1)=1;
            A(3*count-1,count)=-1;
            B(3*count-1,1)=1;
            A(3*count,count+1)=-1;
            B(3*count,1)=-1;
            initial_value(count)=parameters_matrix{i,j}(1);
            initial_value(count+1)=parameters_matrix{i,j}(2);
            lowbound(count)=-1;
            upbound(count)=1;
            count=count+1;
        end
    end
end
A( ~any(A')',: )=[];
B( ~any(B,2)')=[];
%%
%全局最优算法
eval(variance_str);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',initial_value,...
                             'objective',ln_lik,'Aineq',A,'bineq',B,'lb',initial_value-abs(initial_value)*0.02,'ub',initial_value+abs(initial_value)*0.02);
% problem.solver = 'simulannealbnd';
% % [jie,fval] = simulannealbnd(ln_lik,initial_value,lowbound,upbound);
% [jie,fval] = simulannealbnd(ln_lik,initial_value,initial_value-abs(initial_value)*0.02,initial_value+abs(initial_value)*0.02);
% para_maxtrix=zeros(dimension,dimension);
% count=1;
% for j=1:dimension-1
%     for i=j+1:dimension
%         para_maxtrix(i,j)=jie(count);
%         count=count+1;
%     end
% end
% likk=lik_vine(index_matrix,copulafamily_matrix,para_maxtrix,variance);
[xmin,fmin,flag,outpt,allmins] = run(gs,problem);
outputArg1 = xmin;
outputArg2 = -fmin;
end

