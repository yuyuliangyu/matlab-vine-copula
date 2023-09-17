function [outputArg1] = copula_diff(inputArg1,inputArg2,varargin)
syms u v a nu;
%COPULA_DIFF  copula(u，v)函数的v的偏导数
% family={'Gaussian','Clayton','Frank','Gumbel','t'};
% information=varargin{1,1};
%%
if size(varargin,2)==2
    u=inputArg1;
    v=inputArg2;
    str=varargin{1};
    a=varargin{2};
    % outputArg1=[inputArg1,inputArg2];
    if strcmp(str,'Gumbel')
        %gumbel  C = exp(-((-log(u))^a + (-log(v))^a)^(1/a));    d=diff(C,v,1);
        outputArg1=(exp(-((-log(u)).^a + (-log(v)).^a).^(1/a)).*(-log(v)).^(a - 1).*((-log(u)).^a + (-log(v)).^a).^(1/a - 1))./v;
    elseif strcmp(str,'Clayton')
        %clayton C = (u^(-a) + v^(-a) - 1)^(-1/a);
        outputArg1=1./(v.^(a + 1).*(1./u.^a + 1./v.^a - 1).^(1/a + 1));
    elseif strcmp(str,'Frank')
        %frank   C = -(1/a)*log(1 + (exp(-a*u)-1)*(exp(-a*v)-1)/(exp(-a)-1));
        outputArg1=(exp(-a.*v).*(exp(-a.*u) - 1))./((exp(-a) - 1).*(((exp(-a.*u) - 1).*(exp(-a.*v) - 1))./(exp(-a) - 1) + 1));
    elseif strcmp(str,'Gaussian')
        outputArg1=(  copulacdf(str,[u v+0.001],a)  -  copulacdf(str,[u v],a)  )./(0.001);
    elseif strcmp(str,'t')
        outputArg1=(  copulacdf(str,[u v+0.001],a(1),a(2))  -  copulacdf(str,[u v],a(1),a(2))  )./(0.001);
    end
    outputArg1=double(outputArg1);
% elseif size(varargin,2)==1
%     outputArg1=1;
elseif length(varargin)==3 && varargin{3}=="inverse"
    v=inputArg2;
    y=inputArg1;
    str=varargin{1};
    a=varargin{2};
    % outputArg1=[inputArg1,inputArg2];
    if strcmp(str,'Gumbel')
        %gumbel  C = exp(-((-log(u))^a + (-log(v))^a)^(1/a));    d=diff(C,v,1);
        outputArg1=(exp(-((-log(u)).^a + (-log(v)).^a).^(1/a)).*(-log(v)).^(a - 1).*((-log(u)).^a + (-log(v)).^a).^(1/a - 1))./v;
        u=vpasolve(y==(exp(-((-log(u)).^a + (-log(v)).^a).^(1/a)).*(-log(v)).^(a - 1).*((-log(u)).^a + (-log(v)).^a).^(1/a - 1))./v,u,[0 1]);
        u=vpasolve(y==u.*u+0.5*v.*v,u,[0 1]);
    elseif strcmp(str,'Clayton')
        %clayton C = (u^(-a) + v^(-a) - 1)^(-1/a);
        outputArg1=1./(v.^(a + 1).*(1./u.^a + 1./v.^a - 1).^(1/a + 1));
    elseif strcmp(str,'Frank')
        %frank   C = -(1/a)*log(1 + (exp(-a*u)-1)*(exp(-a*v)-1)/(exp(-a)-1));
        outputArg1=(exp(-a.*v).*(exp(-a.*u) - 1))./((exp(-a) - 1).*(((exp(-a.*u) - 1).*(exp(-a.*v) - 1))./(exp(-a) - 1) + 1));
    elseif strcmp(str,'Gaussian')
        vpasolve(outputArg1==(  copulacdf(str,[u v+0.001],a)  -  copulacdf(str,[u v],a)  )./(0.001),u,[0 1]);
    elseif strcmp(str,'t')
        outputArg1=(  copulacdf(str,[u v+0.001],a(1),a(2))  -  copulacdf(str,[u v],a(1),a(2))  )./(0.001);
    end
    outputArg1=double(outputArg1);
end
end
% syms x;
% p=[1 4 0;-3 -2 -3;2 0 0]
% pt=transpose(p)
% solve(pt*x*p-x==-eye(3))


