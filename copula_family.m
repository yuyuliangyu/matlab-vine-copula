function [outputArg1] = copula_family(varargin)
%COPULA_FAMILY 此处显示有关此函数的摘要
%   此处显示详细说明
syms u v a nu;
%%
C_gumbel_cdf= exp(-((-log(u)).^a + (-log(v)).^a).^(1./a));
%gumbel的分布函数,如下同理
C_gumbel_pdf=diff(diff(C_gumbel_cdf,v,1),u,1);
%gumbel的密度函数,如下同理
C_gumbel_diff=diff(C_gumbel_cdf,v,1);
%对v求偏导,如下同理
%%
C_clayton_cdf= (u.^(-a) + v.^(-a) - 1).^(-1./a);
C_clayton_pdf=diff(diff(C_clayton_cdf,v,1),u,1);
C_clayton_diff=diff(C_clayton_cdf,v,1);
%%
C_frank_cdf= -(1./a).*log(1 + (exp(-a.*u)-1).*(exp(-a.*v)-1)./(exp(-a)-1));
C_frank_pdf=diff(diff(C_frank_cdf,v,1),u,1);
C_frank_diff=diff(C_frank_cdf,v,1);
%%
C_gaussian_pdf=( (1-a.^2).^( -(1./2) ) ).*( exp( -(1./2).*( u.^2+v.^2-2.*a.*u.*v ).*(1-a.^2) ) ).*( exp( (1./2).*(u.^2+v.^2) ) );
%%
%注意，这里的u.v是逆累计分布，所以在计算t分布时候，需要对u.v进行逆累积分布变换
C_t_pdf=(1)./( sqrt(1-a.^2) ).*( gamma( (nu+2)./2 ).*gamma( nu./2 ) )......
    ./( (gamma( (nu+1)./2) ).^2 ).*( 1+(u.^2+v.^2-2.*a.*u.*v )......
    ./( nu.*(1-a.^2) ) ).^( -(nu+2)./2 )./( (1+(u.^2)./(nu) ).^( -(nu+1)./2 ).*( 1+(v.^2)./(nu) ).^( -(nu+1)./2) );
%%
copula_family=cell(1,3);
copula_family{1}={C_gumbel_cdf,C_gumbel_pdf,C_gumbel_diff};
copula_family{2}={C_clayton_cdf,C_clayton_pdf,C_clayton_diff};
copula_family{3}={C_frank_cdf,C_frank_pdf,C_frank_diff};
%%
copulafamily=varargin{1};
%%
if strcmp(copulafamily,'Gumbel')
    outputArg1=copula_family{1};  
elseif strcmp(copulafamily,'Clayton')
    outputArg1=copula_family{2};
elseif strcmp(copulafamily,'Frank')
    outputArg1=copula_family{3};
elseif strcmp(copulafamily,'Gaussian')
    outputArg1={C_gaussian_pdf};
elseif strcmp(copulafamily,'t')
    outputArg1={C_t_pdf};
end

end

