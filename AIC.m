function [outputArg1] = AIC(inputArg1,inputArg2)
%AIC 此处显示有关此函数的摘要
%   此处显示详细说明
family={'Gaussian','Clayton','Frank','Gumbel','t'};
family_aicbic=zeros(2,length(family));
for i=1:length(family)-1
rhohat = copulafit(family{i},inputArg1);
log_L= sum( log( copulapdf(family{i},inputArg1,rhohat(1,end)) ) );
[aic,bic]= aicbic(log_L,1,size(inputArg1,1));%aicbic输入的分别是似然函数，参数数量，数据规模
family_aicbic(1,i)=aic;
family_aicbic(2,i)=bic;
end
%%
%t-copula单独列出来
[rhohat,nu] = copulafit('t',inputArg1);
log_L= sum( log( copulapdf('t',inputArg1,rhohat(1,end),nu) ) );
[aic,bic]= aicbic(log_L,2,size(inputArg1,1));
family_aicbic(:,end)=[aic,bic]';
min_aicbic=min(family_aicbic');
%%
%判断最小的aic或者bic是那一个类型的copula
if strcmp(inputArg2,'aic')%判断是否等于
    [a b]=find(family_aicbic(1,:)==min_aicbic(1))
    min_aicbic=min_aicbic(1);
elseif strcmp(inputArg2,'bic')%判断是否等于
    [a b]=find(family_aicbic(2,:)==min_aicbic(2))
    min_aicbic=min_aicbic(2);
end
%%
if ~strcmp(family{b}, 't')%判断是否等于t
    rhohat=copulafit(family{b},inputArg1);
    y = copulacdf(family{b},inputArg1,rhohat(1,end));
    outputArg1 ={family{b},rhohat(1,end),min_aicbic,y};
else
    [rhohat,nu]=copulafit(family{b},inputArg1);
    y = copulacdf('t',inputArg1,rhohat(1,end),nu);
%   y = copulacdf('t',inputArg1,copulafit('t',inputArg1),rhohat(1,end),nu);
    outputArg1 ={family{b},[rhohat(1,end),nu],min_aicbic,y};
end
end

