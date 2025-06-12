function [DataDivision,numFeatures,numClasses]=DataLstm(inputs,labels,DailyYield,CLOSE,DATE,horizon,DataPercentage)
%% تنظیم پارامترهای تقسیم داده های آموزش ، تست و اعتبارسنجی

DailyYield=DailyYield';
CLOSE=CLOSE';
DATE=DATE';

% افق زمانی
% horizon=1;

%% 
inputs=inputs(:,1:end-horizon);
labels=labels(:,1+horizon:end);
DailyYield=DailyYield(1,1+horizon:end);
CLOSE=CLOSE(1,1+horizon:end);
DATE=DATE(1,1+horizon:end);

x100=1;

nn=size(inputs,2);
%% داده های ورودی

n_train=round(nn*DataPercentage(1)); 
n_Val=round(nn*sum(DataPercentage(1:2)));
Xtrain_Matrix=inputs(:,1:n_train);
XValidation_Matrix=inputs(:,n_train+1:n_Val);
Xtest_Matrix=inputs(:,n_Val+1:end);
x100=1;
% نرمالیزه کردن
% Min_Xtrain_Matrix=min(Xtrain_Matrix,[],2);
% Max_Xtrain_Matrix=max(Xtrain_Matrix,[],2);
% Range_Xtrain_Matrix=Max_Xtrain_Matrix-Min_Xtrain_Matrix;
% Range_Xtrain_Matrix(Range_Xtrain_Matrix==0)=1;
% 
% Xtrain_Matrix_Normalization=(Xtrain_Matrix-Min_Xtrain_Matrix)./Range_Xtrain_Matrix;
% XValidation_Matrix_Normalization=(XValidation_Matrix-Min_Xtrain_Matrix)./Range_Xtrain_Matrix;
% Xtest_Matrix_Normalization=(Xtest_Matrix-Min_Xtrain_Matrix)./Range_Xtrain_Matrix;


[Xtrain_Matrix_Normalization,NormalizedParameters] = mapminmax(Xtrain_Matrix);
XValidation_Matrix_Normalization= mapminmax('apply',XValidation_Matrix,NormalizedParameters);
Xtest_Matrix_Normalization= mapminmax('apply',Xtest_Matrix,NormalizedParameters);

%% ایجاد برچسب  شبکه
% تبدیل عددی به رشته ایی
t=cell(0,1);
t=string(t);
for i=1:length(labels)
    if labels(1,i)==1
        t(1,i)='buy';
    else
        t(1,i)='sell';
    end
end

% دسته بندی
t=categorical(t);
% tt=t(:,1:Tr);
Ytrain_Matrix=t(:,1:n_train);
YValidation_Matrix=t(:,n_train+1:n_Val);
Ytest_Matrix=t(:,n_Val+1:end);

Returns_Ytrain=DailyYield(1,1:n_train);
Returns_YValidation=DailyYield(1,n_train+1:n_Val);
Returns_Ytest=DailyYield(1,n_Val+1:end);


CLOSE_Ytrain=CLOSE(1,1:n_train);
CLOSE_YValidation=CLOSE(1,n_train+1:n_Val);
CLOSE_Ytest=CLOSE(1,n_Val+1:end);

DATE_Ytrain=DATE(1,1:n_train);
DATE_YValidation=DATE(1,n_train+1:n_Val);
DATE_Ytest=DATE(1,n_Val+1:end);

x100=1;
% جدا کردن داده برای درخت تصافی و تست درخت تصمیم تابع هدف الگوریتم ژنتیک

% تعداد طبقات 
classes = categories(t);

%% تعریف ساختار شبکه
numFeatures = height(Xtrain_Matrix);
numClasses =height(classes);

DataDivision=struct('Xtrain_Matrix_Normalization',Xtrain_Matrix_Normalization,...
    'XValidation_Matrix_Normalization',XValidation_Matrix_Normalization,...
    'Xtest_Matrix_Normalization',Xtest_Matrix_Normalization,...
    'Ytrain_Matrix',Ytrain_Matrix,...
    'YValidation_Matrix',YValidation_Matrix,...
    'Ytest_Matrix',Ytest_Matrix,...
    'Returns_Ytrain',Returns_Ytrain,...
    'Returns_YValidation',Returns_YValidation,...
    'Returns_Ytest',Returns_Ytest,...
    'CLOSE_Ytrain',CLOSE_Ytrain,...
    'CLOSE_YValidation',CLOSE_YValidation,...
    'CLOSE_Ytest',CLOSE_Ytest,...
    'DATE_Ytrain',DATE_Ytrain,...
    'DATE_YValidation',DATE_YValidation,...
    'DATE_Ytest',DATE_Ytest,...
    'NormalizedParameters',NormalizedParameters);
    % 'Min_Normalize',Min_Xtrain_Matrix,...
    % 'Max_Normalize',Max_Xtrain_Matrix,...
    % 'Range_Normalize',Range_Xtrain_Matrix);

x1000=1;
end