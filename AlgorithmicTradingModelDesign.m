function [Output,CrossValidationData]=AlgorithmicTradingModelDesign()

clear
clc
tic

disp('کد نماد مورد نظر که قصد دارید با آن نام نتایج نهایی مدل ذخیره شود را وارد نمائید');
Symbol_Code=input("code stock for example 'IRO1KSIM0002' : " ,"s");
%% بارگذاری داده های قیمتی سهام
% انتخاب فابل اکسل برای بارگذاری داده 

location=uigetdir('C:\', 'مسیر پوشه حاوی فایل اکسل را انتخاب کنید');
% دریافت نام فایل از کاربر (با یا بدون پسوند)
fileName = input('نام فایل اکسل (با یا بدون پسوند) را وارد کنید: ', 's');

% لیست پسوندهای احتمالی برای بررسی (به ترتیب اولویت)
extensions = {'.xlsx', '.xls', '.csv'};
% بررسی آیا نام فایل از قبل پسوند دارد یا خیر
[~, ~, ext] = fileparts(fileName);
if isempty(ext)
    fileFound = false;
    for i=1:length(extensions)
        fullFileName=fullfile(location,[fileName,extensions{i}]);
        if exist(fullFileName,"file")==2
            filePath = fullFileName;
            fileFound = true;
            break
        end
    end
    if ~fileFound
        error('فایل با نام "%s" و پسوندهای %s یافت نشد!', fileName, strjoin(extensions, ', '));
    end
else
    % اگر پسوند دارد، مستقیماً مسیر را بساز
    filePath = fullfile(location, fileName);
    if ~exist(filePath, 'file')
        error('فایل "%s" در مسیر انتخابی وجود ندارد!', fileName);
    end
end

disp(['فایل یافت شد: ', filePath]);
% ادامه پردازش (مثلاً خواندن فایل)
if contains(filePath, '.xlsx') || contains(filePath, '.xls')
    data = readtable(filePath,"VariableNamingRule","preserve");
elseif contains(filePath, '.csv')
    data = readtable(filePath, 'Delimiter', ',','VariableNamingRule','preserve');
else
    error('فرمت فایل پشتیبانی نمیشود!');
end

OPEN=data.("<OPEN>");
HIGH=data.("<HIGH>");
LOW=data.("<LOW>");
CLOSE=data.("<CLOSE>");
VOL=data.("<VOL>");
DATE_Shamsi=data.("<DTYYYYMMDD>");
if ~isnumeric(OPEN)
    OPEN=str2double(OPEN);
    HIGH=str2double(HIGH);
    LOW=str2double(LOW);
    CLOSE=str2double(CLOSE);
    VOL=str2double(VOL);
    DATE=str2double(DATE_Shamsi);
end

DailyYield=[0;price2ret(CLOSE)];
% افزایش یا کاهش بازدهی
% DecreaseOrIncreaseInReturns=1*(DailyYield>0);

%% فاز دوم
% استخراج ویژگی ها و مهندسی ویژگی 
% محاسبه اندیکاتورها و استراتژی معاملاتی سنتی
[IndicatorsData,TraditionalStrategy,GenFeature]=technicalIndicatorsForAlgoTrading(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
StockPriceData=table(CLOSE,HIGH,LOW,OPEN,VOL);
IndicatorsData=[StockPriceData,IndicatorsData];
x100=1;
% مانا کردن داده اندیکاتورها مرتبه اول
stationary1_IndicatorsData=diff(IndicatorsData);
stationary1_IndicatorsData=table2array(stationary1_IndicatorsData);
stationary1_IndicatorsData=[zeros(1,size(stationary1_IndicatorsData,2));stationary1_IndicatorsData];

var1=IndicatorsData.Properties.VariableNames;
for i=1:size(var1,2)
    var2{1,i}=strcat('diff1_',var1{1,i});
end

stationary1_IndicatorsData=array2table(stationary1_IndicatorsData,'VariableNames',var2);

Initial_input_Table=[IndicatorsData,stationary1_IndicatorsData,TraditionalStrategy];
Initial_input_VariableNames=Initial_input_Table.Properties.VariableNames;

Initial_input_DataArray=table2array(Initial_input_Table);
RowSum=sum(Initial_input_DataArray,2);

% مشخص کردن سطرهایی که دارای nan هستنذ
nan_row=any(isnan(RowSum),2);
x100=1;
Initial_input_Table=Initial_input_Table(~nan_row,:);
% حذف سطرهای دارای nan
OPEN=OPEN(~nan_row,:);
HIGH=HIGH(~nan_row,:);
LOW=LOW(~nan_row,:);
CLOSE=CLOSE(~nan_row,:);
VOL=VOL(~nan_row,:);
DATE=DATE(~nan_row,:);
DailyYield=DailyYield(~nan_row,:);
Initial_input_DataArray=Initial_input_DataArray(~nan_row,:);

%%%%%
%% محاسبات مربوط به برچسب ها
x100=1;
% تبدیل به ماتریس استراتژی های سنتی
TraditionalStrategy_DataArray=table2array(TraditionalStrategy);
TraditionalStrategy_DataArray=TraditionalStrategy_DataArray(~nan_row,:);
% تبدیل داده بصورتی که سطرها نشانده ویژگی و ستون نمونه باشد
% TraditionalStrategy_DataArray=TraditionalStrategy_DataArray';
% برچسب زنی داده ها که بعنوان هدف سیستم می باشد
[x,fval,pop]=GaBineryForLabel(TraditionalStrategy_DataArray',CLOSE,DATE);

x=x';
labels=[x ~x];
% سطرها بعنوان کلاس ها و ستونها بعنوان نمونه ها
labels=labels';  % Row matrix

x100=1;
%%  جداسازی داده آموزش ، اعتبارسنجی و تست


DataPercentage=[0.75 0.60 0.45  0.30
    0.15 0.15 0.15 0.15
    0.10 0.25 0.40 0.55];

DataPercentage=[0.75 0.60
    0.15 0.15
    0.10 0.25];

CrossValidationData=cell(size(DataPercentage,2),1);
% افق زمانی پیش بینی
horizon=1;
CrossValidationModel="yes";
%% اعتبارسنجی متقاطع
for fold=1:size(DataPercentage,2)
    [CrossValidationData1,FirstSelectedFeatures,Transformers,SecondSelectedFeatures]=Cross_validation(Initial_input_DataArray',labels,DailyYield,CLOSE,DATE,horizon,Initial_input_VariableNames,CrossValidationModel,DataPercentage(:,fold));
    CrossValidationData{fold,1}=CrossValidationData1;

end
x100=1;
%% ....................................................................................................................................................
%% بهینه سازی پارامترهای شبکه عصبی و انتخاب بهترین ساختار شبکه
ForTheModel="BayesianModel";
results=BayesianOptimizationNeuralNetwork(CrossValidationData,ForTheModel);
% پارامترهای شبکه
% numWindows=results.XAtMinObjective.numWindows;
numLayers=results.XAtMinObjective.numLayers;
numHiddenUnits1=results.XAtMinObjective.numHiddenUnits1;
numHiddenUnits2=results.XAtMinObjective.numHiddenUnits2;
numHiddenUnits3=results.XAtMinObjective.numHiddenUnits3;
numHiddenUnits4=results.XAtMinObjective.numHiddenUnits4;
numHiddenUnits5=results.XAtMinObjective.numHiddenUnits5;

probability1=results.XAtMinObjective.probability1;
probability2=results.XAtMinObjective.probability2;
probability3=results.XAtMinObjective.probability3;
probability4=results.XAtMinObjective.probability4;
probability5=results.XAtMinObjective.probability5;
probability6=results.XAtMinObjective.probability6;

numHiddenUnitsFC1=results.XAtMinObjective.numHiddenUnitsFC1;
numEpochs=results.XAtMinObjective.numEpochs;
MiniBatchSize=results.XAtMinObjective.MiniBatchSize;
gradientThreshold=results.XAtMinObjective.gradientThreshold;
L2Value=results.XAtMinObjective.L2Value;
LearnRate=results.XAtMinObjective.LearnRate;


x100=1;

%% ....................................................................................................................................................
%% مرحله طراحی نهایی مدل
DataPercentage=[0.90 % درصد داده آموزش
    0  % درصد داده اعتبارسنجی
    0.10  % درصد داده تست
    ];

CrossValidationData=cell(size(DataPercentage,2),1);
CrossValidationModel="no";

[CrossValidationData2,FirstSelectedFeatures,Transformers,SecondSelectedFeatures]=Cross_validation(Initial_input_DataArray',labels,DailyYield,CLOSE,DATE,horizon,Initial_input_VariableNames,CrossValidationModel,DataPercentage);
CrossValidationData{1,1}=CrossValidationData2;
%%
%%%%%%%%

x100=1;

%% طراحی مدل نهایی و ارزیابی روی داده تست
ForTheModel="FinalModel";

x100=1;
[ACC,NetBest]=TrainNeuralNetwork(CrossValidationData,ForTheModel,numLayers,numHiddenUnits1,numHiddenUnits2,numHiddenUnits3,numHiddenUnits4,numHiddenUnits5,...
    probability1,probability2,probability3,probability4,probability5,probability6,numHiddenUnitsFC1,numEpochs,MiniBatchSize,gradientThreshold,L2Value,LearnRate);

x100=1;
% مدل نهای شبکه عصبی
Neural_Network_Model=NetBest.net;
% نتایج داده تست
Test_Data_Results=struct('Annual_Sharpe_Ratio_Test_Data',NetBest.SharpeRatio_YPredTest,...)
    'Annual_Return_On_Test_Data',NetBest.Annual_Return_YPredTest,...
    'max_dd_Test_Data',NetBest.max_dd_YPredTest,...
    'SharpeRatio_Complete_Test_Data',NetBest.SharpeRatio_Complete_YPredTest,...
    'Annual_Return_Complete_Test_Data',NetBest.Annual_Return_Complete_YPredTest,...
    'YPred_test',NetBest.YPred_test,...
    'SuccessRate',NetBest.SuccessRate,...
    'Profit_LossRatio',NetBest.Profit_LossRatio,...
    'NumberOfNegativeTrades',NetBest.NumberOfNegativeTrades,...
    'Kelly',NetBest.Kelly,...
    'NumberTrades',NetBest.NumberTrades,...
    'transaction_Table',NetBest.transaction_Table,...
    'Date_Of_Start_Of_Test_Data',NetBest.Date_Of_Start_Of_Test_Data,...
    'End_Date_Of_Test_Data',NetBest.End_Date_Of_Test_Data);

% داده ها برای ذخیره سازی
Output=struct('First_Choice_Features',FirstSelectedFeatures,...
    'Second_Choice_Features',SecondSelectedFeatures,...
    'Feature_Engineering',Transformers,...
    'Neural_Network_Model',Neural_Network_Model,...
    'Min_Normalize',CrossValidationData{1,1}.Min_Normalize,...
    'Max_Normalize',CrossValidationData{1,1}.Max_Normalize,...
    'Range_Normalize',CrossValidationData{1,1}.Range_Normalize,...
    'Test_Data_Results',Test_Data_Results);

% کد ذخیره نتایج را وارد کنید برای مثال IRO1FKHZ0003



FolderPath='D:\StockNeuralNetworkModels';

g=exist(fullfile(FolderPath, [Symbol_Code '.mat']));

if g==0
    save(fullfile(FolderPath, [Symbol_Code '.mat']), 'Output');
else
    delete(fullfile(FolderPath, [Symbol_Code '.mat']));
    save(fullfile(FolderPath, [Symbol_Code '.mat']), 'Output');
end

toc
x100=1;

end








function [CrossValidationData,FirstSelectedFeatures,Transformers,SecondSelectedFeatures]=Cross_validation(Initial_input_DataArray,labels,DailyYield,CLOSE,DATE,horizon,Initial_input_VariableNames,CrossValidationModel,DataPercentage)

[DataDivision,numFeatures,numClasses]=DataLstm(Initial_input_DataArray,labels,DailyYield,CLOSE,DATE,horizon,DataPercentage);

% داده های آموزش ، اعتبارسنجی و تست بصورت سطرها نمونه ها و ستون ها ویژگی ها
Xtrain_VariableNames=Initial_input_VariableNames;
Xtrain_DataArray=DataDivision.Xtrain_Matrix_Normalization';
switch CrossValidationModel
    case "yes"
        XValidation_DataArray=DataDivision.XValidation_Matrix_Normalization';
    case "no"
        XValidation_DataArray=[];
end
Xtest_DataArray=DataDivision.Xtest_Matrix_Normalization';
Ytrain_DataArray=DataDivision.Ytrain_Matrix';
switch CrossValidationModel
    case "yes"
        YValidation_DataArray=DataDivision.YValidation_Matrix';
    case "no"
        YValidation_DataArray=[];
end
Ytest_DataArray=DataDivision.Ytest_Matrix';

Returns_Ytrain=DataDivision.Returns_Ytrain;
switch CrossValidationModel
    case "yes"
        Returns_YValidation=DataDivision.Returns_YValidation;
    case "no"
        Returns_YValidation=[];
end
Returns_Ytest=DataDivision.Returns_Ytest;

CLOSE_Ytrain=DataDivision.CLOSE_Ytrain;
switch CrossValidationModel
    case "yes"
        CLOSE_YValidation=DataDivision.CLOSE_YValidation;
    case "no"
        CLOSE_YValidation=[];
end
CLOSE_Ytest=DataDivision.CLOSE_Ytest;

DATE_Ytrain=DataDivision.DATE_Ytrain;
switch CrossValidationModel
    case "yes"
        DATE_YValidation=DataDivision.DATE_YValidation;
    case "no"
        DATE_YValidation=[];
end
DATE_Ytest=DataDivision.DATE_Ytest;

Min_Normalize=DataDivision.Min_Normalize;
Max_Normalize=DataDivision.Max_Normalize;
Range_Normalize=DataDivision.Range_Normalize;

x100=1;
%% انتخاب ویژگی اولیه
% تبدیل داده آموزش بصورت جدول
Xtrain_Table=array2table(Xtrain_DataArray,"VariableNames",Xtrain_VariableNames);

Ytrain_Table=array2table(Ytrain_DataArray,"VariableNames",{'Label'});
% Label_Returns_Ytrain=array2table(Returns_Ytrain',"VariableNames",{'Label'});

[idx,scores]=fscmrmr(Xtrain_Table,Ytrain_Table);

% ویژگی های انتخابی
SelectedPrimaryFeature=idx(1:150);
% مشخص کردن ویژگی های انتخابی
FirstSelectedFeatures=zeros(1,size(Xtrain_DataArray,2));
FirstSelectedFeatures(SelectedPrimaryFeature)=1;  % باید ذخیره شود

% ویژگی های انتخاب شده در مرحله اول
Xtrain_VariableNames=Xtrain_VariableNames(1,FirstSelectedFeatures==1);
Xtrain_DataArray=Xtrain_DataArray(:,FirstSelectedFeatures==1);
switch CrossValidationModel
    case "yes"
        XValidation_DataArray=XValidation_DataArray(:,FirstSelectedFeatures==1);
    case "no"
        XValidation_DataArray=[];
end
Xtest_DataArray=Xtest_DataArray(:,FirstSelectedFeatures==1);

x100=1;
%% محاسبه تبدیل موجک
% % تبدیل موجک روی داده های آموزش
% Xtrain_DataArray_DWT=[];
% VariableNames_DWT=cell(0,0);
% for ij=1:size(Xtrain_DataArray,2)
%     DWT=[];
%     DWT=modwt(Xtrain_DataArray(:,ij),'db4',4);
%     DWT=DWT';
%     SI_DWT=size(DWT,2);
%     Name_DWT=arrayfun(@(x) sprintf('%s_DWT_%d',Xtrain_VariableNames{1,ij},x),1:SI_DWT,'UniformOutput',false);
%     VariableNames_DWT=[VariableNames_DWT,Name_DWT];
%
%     Xtrain_DataArray_DWT=[Xtrain_DataArray_DWT,DWT];
% end
%
% % تبدیل موجک روی داده های اعتبارسنجی
% XValidation_DataArray_DWT=[];
% for ij=1:size(XValidation_DataArray,2)
%     DWT=[];
%     DWT=modwt(XValidation_DataArray(:,ij),'db4',4);
%     DWT=DWT';
%     XValidation_DataArray_DWT=[XValidation_DataArray_DWT,DWT];
% end
%
% % تبدیل موجک روی داده های تست
% Xtest_DataArray_DWT=[];
% for ij=1:size(Xtest_DataArray,2)
%     DWT=[];
%     DWT=modwt(Xtest_DataArray(:,ij),'db4',4);
%     DWT=DWT';
%     Xtest_DataArray_DWT=[Xtest_DataArray_DWT,DWT];
% end
%
% % اضافه کردن ویژگی جدید با تبدیل موجک به داده ها
% Xtrain_VariableNames=[Xtrain_VariableNames,VariableNames_DWT];
% Xtrain_DataArray=[Xtrain_DataArray,Xtrain_DataArray_DWT];
% XValidation_DataArray=[XValidation_DataArray,XValidation_DataArray_DWT];
% Xtest_DataArray=[Xtest_DataArray,Xtest_DataArray_DWT];
%
x100=1;

%% انتخاب ویژگی دوم
% % تبدیل داده آموزش بصورت جدول
% Xtrain_Table=array2table(Xtrain_DataArray,"VariableNames",Xtrain_VariableNames);
%
% Ytrain_Table=array2table(Ytrain_DataArray,"VariableNames",{'Label'});
% % Label_Returns_Ytrain=array2table(Returns_Ytrain',"VariableNames",{'Label'});
%
% [idx,scores]=fscmrmr(Xtrain_Table,Ytrain_Table);
% x100=1;
% % ویژگی های انتخابی
% SelectedPrimaryFeature=idx(1:150);
% % مشخص کردن ویژگی های انتخابی
% FirstSelectedFeatures=zeros(1,size(Xtrain_DataArray,2));
% FirstSelectedFeatures(SelectedPrimaryFeature)=1;  % باید ذخیره شود
%
% % ویژگی های انتخاب شده در مرحله اول
% Xtrain_VariableNames=Xtrain_VariableNames(1,FirstSelectedFeatures==1);
% Xtrain_DataArray=Xtrain_DataArray(:,FirstSelectedFeatures==1);
% XValidation_DataArray=XValidation_DataArray(:,FirstSelectedFeatures==1);
% Xtest_DataArray=Xtest_DataArray(:,FirstSelectedFeatures==1);


%% مهندسی ویژگی

Xtrain_Table=array2table(Xtrain_DataArray,"VariableNames",Xtrain_VariableNames);
switch CrossValidationModel
    case "yes"
        XValidation_Table=array2table(XValidation_DataArray,"VariableNames",Xtrain_VariableNames);
    case "no"
        XValidation_Table=[];
end
Xtest_Table=array2table(Xtest_DataArray,"VariableNames",Xtrain_VariableNames);

Xtrain_GenFeature=[Xtrain_Table,Ytrain_Table];
[Transformers,NewTbl]=gencfeatures(Xtrain_GenFeature,"Label",5000,"TransformedDataStandardization","none","CategoricalEncodingLimit",1);

Xtrain_Table=transform(Transformers,Xtrain_Table);
switch CrossValidationModel
    case "yes"
        XValidation_Table=transform(Transformers,XValidation_Table);
    case "no"
        XValidation_Table=[];
end
Xtest_Table=transform(Transformers,Xtest_Table);

% داده های جدید
Xtrain_VariableNames=Xtrain_Table.Properties.VariableNames;
Xtrain_DataArray=table2array(Xtrain_Table);

switch CrossValidationModel
    case "yes"
        XValidation_DataArray=table2array(XValidation_Table);
    case "no"
        XValidation_DataArray=[];
end
Xtest_DataArray=table2array(Xtest_Table);


%% انتخاب ویژگی دوم
x100=1;
% [idx,scores]= fsrmrmr(Xtrain_Matrix_GenFeature_mrmr,Ytrain_Matrix_mrmr);
[idx,scores]=fscmrmr(Xtrain_Table,Ytrain_Table);
x100=1;
idx=idx(1:50);

SecondSelectedFeatures=zeros(1,size(Xtrain_DataArray,2));
SecondSelectedFeatures(idx)=1;  % باید ذخیره شود


Xtrain_VariableNames=Xtrain_VariableNames(1,SecondSelectedFeatures==1);
Xtrain_DataArray=Xtrain_DataArray(:,SecondSelectedFeatures==1);

switch CrossValidationModel
    case "yes"
        XValidation_DataArray=XValidation_DataArray(:,SecondSelectedFeatures==1);
    case "no"
        XValidation_DataArray=[];
end
Xtest_DataArray=Xtest_DataArray(:,SecondSelectedFeatures==1);

x100=1;


%% کاهش بعد داده ها
% کاهش بعد داده
% results=BayesianOptimizationForEncoder(Xtrain_DataArray');
% MaxEpochs=results.XAtMinObjective.MaxEpochs;
% hiddenSize=results.XAtMinObjective.hiddenSize;
% L2Regularization=results.XAtMinObjective.L2Regularization;
% Sparsity_Regularization=results.XAtMinObjective.Sparsity_Regularization;
% Sparsity_Proportion=results.XAtMinObjective.Sparsity_Proportion;
% % اجرای مدل کاهش بعد با پارامترهای تنظیم شده در مرحله قبل
% [mseError_best,autoenc_best]=encoders(Xtrain_DataArray',MaxEpochs,hiddenSize,L2Regularization,Sparsity_Regularization,Sparsity_Proportion);
% % mseError_best=inf;
% x100=1;
% for i=1:20
%     [mseError,autoenc]=encoders(Xtrain_DataArray',MaxEpochs,hiddenSize,L2Regularization,Sparsity_Regularization,Sparsity_Proportion);
%     if mseError<mseError_best
%         mseError_best=mseError
%         autoenc_best=autoenc;
%     end
% end
%
%
%
% Xtrain_DataArray=encode(autoenc_best,Xtrain_DataArray');
% XValidation_DataArray=encode(autoenc_best,XValidation_DataArray');
% Xtest_DataArray=encode(autoenc_best,Xtest_DataArray');


Data_Encoder=struct('Xtrain_Matrix_Encoder',Xtrain_DataArray',...
    'XValidation_Matrix_Encoder',XValidation_DataArray',...
    'Xtest_Matrix_Encoder',Xtest_DataArray',...
    'Ytrain_Matrix',Ytrain_DataArray', ...
    'YValidation_Matrix',YValidation_DataArray',...
    'Ytest_Matrix',Ytest_DataArray',...
    'Returns_Ytrain',Returns_Ytrain,...
    'Returns_YValidation',Returns_YValidation,...
    'Returns_Ytest',Returns_Ytest,...
    'CLOSE_Ytrain',CLOSE_Ytrain,...
    'CLOSE_YValidation',CLOSE_YValidation,...
    'CLOSE_Ytest',CLOSE_Ytest,...
    'DATE_Ytrain',DATE_Ytrain,...
    'DATE_YValidation',DATE_YValidation,...
    'DATE_Ytest',DATE_Ytest,...
    'Min_Normalize',Min_Normalize,...
    'Max_Normalize',Max_Normalize,...
    'Range_Normalize',Range_Normalize);

x100=1;
CrossValidationData=Data_Encoder;











end
