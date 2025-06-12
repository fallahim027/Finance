function [ModelOutputResults,BuyingPosition,Forecast]=PredictionWithNeuralNetwork()

clear
clc
tic 

%%  ساختار خروجی نتایج

ModelOutputResults=cell(0,0);

ModelOutputResults{1,2}="کد نماد";
ModelOutputResults{2,2}="نماد";

ModelOutputResults{3,1}="بررسی وضعیت سهام";
ModelOutputResults{4,1}="بررسی وضعیت سهام";

ModelOutputResults{3,2}="موقعیت در روز آتی";
ModelOutputResults{4,2}="تعداد رخداد موقعیت";

ModelOutputResults{6,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{7,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{8,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{9,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{10,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{11,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{12,1}="آماره مدل و مدیریت ریسک";
ModelOutputResults{13,1}="آماره مدل و مدیریت ریسک";


ModelOutputResults{6,2}="نسبت شارپ سالیانه ";
ModelOutputResults{7,2}="بازدهی سالیانه ";
ModelOutputResults{8,2}="بیشینه افت سرمایه ";
ModelOutputResults{9,2}="درصد موفقیت ";
ModelOutputResults{10,2}="نسبت سود به زیان ";
ModelOutputResults{11,2}="تعداد معاملات منفی متوالی";
ModelOutputResults{12,2}="حداکثر سرمایه گذاری با شاخص Kelly";
ModelOutputResults{13,2}="تعداد معاملات";


ModelOutputResults{15,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{16,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{17,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{18,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{19,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{20,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{21,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{22,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{23,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{24,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{25,1}="در زمان آزمایش و طراحی مدل";
ModelOutputResults{26,1}="در زمان آزمایش و طراحی مدل";


ModelOutputResults{15,2}="تاریخ شروع داده آزمایش";
ModelOutputResults{16,2}="تاریخ پایان داده آزمایش";
ModelOutputResults{17,2}="نسبت شارپ سالیانه معاملات";
ModelOutputResults{18,2}="بازدهی سالیانه معاملات ";
ModelOutputResults{19,2}="بیشینه افت سرمایه معاملات";
ModelOutputResults{20,2}="درصد موفقیت معاملات";
ModelOutputResults{21,2}="نسبت سود به زیان معاملات";
ModelOutputResults{22,2}="تعداد معاملات منفی متوالی";
ModelOutputResults{23,2}="حداکثر سرمایه گذاری با شاخص Kelly";
ModelOutputResults{24,2}="تعداد معاملات";
ModelOutputResults{25,2}="نسبت شارپ سالیانه کل داده آزمایش";
ModelOutputResults{26,2}="بازدهی سالیانه کل داده آزمایش";

%% خروجی هایی که در موقعیت خرید و نگهداری هستند
BuyingPosition=cell(0,0);
BuyingPosition{1,2}="کد نماد";
BuyingPosition{2,2}="نماد";

BuyingPosition{3,1}="بررسی وضعیت سهام";
BuyingPosition{4,1}="بررسی وضعیت سهام";

BuyingPosition{3,2}="موقعیت در روز آتی";
BuyingPosition{4,2}="تعداد رخداد موقعیت";

BuyingPosition{6,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{7,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{8,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{9,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{10,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{11,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{12,1}="آماره مدل و مدیریت ریسک";
BuyingPosition{13,1}="آماره مدل و مدیریت ریسک";


BuyingPosition{6,2}="نسبت شارپ سالیانه ";
BuyingPosition{7,2}="بازدهی سالیانه ";
BuyingPosition{8,2}="بیشینه افت سرمایه ";
BuyingPosition{9,2}="درصد موفقیت ";
BuyingPosition{10,2}="نسبت سود به زیان ";
BuyingPosition{11,2}="تعداد معاملات منفی متوالی";
BuyingPosition{12,2}="حداکثر سرمایه گذاری با شاخص Kelly";
BuyingPosition{13,2}="تعداد معاملات";


%% مسیر لود داده

[file,location]=uigetfile('*.*','فایل اکسل نام نمادها را انتخاب کنید');

file=fullfile(location,file);
% فراخوانی داده نام نمادها
NameSymbol_Table=readtable(file,"VariableNamingRule","preserve");

% مسیر خروجی های بک تست مدل
location_Backtest='D:\';
file_Backtest='StockNeuralNetworkModels';
location_Backtest=fullfile(location_Backtest,file_Backtest);
% مسیر داده های اکسل دده های تاریخی سهام
location_DataStock='C:\Users\USER\Documents\TseClient 2.0\Adjusted';

cc='.mat';
ex='.xls';

ColumnNumber=2;
ColumnNumberBuyingPosition=2;

x100=1;
for i=1%1:size(NameSymbol_Table,1)
    CodeSymbol=NameSymbol_Table{i,1};
    CodeSymbol=string(CodeSymbol);

    NameSymbol=NameSymbol_Table{i,2};
    NameSymbol=string(NameSymbol);

    NameExcel=NameSymbol_Table{i,3};
    NameExcel=string(NameExcel);
    % نام فایل اکسل با پسوند
    Path_Excel_with_Extension=strcat(NameExcel,ex);
    CodeSymbol_Backtest=strcat(CodeSymbol,cc);

    % مسیر و نام فایل اکسل
    ExcelDataPath=fullfile(location_DataStock,Path_Excel_with_Extension);
    % فایل داده های بک تست مدل سهام
    File_CodeSymbol_Backtest=fullfile(location_Backtest,CodeSymbol_Backtest);
    if exist(File_CodeSymbol_Backtest)==2 && exist(ExcelDataPath)==2
        ColumnNumber=ColumnNumber+1;
        % بازخوانی پارامترهای شبکه عصبی و نتایج بک تست مدل سهام
        ModelParameters=load(File_CodeSymbol_Backtest);

        % فراخوانی داده تاریخی سهام
        
        Data=readtable(ExcelDataPath,VariableNamingRule="preserve");

        [Forecast,Data_Table]=Backtest(Data,ModelParameters);
        % انتخاب داده های از تست مدل به بعد
        Data_Table1=table2array(Data_Table);
        TestDate=1*(Data_Table.DATE>=ModelParameters.Output.Test_Data_Results.Date_Of_Start_Of_Test_Data);
        Data_Table1=Data_Table1(TestDate==1,:);
        Data_Table_VaribleName=Data_Table.Properties.VariableNames;
        Data_Table1=array2table(Data_Table1,"VariableNames",Data_Table_VaribleName);

        % برچسب های متناظر با داده های به بعد
        Forecast=Forecast';
        SizeDataTable1=size(Data_Table1,1);
        SizeForecast=size(Forecast,1);
        Forecast1=Forecast(SizeForecast-SizeDataTable1:end-1,1);

        % محاسبه ابزارهای مدیریت ریسک
        Returns=(Data_Table1.DailyYield)';
        YPred=Forecast1';
        [SharpeRatio,max_dd,Annual_Return,SharpeRatio_Complete,Annual_Return_Complete]=SharpeRatioForNeuralNetworkEvaluation(YPred,Returns);

        CLOSE_Ytest=(Data_Table1.CLOSE)';
        DATE_Ytest=(Data_Table1.DATE)';
        [RiskManag]=RiskManagement(CLOSE_Ytest,DATE_Ytest,YPred);

        % بررسی وضعیت خرید یا فروش
        Forecast=string(Forecast);
        Nu_Event=zeros(size(Forecast,1),1);
        Nu_Event(1,1)=1;
        for j=2:size(Forecast,1)
            if strcmp(Forecast(j,1),"sell") && strcmp(Forecast(j-1,1),"sell")
                Nu_Event(j,1)=Nu_Event(j-1)+1;
                Status="فعلا خرید نکنید";
            elseif strcmp(Forecast(j,1),"sell") && strcmp(Forecast(j-1,1),"buy")
                Nu_Event(j,1)=1;
                Status="اگر سهام دارید بفروشید";
            elseif strcmp(Forecast(j,1),"buy") && strcmp(Forecast(j-1,1),"buy")
                Nu_Event(j,1)=Nu_Event(j-1)+1;
                Status="خرید خود را نگهدارید";
            elseif strcmp(Forecast(j,1),"buy") && strcmp(Forecast(j-1,1),"sell")
                Nu_Event(j,1)=1;
                Status="اقدام به خرید کنید";
            end
        end
        NumberOfEvents=Nu_Event(end,1);



        % نمایش خروجی ها

        ModelOutputResults{1,ColumnNumber}=CodeSymbol;
        ModelOutputResults{2,ColumnNumber}=NameSymbol;
        ModelOutputResults{3,ColumnNumber}=Status;
        ModelOutputResults{4,ColumnNumber}=NumberOfEvents;
        ModelOutputResults{6,ColumnNumber}=round(SharpeRatio,2);
        ModelOutputResults{7,ColumnNumber}=round(Annual_Return,2);
        ModelOutputResults{8,ColumnNumber}=round(max_dd,2);
        ModelOutputResults{9,ColumnNumber}=round(RiskManag.SuccessRate,2);
        ModelOutputResults{10,ColumnNumber}=round(RiskManag.Profit_LossRatio,2);
        ModelOutputResults{11,ColumnNumber}=round(RiskManag.NumberOfNegativeTrades,2);
        ModelOutputResults{12,ColumnNumber}=round(RiskManag.Kelly,2);
        ModelOutputResults{13,ColumnNumber}=round(RiskManag.NumberTrades,2);

        ModelOutputResults{15,ColumnNumber}=ModelParameters.Output.Test_Data_Results.Date_Of_Start_Of_Test_Data;
        ModelOutputResults{16,ColumnNumber}=ModelParameters.Output.Test_Data_Results.End_Date_Of_Test_Data;
        ModelOutputResults{17,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.Annual_Sharpe_Ratio_Test_Data,2);
        ModelOutputResults{18,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.Annual_Return_On_Test_Data,2);
        ModelOutputResults{19,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.max_dd_Test_Data,2);
        ModelOutputResults{20,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.SuccessRate,2);
        ModelOutputResults{21,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.Profit_LossRatio,2);
        ModelOutputResults{22,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.NumberOfNegativeTrades,2);
        ModelOutputResults{23,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.Kelly,2);
        ModelOutputResults{24,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.NumberTrades,2);
        ModelOutputResults{25,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.SharpeRatio_Complete_Test_Data,2);
        ModelOutputResults{26,ColumnNumber}=round(ModelParameters.Output.Test_Data_Results.Annual_Return_Complete_Test_Data,2);

        if (strcmp(Status,"اقدام به خرید کنید") || strcmp(Status,"خرید خود را نگهدارید") || strcmp(Status,"اگر سهام دارید بفروشید" )) && ModelParameters.Output.Test_Data_Results.Kelly > 0.60 && ModelParameters.Output.Test_Data_Results.Annual_Sharpe_Ratio_Test_Data > 1.5
            ColumnNumberBuyingPosition=ColumnNumberBuyingPosition+1;
            BuyingPosition{1,ColumnNumberBuyingPosition}=CodeSymbol;
            BuyingPosition{2,ColumnNumberBuyingPosition}=NameSymbol;
            BuyingPosition{3,ColumnNumberBuyingPosition}=Status;
            BuyingPosition{4,ColumnNumberBuyingPosition}=NumberOfEvents;
            BuyingPosition{6,ColumnNumberBuyingPosition}=round(SharpeRatio,2);
            BuyingPosition{7,ColumnNumberBuyingPosition}=round(Annual_Return,2);
            BuyingPosition{8,ColumnNumberBuyingPosition}=round(max_dd,2);
            BuyingPosition{9,ColumnNumberBuyingPosition}=round(RiskManag.SuccessRate,2);
            BuyingPosition{10,ColumnNumberBuyingPosition}=round(RiskManag.Profit_LossRatio,2);
            BuyingPosition{11,ColumnNumberBuyingPosition}=round(RiskManag.NumberOfNegativeTrades,2);
            BuyingPosition{12,ColumnNumberBuyingPosition}=round(RiskManag.Kelly,2);
            BuyingPosition{13,ColumnNumberBuyingPosition}=round(RiskManag.NumberTrades,2);



        end


x100=1;










    end
end


% writecell(ModelOutputResults,'C.xls');


toc
end




function [Forecast,Data_Table]=Backtest(Data,ModelParameters)

%% استخراج داده های قیمتی سهام
OPEN=Data.("<OPEN>");
HIGH=Data.("<HIGH>");
LOW=Data.("<LOW>");
CLOSE=Data.("<CLOSE>");
VOL=Data.("<VOL>");
DATE_Shamsi=Data.("<DTYYYYMMDD>");
if ~isnumeric(OPEN)
    OPEN=str2double(OPEN);
    HIGH=str2double(HIGH);
    LOW=str2double(LOW);
    CLOSE=str2double(CLOSE);
    VOL=str2double(VOL);
    DATE=str2double(DATE_Shamsi);
end
DailyYield=[0;price2ret(CLOSE)];


%% استخراج پارامترهای مدل و پیش پردازش مدل
% پارامترهای نرمالیزه کردن داده ها
Max_Normalize=ModelParameters.Output.Max_Normalize';
Min_Normalize=ModelParameters.Output.Min_Normalize';
Range_Normalize=ModelParameters.Output.Range_Normalize';

% انتخاب ویژگی مرحله اول
SelectedPrimaryFeature=ModelParameters.Output.First_Choice_Features;
% مدل مهندسی ویژگی
Transformers=ModelParameters.Output.Feature_Engineering;
% انتخاب ویژگی مرحله دوم
SecondSelectedFeatures=ModelParameters.Output.Second_Choice_Features;
% مدل شبکه عصبی
Net=ModelParameters.Output.Neural_Network_Model;

%%

% فاز اول محاسبه اندیکاتورهای تکنیکال و استراتژی های سنتی
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

Data_Table=table(DATE,OPEN,HIGH,LOW,CLOSE,VOL,DailyYield);

% نرمالیزه کردن داده ها
Initial_input_Table=(Initial_input_Table-Min_Normalize)./Range_Normalize;

% ویژگی های انتخاب شده در مرحله اول

ForecastData=Initial_input_Table(:,SelectedPrimaryFeature==1);

% مهندسی ویژگی خودکار
ForecastData=transform(Transformers,ForecastData);

% ویژگی های انتخاب شده در مرحله دوم
ForecastData=ForecastData(:,SecondSelectedFeatures==1);


% آماده سازی داده ها برای ورود به شبکه
ForecastData_NameVarible=ForecastData.Properties.VariableNames;
ForecastData=table2array(ForecastData);

ForecastData=ForecastData';
x100=1;
% پیش بینی با مدل شبکه عصبی
Forecast=classify(Net,ForecastData);








end






function [SharpeRatio,max_dd,Annual_Return,SharpeRatio_Complete,Annual_Return_Complete]=SharpeRatioForNeuralNetworkEvaluation(YPred,Returns)

YPred_test_double=double(YPred);
YPred_test_double(YPred_test_double==2)=0;

Returns_Ytest_adj=Returns;
% هزینه معاملاتی بعلاوه مقداری زیان غیرعادی
TransactionCost=-0.03;
n_active=length(YPred_test_double);
x100=1;
for i=2:n_active
    if YPred_test_double(1,i-1)==0 && YPred_test_double(1,i)==1
        Returns_Ytest_adj(1,i)=TransactionCost;
    end
end
% تعداد کل روز معاملاتی سهام
n1=numel(Returns_Ytest_adj);
% تعداد سال
n_years=n1/252;

if sum(YPred_test_double)>10
    array=Returns_Ytest_adj(YPred_test_double==1);
    n_active=width(array);
    n_active2=width(Returns);
    % متوسط تعداد روزهای معاملاتی در سال
    n_active_Years=n_active/n_years;
    n_active_Years2=n_active2/n_years;
    DR = exp(mean(log(1 + array))) - 1; % Average daily return
    DR2 = exp(mean(log(1 + Returns))) - 1; % Average daily return

    AR=((1+DR)^n_active_Years)-1;                % Annualised Rate
    AR2=((1+DR2)^n_active_Years2)-1;                % Annualised Rate

    Annual_Return=AR;
    Annual_Return_Complete=AR2;

    Dstd=std(array);                % standard deviation day
    Dstd2=std(Returns);

    SharpeRatio=sqrt(n_active_Years)*(DR./Dstd);
    SharpeRatio_Complete=sqrt(n_active_Years2)*(DR2./Dstd2);


    SharpeRatio(isnan(SharpeRatio))=0;
    SharpeRatio(isinf(SharpeRatio))=0;

    SharpeRatio_Complete(isnan(SharpeRatio_Complete))=0;
    SharpeRatio_Complete(isinf(SharpeRatio_Complete))=0;

    SharpeRatio=real(SharpeRatio);
    SharpeRatio_Complete=real(SharpeRatio_Complete);

    Annual_Return=real(Annual_Return);
    Annual_Return_Complete=real(Annual_Return_Complete);

    cum_ret = cumprod(1 + array);
    max_dd = maxdrawdown(cum_ret);
    max_dd=real(max_dd);
else
    SharpeRatio=0;
    max_dd=0;
    Annual_Return=0;
    SharpeRatio_Complete=0;
    Annual_Return_Complete=0;
end



x100=1;





end




function[RiskManag]=RiskManagement(CLOSE_Ytest,DATE_Ytest,YPred_test)

CLOSE=CLOSE_Ytest;
DATE=DATE_Ytest;
YPred_test_double=double(YPred_test);
YPred_test_double(YPred_test_double==2)=0;
% سیگنال خرید
buy_signal=YPred_test_double;
% سیگنال فروش
sell_signal=~YPred_test_double;

P_n_L= Buy_only_trade_execution_algo2(CLOSE,DATE,buy_signal,sell_signal);
enter_long_price=P_n_L.enter_long_price;
exit_long_price=P_n_L.exit_long_price;

ret = ((log(exit_long_price) - log(enter_long_price))-0.03)';

%% محاسبه درصد موفقیت
if isempty(ret)
    SuccessRate=0;
    NumberOfNegativeTrades=0;
    Kelly=0;
    Profit_LossRatio=0;
    NumberTrades=0;
else
    %% تعداد معاملات
    NumberTrades=numel(ret);
    % جدا کردن بازدهی های مثبت
    positive_numbers=ret(ret>0);
    % جدا کردن بازدهی های منفی
    Negative_numbers=ret(ret<=0);
    % میانگین بازدهی های مثبت
    if isempty(positive_numbers)
        W=0;
    else
        W=mean(positive_numbers);
    end
    % میانگین بازدهی های منفی
    if isempty(Negative_numbers)
        L = 1e-6; % مقدار حداقلی مثبت
    else
        L=mean(Negative_numbers);
    end
    % محاسبه درصد موفقیت
    P=mean(ret>0);
    SuccessRate=P;
    SuccessRate(isnan(SuccessRate))=0;
    SuccessRate(isinf(SuccessRate))=0;

    % نسبت سود به ضرر
    Profit_LossRatio=W/abs(L);
    %% مشخص کردن معاملات منفی
    is_negative =1*(ret < 0);
    % حداکثر منفی متوالی
    NumberOfNegativeTrades=0;
    NOF=0;

    for ii=1:size(is_negative,1)
        if ii==1 && is_negative(ii)==1
            NOF=NOF+1;
            if NOF>NumberOfNegativeTrades
                NumberOfNegativeTrades=NOF;
            end
        elseif ii>1
            if is_negative(ii)==1
                NOF=NOF+1;
                if NOF>NumberOfNegativeTrades
                    NumberOfNegativeTrades=NOF;
                end
            elseif is_negative(ii)==0
                NOF=0;
            end
        end
    end

    NumberOfNegativeTrades(isnan(NumberOfNegativeTrades))=0;
    NumberOfNegativeTrades(isinf(NumberOfNegativeTrades))=0;

    %% حداکثر سرمایه گذاری
     % 'محاسبه حداکثر حجم برمبنای رابطه کلی'
    Kelly=((Profit_LossRatio*SuccessRate)-(1-SuccessRate))/Profit_LossRatio;
end


 ret1=ret';
transaction=[P_n_L.enter_long_time;P_n_L.exit_long_time;P_n_L.enter_long_price;P_n_L.exit_long_price;((log(exit_long_price) - log(enter_long_price)))];
x100=1;
transaction_Table=table(transaction,'RowNames',{'تاریخ ورود','تاریخ خروج','قیمت ورود','فیمت خروج','بازدهی'}); 

RiskManag=struct('SuccessRate',SuccessRate,...
    'Profit_LossRatio',Profit_LossRatio,...
    'NumberOfNegativeTrades',NumberOfNegativeTrades,...
    'Kelly',Kelly,...
    'NumberTrades',NumberTrades,...
    'transaction_Table',transaction_Table);

end


function P_n_L= Buy_only_trade_execution_algo2(x,D,Buy_Signal,Sell_Signal)

%          BUY ONLY - this algorithm is used as a comparison method
%          against buy and hold strategies to assess strategy
%          performance.
%
%          P_n_L       = a generated structure that holds all trade data
%                        attributed to each executed trade
%          x           = vector of closing prices of the security being 
%                        traded
%          D           = vector of dates corresponding to the prices of 
%                        vector x
%          Buy_Signal  = Trading signal, set of instructions to enter a
%                        long position in the underlying asset
%          Sell_Signal = Trading signal, set of instructions to exit a long
%                        position in the underlying asset

                %% Create Profit and Loss structure

%          This structure is used to hold all executed trading signals and
%          related trade prices and times

                        P_n_L = struct ('enter_long_price', [],...
                                        'enter_long_time',[],...
                                        'exit_long_price',[],...
                                        'exit_long_time', [],...
                                        'enter_long_i',[],...
                                        'exit_long_i',[]);
                                    
                    %% GENERATING TRADING SIGNALS
                    
P_n_L.long_count = 0;   % Count used to track trades

for i = 1:(length(x)-1)
    
    if (Buy_Signal(i))  % Buy signal generated
           
        if (P_n_L.long_count == 1 )  
            
                        % Long position currently in play 
                        
            P_n_L.long_count = 1;    
            
        elseif (P_n_L.long_count == 0) 
            
                        % No position, long position initiated and logged. 
                        % Long count = 1
                                       
            P_n_L.enter_long_price = [P_n_L.enter_long_price x(i)];                                           
            P_n_L.enter_long_time = [P_n_L.enter_long_time D(i)];
            P_n_L.enter_long_i = [P_n_L.enter_long_i i];
            P_n_L.long_count = 1;    
            
        end       
        
    end
    
    if (Sell_Signal(i)) % Sell signal generated
        
        if (P_n_L.long_count == 1)
            
                        % Long position currently in play, position is sold
                        % and logged
                        % Long count = 0
                        
            P_n_L.exit_long_price = [P_n_L.exit_long_price x(i)];                        
            P_n_L.exit_long_time = [P_n_L.exit_long_time D(i)];
            P_n_L.exit_long_i = [P_n_L.exit_long_i i];
            P_n_L.long_count = 0;

        elseif (P_n_L.long_count == 0) 
            
                        % No position 
                        
            P_n_L.long_count = 0;
            
        end          
    end      
    
end
 

if (length(P_n_L.exit_long_price) ~= length(P_n_L.enter_long_price))    
    d = length(P_n_L.exit_long_price);    
    P_n_L.enter_long_price = P_n_L.enter_long_price(1:d);    
    P_n_L.enter_long_time = P_n_L.enter_long_time(1:d); 
    P_n_L.enter_long_i = P_n_L.enter_long_i(1:d);
    P_n_L.exit_long_i = P_n_L.exit_long_i(1:d);
end    

    d = length(P_n_L.exit_long_price);
    P_n_L.enter_long_i = P_n_L.enter_long_i(1:d);
    P_n_L.exit_long_i = P_n_L.exit_long_i(1:d);


    
end

