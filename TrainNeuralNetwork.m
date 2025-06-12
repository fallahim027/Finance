function [ACC,NetBest]=TrainNeuralNetwork(CrossValidationData,ForTheModel,numLayers,numHiddenUnits1,numHiddenUnits2,numHiddenUnits3,numHiddenUnits4,numHiddenUnits5,probability1,probability2,probability3,probability4,probability5,probability6,numHiddenUnitsFC1,numEpochs,MiniBatchSize,gradientThreshold,L2Value,LearnRate)

% NetBest=cell(4,1);

% NetBest=struct('net',{},...
%     'SharpeRatio_YPredTest',[],...
%     'max_dd_YPredTest',[],...
%     'Annual_Return_YPredTest',[],...
%     'SharpeRatio_Complete_YPredTest',[],...
%     'Annual_Return_Complete_YPredTest',[]);

% switch numWindows
%     case 1
%         Windows=7;
%     case 2
%         Windows=10;
%     case 3
%         Windows=14;
%     case 4
%         Windows=21;
%     case 5
%         Windows=30;
%     case 6
%         Windows=60;
%     case 7
%         Windows=90;
% end

Windows=30;

%% طراحی شبکه عصبی
n=size(CrossValidationData,1);
AccuracyesValidation=zeros(n,1);
SharpeRatioObjective=zeros(n,1);
AccuracyesTest=zeros(n,1);
SelectionCriteria=zeros(n,1);
AccuracyesTrading=zeros(n,1);

% پاک کردن حافظه gpu
gpuDevice().reset;

mo=2;
switch mo
    case 1
        models="Window";
    case 2
        models="NoWindow";
end
x100=1;
for i=1:size(CrossValidationData,1)
    % پاک کردن حافظه gpu
    gpuDevice().reset;

    Xtrain_Matrix_Encoder=CrossValidationData{i}.Xtrain_Matrix_Encoder;
    XValidation_Matrix_Encoder=CrossValidationData{i}.XValidation_Matrix_Encoder;
    Xtest_Matrix_Encoder=CrossValidationData{i}.Xtest_Matrix_Encoder;
    Ytrain_Matrix=CrossValidationData{i}.Ytrain_Matrix;
    YValidation_Matrix=CrossValidationData{i}.YValidation_Matrix;
    Ytest_Matrix=CrossValidationData{i}.Ytest_Matrix;

    XtrainGPU=Xtrain_Matrix_Encoder;
    % تبدیل داده ها بصورت پنجره ایی

    if strcmp(models,"Window")
        [Xtrain_Window,XValidation_Window,Xtest_Window,Ytrain_Window,YValidation_Window,Ytest_Window]=DataWindowing(Windows,Xtrain_Matrix_Encoder,XValidation_Matrix_Encoder,Xtest_Matrix_Encoder,Ytrain_Matrix,YValidation_Matrix,Ytest_Matrix);
        Xtrain_Matrix_Encoder=Xtrain_Window;
        XValidation_Matrix_Encoder=XValidation_Window;
        Xtest_Matrix_Encoder=Xtest_Window;
        Ytrain_Matrix=Ytrain_Window;
        YValidation_Matrix=YValidation_Window;
        Ytest_Matrix=Ytest_Window;
        XtrainGPU=Xtrain_Matrix_Encoder;
    end

    % تعداد ویژگی ها
    numFeatures=size(CrossValidationData{i}.Xtrain_Matrix_Encoder,1);
    % تعداد کلاس ها
    classes = categories(CrossValidationData{i}.Ytrain_Matrix);
    numClasses=height(classes);
    ClassWeights=[9,1];
x100=1;
    switch numLayers
        case 1
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits2,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability2)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                % reluLayer()
                swishLayer()
                dropoutLayer(probability3)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];
        case 2
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits2,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability2)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits3,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability3)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability4)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];
        case 3
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits2,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability2)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits3,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability3)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits4,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability4)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability5)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];





        case 4
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits2,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability2)
                lstmLayer(numHiddenUnits3,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability3)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability4)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];
        case 5
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits2,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability2)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits3,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability3)
                lstmLayer(numHiddenUnits4,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability4)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability5)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];
        case 6
            %% CNN_biLSTM 1
            layers = [ sequenceInputLayer(numFeatures,Normalization="none",Name="input");
                % لایه اول کانولیشن
                convolution1dLayer(3,numHiddenUnits1,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability1)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits2,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability2)
                maxPooling1dLayer(2,"Padding","same")

                convolution1dLayer(3,numHiddenUnits3,'Padding','same','DilationFactor',2)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability3)
                maxPooling1dLayer(2,"Padding","same")

                flattenLayer()
                lstmLayer(numHiddenUnits4,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability4)
                lstmLayer(numHiddenUnits5,'OutputMode','sequence','InputWeightsInitializer','narrow-normal')
                layerNormalizationLayer()
                dropoutLayer(probability5)

                fullyConnectedLayer(numHiddenUnitsFC1)
                batchNormalizationLayer()
                swishLayer()
                dropoutLayer(probability6)

                fullyConnectedLayer(numClasses)
                softmaxLayer()
                classificationLayer('Classes',classes,'ClassWeights',ClassWeights)];

    end
    %
    % XValidation_Window2=XValidation_Window;
    % YValidation_Window2=YValidation_Window;
    % Xtest_Window2=Xtest_Window;
    % Ytest_Window2=Ytest_Window;
    switch ForTheModel
        case "BayesianModel"
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','gpu',...
                'Shuffle','every-epoch',...
                'MaxEpochs',numEpochs,...
                'MiniBatchSize',MiniBatchSize,...
                'ValidationData',{XValidation_Matrix_Encoder,YValidation_Matrix},...
                'ValidationFrequency',15,...
                'ValidationPatience',5, ...
                'GradientThreshold',gradientThreshold, ...
                'L2Regularization',L2Value,...
                'Verbose',0,...
                'Acceleration','auto',...
                'PreprocessingEnvironment','serial',...
                'InitialLearnRate',LearnRate,...
                'LearnRateSchedule','piecewise',...
                'LearnRateDropFactor',0.1,...
                'LearnRateDropPeriod',round(numEpochs/3),...
                'OutputNetwork','last-iteration',...
                'SequenceLength','shortest');

        case "FinalModel"
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','gpu',...
                'Shuffle','every-epoch',...
                'MaxEpochs',numEpochs,...
                'GradientThreshold',gradientThreshold, ...
                'L2Regularization',L2Value,...
                'Verbose',0,...
                'Acceleration','auto',...
                'PreprocessingEnvironment','serial',...
                'InitialLearnRate',LearnRate,...
                'LearnRateSchedule','piecewise',...
                'LearnRateDropFactor',0.1,...
                'LearnRateDropPeriod',round(numEpochs/3),...
                'OutputNetwork','last-iteration',...
                'SequenceLength','shortest');
    end

    x100=1;
    switch ForTheModel
        case "BayesianModel"
            net = trainNetwork(XtrainGPU,Ytrain_Matrix,layers,options);
x100=1;
            YPred_Validation=classify(net,XValidation_Matrix_Encoder);
            if strcmp(models,"Window")
                YPred_Validation=YPred_Validation';
            end



            AccValidation=sum(YPred_Validation==YValidation_Matrix)./numel(YValidation_Matrix);

            AccuracyesValidation(i,1)=AccValidation;

            % نسبت شارپ
            Returns_YValidation=CrossValidationData{i, 1}.Returns_YValidation;  
            [SharpeRatio,max_dd,Annual_Return,SharpeRatio_Complete,Annual_Return_Complete]=SharpeRatioForNeuralNetworkEvaluation(YPred_Validation,Returns_YValidation);
            SharpeRatioObjective(i,1)=SharpeRatio;

            ff=['SharpeRatio_Validation ',num2str(SharpeRatio),'  max_dd_Validation ',num2str(max_dd),...
                '  Annual_Return_Validation ',num2str(Annual_Return),' SharpeRatio_Complete ',num2str(SharpeRatio_Complete),' Annual_Return_Complete ',num2str(Annual_Return_Complete)];
            disp(ff)

            %% درصد موفقیت داده های تست

            YPred_test = classify(net,Xtest_Matrix_Encoder);
            if strcmp(models,"Window")
                YPred_test=YPred_test';
            end
            AccTest = sum(YPred_test == Ytest_Matrix)./numel(Ytest_Matrix);
x100=1;
            AccuracyesTest(i,1)=AccTest;
            % نسبت شارپ
            Returns_Ytest=CrossValidationData{i, 1}.Returns_Ytest; 
            [SharpeRatio,max_dd,Annual_Return,SharpeRatio_Complete,Annual_Return_Complete]=SharpeRatioForNeuralNetworkEvaluation(YPred_test,Returns_Ytest);
            SharpeRatio_YPredTest=SharpeRatio;
            max_dd_YPredTest=max_dd;
            Annual_Return_YPredTest=Annual_Return;

            %% درصد موفقیت داده های آموزش

            YPred_Ytrain = classify(net,XtrainGPU);
            if strcmp(models,"Window")
                YPred_Ytrain=YPred_Ytrain';
            end
            AccTrading = sum(YPred_Ytrain == Ytrain_Matrix)./numel(Ytrain_Matrix);

            AccuracyesTrading(i,1)=AccTrading;

        case "FinalModel"
            net = trainNetwork(XtrainGPU,Ytrain_Matrix,layers,options);

            YPred_test = classify(net,Xtest_Matrix_Encoder);
            AccTest = sum(YPred_test == Ytest_Matrix)./numel(Ytest_Matrix);

            Returns_Ytest=CrossValidationData{i, 1}.Returns_Ytest;
            [SharpeRatio,max_dd,Annual_Return,SharpeRatio_Complete,Annual_Return_Complete]=SharpeRatioForNeuralNetworkEvaluation(YPred_test,Returns_Ytest);
            SharpeRatio_YPredTest=SharpeRatio;
            max_dd_YPredTest=max_dd;
            Annual_Return_YPredTest=Annual_Return;
            SharpeRatio_Complete_YPredTest=SharpeRatio_Complete;
            Annual_Return_Complete_YPredTest=Annual_Return_Complete;


            % سایر ابزارهای مدیریت ریسک
            CLOSE_Ytest=CrossValidationData{i, 1}.CLOSE_Ytest;
            DATE_Ytest=CrossValidationData{i, 1}.DATE_Ytest;

            [RiskManag]=RiskManagement(CLOSE_Ytest,DATE_Ytest,YPred_test);
            SuccessRate=RiskManag.SuccessRate;
            Profit_LossRatio=RiskManag.Profit_LossRatio;
            NumberOfNegativeTrades=RiskManag.NumberOfNegativeTrades;
            Kelly=RiskManag.Kelly;
            transaction_Table=RiskManag.transaction_Table;
            NumberTrades=RiskManag.NumberTrades;

            x100=1;
    end

x100=1;
end
switch ForTheModel
    case "BayesianModel"
        % ACC=-mean(AccuracyesValidation);
        ACC=-mean(SharpeRatioObjective);
        a1=AccuracyesValidation'
        a2=AccuracyesTest'
    case "FinalModel"
        ACC=AccTest
        NetBest=struct('net',net,...
            'SharpeRatio_YPredTest',SharpeRatio_YPredTest,...
            'max_dd_YPredTest',max_dd_YPredTest,...
            'Annual_Return_YPredTest',Annual_Return_YPredTest,...
            'SharpeRatio_Complete_YPredTest',SharpeRatio_Complete_YPredTest,...
            'Annual_Return_Complete_YPredTest',Annual_Return_Complete_YPredTest,...
            'YPred_test',YPred_test,...
            'SuccessRate',SuccessRate,...
            'Profit_LossRatio',Profit_LossRatio,...
            'NumberOfNegativeTrades',NumberOfNegativeTrades,...
            'Kelly',Kelly,...
            'NumberTrades',NumberTrades,...
            'transaction_Table',transaction_Table,...
            'Date_Of_Start_Of_Test_Data',DATE_Ytest(1),...
            'End_Date_Of_Test_Data',DATE_Ytest(end));
end


x100=1;
end

%%
function [Xtrain_Window,XValidation_Window,Xtest_Window,Ytrain_Window,YValidation_Window,Ytest_Window]=DataWindowing(Windows,Xtrain_Matrix_Encoder,XValidation_Matrix_Encoder,Xtest_Matrix_Encoder,Ytrain_Matrix,YValidation_Matrix,Ytest_Matrix)

% Windows=30;
% داده آموزش پنجره ایی
Xtrain_Window=cell(size(Xtrain_Matrix_Encoder,2)-Windows+1,1);
for i=Windows:size(Xtrain_Matrix_Encoder,2)
    Xtrain_Window{i-Windows+1,1}=Xtrain_Matrix_Encoder(:,i-Windows+1:i);
end
Ytrain_Window=Ytrain_Matrix(1,Windows:end);

% داده اعتبارسنجی پنجره هایی
XValidation_Window=cell(size(XValidation_Matrix_Encoder,2)-Windows+1,1);
for i=Windows:size(XValidation_Matrix_Encoder,2)
    XValidation_Window{i-Windows+1,1}=XValidation_Matrix_Encoder(:,i-Windows+1:i);
end
YValidation_Window=YValidation_Matrix(1,Windows:end);

% داده تست پنجره ایی
Xtest_Window=cell(size(Xtest_Matrix_Encoder,2)-Windows+1,1);
for i=Windows:size(Xtest_Matrix_Encoder,2)
    Xtest_Window{i-Windows+1,1}=Xtest_Matrix_Encoder(:,i-Windows+1:i);
end
Ytest_Window=Ytest_Matrix(1,Windows:end);

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

