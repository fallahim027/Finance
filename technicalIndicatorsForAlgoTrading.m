function[IndicatorsData,TraditionalStrategy,GenFeature]=technicalIndicatorsForAlgoTrading(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

% Movavg Exponential indicator
[MovExponential,signalMovExponential]=MovavgExponentialindicator(CLOSE);

[MovSimple,SignalMovSimple]=MovavgSimpleindicator(CLOSE);

[SignalTwoMovSimple]=TwoMovavgSimpleindicator(CLOSE);

[SignalAO,AO]=AwesomeOscillatorIndicator(HIGH,LOW,CLOSE);

[SignaltwoMovExponential]=TwoMovavgExponentialindicator(CLOSE);

[SignalAC, AC] = ACIndicator(HIGH, LOW, CLOSE);

[SignalThreeMovSimple]=ThreeMovavgSimpleindicator(CLOSE);

[SignalSAR,SAR]=SarIndicator(HIGH,LOW,CLOSE);

[SignalThreeMovExponential]=ThreeMovavgExponentialindicator(CLOSE);

[SignalCci,cci]=CciIndicator(HIGH,LOW,CLOSE);

[SignalCci2]=Cci2Indicator(HIGH,LOW,CLOSE);

[SignalIchimoku, Tenkan, Kijun, SenkouA, SenkouB] = EnhancedIchimokuSimple(HIGH, LOW, CLOSE);

[SignalMACD,MACD,histogram,LineSignal]=Macdindicator(CLOSE);

[SignalADX,SignalADXTrend,MDI_ADX,PDI_ADX,ADX]=AdxIndicator(HIGH,LOW,CLOSE);

[SignalUltimateOscillator,UltimateOscillator]=UltimateOscillatorIndicator(HIGH,LOW,CLOSE);

[SignalUltimateOscillatorTrend]=UltimateOscillatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalUltimateOscillator_AC]=UltimateOscillatorWithACIndicator(HIGH,LOW,CLOSE);

[SignalACTrend]=AcWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalAlligator,Jaw_blue,Teeth_Red,Lips_Green]=AlligatorIndicator(HIGH,LOW);

[SignalAlligatorTrend]=AlligatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalRsi,Rsi]=RsiIndicator(CLOSE);

[SignalRsi2]=Rsi2Indicator(CLOSE);

[SignalMFI,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL);

[SignalMfiMov]=MfiWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalAOTrend]=AwesomeOscillatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalADXTrend2]=AdxWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalWilliams,williams]=RWilliamsIndicator(HIGH,LOW,CLOSE);

[SignalWilliamsWithTrend]=RWilliamsWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalSarTrend]=SarWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalBollinger,Bollinger]=BollingerBondIndicator(CLOSE);

[SignalBollingerlower]=BollingerBandlowerIndicator(CLOSE);

[SignalBollingerMacd]=BollingerBandWithMacdIndicator(CLOSE);

[SignalBollingerMiddle]=BollingerBondMiddleIndicator(CLOSE);

[SignalCciWithTernd]=CciWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalSlowStochosc,FastPercentK,FastPercentD,SlowPercentK,SlowPercentD]=SlowStochasticIndicator(HIGH,LOW,CLOSE);

[SignalFastStochosc,~,~,~,~]=FastStochasticIndicator(HIGH,LOW,CLOSE);

[SignalSlow_StoTrend]=SlowStochasticWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalFastStoTrend]=FastStochasticWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[SignalSlowStoMFI]=SlowStochasticWithMfiIndicator(HIGH,LOW,CLOSE,VOL);

[SignalFastStoMFI]=FastStochasticWithMfiIndicator(HIGH,LOW,CLOSE,VOL);

[SignalTwoMovavgExponentialRsi]=TwoMovavgExponentialWithRsiIndicator(CLOSE);

[Signal2Rsi]=TwoRsiIndicator(CLOSE);

[SignalMacdRsiEma]=RsiMacdEmaIndicator(CLOSE);

[SignalMfiWiliams]=MfiWithRWilliamsIndicator(HIGH,LOW,CLOSE,VOL);

[SignalMfiRsi]=MfiWithRsiIndicator(HIGH,LOW,CLOSE,VOL);

[SignalMacdRsi]=RsiMacdIndicator(CLOSE);

[SignalRsiTrend]=RsiWithTerndIndIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

[OtherTechnicalIndicators]=OtherIndicators(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
%% مهندسی ویژگی
returns=[0;price2ret(CLOSE)];
GenFeature=table(MovExponential.MovExponential1,MovSimple.MovSimple1,AO.AO1,AC.AC1,SAR.SAR1,cci.cci1,...
    Tenkan.Tenkan1,Kijun.Kijun1,SenkouA.SenkouA1,SenkouB.SenkouB1,histogram.histogram1,...
    MACD.MACD1,LineSignal.LineSignal1,MDI_ADX.MDI_ADX1,PDI_ADX.PDI_ADX1,ADX.ADX1,UltimateOscillator.UltimateOscillator1,Jaw_blue.Jaw_blue_Alligator1,...
    Teeth_Red.Teeth_Red_Alligator1,Lips_Green.Lips_Green_Alligator1,Rsi.Rsi1,MFI.MFI1,williams.williams1,Bollinger.Bollinger1,...
    FastPercentK.FastPercentK1,FastPercentD.FastPercentD1,SlowPercentK.SlowPercentK1,SlowPercentD.SlowPercentD1,...
    OtherTechnicalIndicators.ROC1,OtherTechnicalIndicators.dn_aroon1,OtherTechnicalIndicators.up_aroon1,OtherTechnicalIndicators.os_aroon1,...
    OtherTechnicalIndicators.tsi1,OtherTechnicalIndicators.pdi_adx1,OtherTechnicalIndicators.mdi_adx1,OtherTechnicalIndicators.adx1,...
    OtherTechnicalIndicators.t31,OtherTechnicalIndicators.obv1,OtherTechnicalIndicators.cmf1,OtherTechnicalIndicators.force1,...
    OtherTechnicalIndicators.middle_keltner1,OtherTechnicalIndicators.upper_keltner1,OtherTechnicalIndicators.lower_keltner1,...
    OtherTechnicalIndicators.atr1,OtherTechnicalIndicators.VolatilityRatio1,OtherTechnicalIndicators.PP_pivot1,OtherTechnicalIndicators.R1_pivot1,...
    OtherTechnicalIndicators.R2_pivot1,OtherTechnicalIndicators.R3_pivot1,OtherTechnicalIndicators.S1_pivot1,...
    OtherTechnicalIndicators.S2_pivot1,OtherTechnicalIndicators.S3_pivot1,returns,...
    'VariableNames',["MovExponential","MovSimple","AO","AC","SAR","cci","Tenkan","Kijun","SenkouA","SenkouB",...
    "histogramMACDI","MACD","LineSignal","MDI_ADX","PDI_ADX","ADX","UltimateOscillator","Jaw_blue_Alligator",...
    "Teeth_Red_Alligator","Lips_Green_Alligator","Rsi","MFI","williams","Bollinger","FastPercentK","FastPercentD",...
    "SlowPercentK","SlowPercentD","ROC","dn_aroon","up_aroon","os_aroon","tsi","pdi_adx","mdi_adx","adx",...
    "t3","obv","cmf","force","middle_keltner","upper_keltner","lower_keltner","atr","VolatilityRatio1","PP_pivot",...
    "R1_pivot","R2_pivot","R3_pivot","S1_pivot","S2_pivot","S3_pivot","returns"]);

% [T,NewTbl]=genrfeatures(GenFeature,"returns",500);
x100=1;
%% Model output
% returns1=array2table(returns);
% داده های اندیکاتورها
IndicatorsData=[MovExponential,MovSimple,AO,AC,SAR,cci,Tenkan,Kijun,SenkouA,SenkouB,MACD,histogram,LineSignal,MDI_ADX,PDI_ADX,ADX,UltimateOscillator,...
    Jaw_blue,Teeth_Red,Lips_Green,Rsi,MFI,williams,Bollinger,FastPercentK,FastPercentD,SlowPercentK,SlowPercentD,OtherTechnicalIndicators];

% استراتژی های سنتی
TraditionalStrategy=[signalMovExponential,SignalMovSimple,SignalTwoMovSimple,SignalAO,SignaltwoMovExponential,SignalAC,SignalThreeMovSimple,...
    SignalSAR,SignalThreeMovExponential,SignalCci,SignalCci2,SignalIchimoku,SignalMACD,SignalADX,SignalADXTrend,SignalUltimateOscillator,...
    SignalUltimateOscillatorTrend,SignalUltimateOscillator_AC,SignalACTrend,SignalAlligator,SignalAlligatorTrend,SignalRsi,SignalRsi2,...
    SignalMFI,SignalMfiMov,SignalAOTrend,SignalADXTrend2,SignalWilliams,SignalWilliamsWithTrend,SignalSarTrend,SignalBollinger,...
    SignalBollingerlower,SignalBollingerMacd,SignalBollingerMiddle,SignalCciWithTernd,SignalSlowStochosc,SignalFastStochosc,SignalSlow_StoTrend,...
    SignalFastStoTrend,SignalSlowStoMFI,SignalFastStoMFI,SignalTwoMovavgExponentialRsi,Signal2Rsi,SignalMacdRsiEma,SignalMfiWiliams,...
    SignalMfiRsi,SignalMacdRsi,SignalRsiTrend];
x100=1;
end
%% توابع اندیکاتورهای تکنیکال و استراتژی های سنتی معاملاتی

%% Movavg Exponential indicator
function[MovExponential,signalMovExponential]=MovavgExponentialindicator(CLOSE)
[row,~]=size(CLOSE);
period=[7 14 21 34 55 89];
p=size(period,2);
MovExponential=zeros(row,p);
idx=1;
for i=period
    MovExponential(:,idx)=movavg(CLOSE,'exponential',i);
    idx=idx+1;
end
signalMovExponential=1*(CLOSE>MovExponential);

% table
[~, num_cols] = size(MovExponential);
NameIndicator=arrayfun(@(x) sprintf('MovExponential%d',x),1:num_cols,'UniformOutput',false);
MovExponential=array2table(MovExponential,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('signalMovExponential%d',x),1:num_cols,'UniformOutput',false);
signalMovExponential=array2table(signalMovExponential,"VariableNames",NameIndicator);

end

%%
function[MovSimple,SignalMovSimple]=MovavgSimpleindicator(CLOSE)
[row,~]=size(CLOSE);
period=[8 13 21 34 55 60];
p=size(period,2);
MovSimple=zeros(row,p);
idx=1;
for i=period
    MovSimple(:,idx)=movavg(CLOSE,'simple',i);
    idx=idx+1;
end
SignalMovSimple=1*(CLOSE>MovSimple);

[~, num_cols] = size(MovSimple);
NameIndicator=arrayfun(@(x) sprintf('MovSimple%d',x),1:num_cols,'UniformOutput',false);
MovSimple=array2table(MovSimple,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SignalMovSimple%d',x),1:num_cols,'UniformOutput',false);
SignalMovSimple=array2table(SignalMovSimple,"VariableNames",NameIndicator);
end

%% two simple simple
function[SignalTwoMovSimple]=TwoMovavgSimpleindicator(CLOSE)
[row,~]=size(CLOSE);
period1=[7 14 21 34];
period2=[63 89];
SignalTwoMovSimple=zeros(row,size(period1,2)*size(period2,2));
idx=1;
for i=period1
    Mov1=movavg(CLOSE,'simple',i);
    for j=period2
        Mov2=movavg(CLOSE,'simple',j);
        SignalTwoMovSimple(:,idx)=1*(Mov1>Mov2);
        idx=idx+1;
    end
end
[~, num_cols] = size(SignalTwoMovSimple);
NameIndicator=arrayfun(@(x) sprintf('SignalTwoMovSimple%d',x),1:num_cols,'UniformOutput',false);
SignalTwoMovSimple=array2table(SignalTwoMovSimple,"VariableNames",NameIndicator);

end

%%

function [SignalAO, AO] = AwesomeOscillatorIndicator(HIGH, LOW, CLOSE)
% ورودیها:
%   short_periods: برداری از دورههای کوتاه (مثال: [3,5,7])
%   long_periods:  برداری از دورههای بلند (مثال: [30,34,40])
% خروجیها:
%   AO: ماتریس با ابعاد [تعداد دادهها × تعداد استراتژیها]
%   SignalAO: ماتریس سیگنالها با ابعاد [تعداد دادهها × تعداد استراتژیها]
short_periods =[3 7 13];   
long_periods = [34 40];   
[row, ~] = size(CLOSE);
MEDIAN_PRICE = (HIGH + LOW) / 2;

% ایجاد تمام ترکیبات ممکن از دورهها
[short_grid, long_grid] = meshgrid(short_periods, long_periods);
periods = [short_grid(:), long_grid(:)];

% حذف ترکیبات نامعتبر (بلندمدت ≤ کوتاهمدت)
valid_indices = periods(:,2) > periods(:,1);
periods = periods(valid_indices, :);

num_strategies = size(periods, 1);
AO = zeros(row, num_strategies);
SignalAO = zeros(row, num_strategies);

% محاسبه AO و سیگنال برای هر استراتژی
for s = 1:num_strategies
    short = periods(s, 1);
    long = periods(s, 2);
    
    % محاسبه میانگینهای متحرک
    SMA_short = movavg(MEDIAN_PRICE, 'simple', short);
    SMA_long = movavg(MEDIAN_PRICE, 'simple', long);
    
    % محاسبه AO
    AO(:, s) = SMA_short - SMA_long;
    
    % تولید سیگنال
    for i = 2:row
        if AO(i, s) < 0 && AO(i-1, s) < AO(i, s)
            SignalAO(i, s) = 1; % سیگنال خرید
        elseif AO(i, s) >= 0 && AO(i-1, s) > AO(i, s)
            SignalAO(i, s) = 0; % سیگنال فروش
        else
            SignalAO(i, s) = SignalAO(i-1, s); % حفظ وضعیت
        end
    end
end

[~, num_cols] = size(AO);
NameIndicator=arrayfun(@(x) sprintf('AO%d',x),1:num_cols,'UniformOutput',false);
AO=array2table(AO,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SignalAO%d',x),1:num_cols,'UniformOutput',false);
SignalAO=array2table(SignalAO,"VariableNames",NameIndicator);
x100=1;
end

%% two simple simple
function[SignaltwoMovExponential]=TwoMovavgExponentialindicator(CLOSE)
[row,~]=size(CLOSE);
period1=[8 13 21 34];
period2=[55 89];
SignaltwoMovExponential=zeros(row,size(period1,2)*size(period2,2));
idx=1;
for i=period1
    Mov1=movavg(CLOSE,'exponential',i);
    for j=period2
        Mov2=movavg(CLOSE,'exponential',j);
        SignaltwoMovExponential(:,idx)=1*(Mov1>Mov2);
        idx=idx+1;
    end
end
[~, num_cols] = size(SignaltwoMovExponential);
NameIndicator=arrayfun(@(x) sprintf('SignaltwoMovExponential%d',x),1:num_cols,'UniformOutput',false);
SignaltwoMovExponential=array2table(SignaltwoMovExponential,"VariableNames",NameIndicator);
end

%%
function [SignalAC, AC] = ACIndicator(HIGH, LOW, CLOSE)
% تنظیمات ثابت برای سادگی
fast_periods = [5 8];    % دورههای سریع
slow_periods = [21 34];  % دورههای کند
smooth_periods = [3 5];  % دورههای هموارسازی

MEDIAN_PRICE = (HIGH + LOW)/2;
num_strategies = length(fast_periods) * length(slow_periods) * length(smooth_periods);
[row, ~] = size(CLOSE);

AC = zeros(row, num_strategies);
SignalAC = zeros(row, num_strategies);
idx = 1;

for fast = fast_periods
    for slow = slow_periods
        % رد ترکیبات نامعتبر
        if slow <= fast 
            continue; 
        end 
        for smooth = smooth_periods
            % محاسبه AC
            AO = movavg(MEDIAN_PRICE, 'simple', fast) - movavg(MEDIAN_PRICE, 'simple', slow);
            AC(:, idx) = AO - movavg(AO, 'simple', smooth);
            
            % تولید سیگنال ساده
            SignalAC(:, idx) = generateSimpleSignal(AC(:, idx), CLOSE);
            idx = idx + 1;
        end
    end
end

AC = AC(:, 1:idx-1);
SignalAC = SignalAC(:, 1:idx-1);

[~, num_cols] = size(AC);
NameIndicator=arrayfun(@(x) sprintf('AC%d',x),1:num_cols,'UniformOutput',false);
AC=array2table(AC,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SignalAC%d',x),1:num_cols,'UniformOutput',false);
SignalAC=array2table(SignalAC,"VariableNames",NameIndicator);

end

function signals = generateSimpleSignal(AC, CLOSE)
signals = zeros(size(AC));
for i = 2:length(AC)
    if AC(i) > 0 && AC(i-1) <= 0 && CLOSE(i) > CLOSE(i-1)
        signals(i) = 1; % خرید
    elseif AC(i) < 0 && AC(i-1) >= 0 && CLOSE(i) < CLOSE(i-1)
        signals(i) =0; % فروش
    else
        signals(i) = signals(i-1); % بدون تغییر
    end
end
end

%%
function[SignalThreeMovSimple]=ThreeMovavgSimpleindicator(CLOSE)

%% three exponential simple
[row,~]=size(CLOSE);
short_periods = [8 13];    % اعداد فیبوناچی کوتاه‌مدت
medium_periods = [21 34]; % اعداد فیبوناچی میان‌مدت
long_periods =[55 89]; % اعداد فیبوناچی بلندمدت
SignalThreeMovSimple=zeros(row,size(short_periods,2)*size(medium_periods,2)*size(long_periods,2));
idx=1;
for i=short_periods
    Mov1=movavg(CLOSE,'simple',i);
    for j=medium_periods
        Mov2=movavg(CLOSE,'simple',j);
        for h=long_periods
            Mov3=movavg(CLOSE,'simple',h);
            SignalThreeMovSimple(:,idx)=1*(Mov1>Mov2 & Mov2>Mov3);
        end
        idx=idx+1;
    end
end
[~, num_cols] = size(SignalThreeMovSimple);
NameIndicator=arrayfun(@(x) sprintf('SignalThreeMovSimple%d',x),1:num_cols,'UniformOutput',false);
SignalThreeMovSimple=array2table(SignalThreeMovSimple,"VariableNames",NameIndicator);

end

%%
function [SignalSAR,SAR]=SarIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

 % Model

 AFstart=[0.01 0.02 0.05];
 AFmax=[0.15  0.3 0.5];
 AFdelta=[0.01 0.02 0.05];
 idx=1;
 SignalSAR=zeros(row,size(AFstart,2)*size(AFmax,2)*size(AFdelta,2));
 SAR=zeros(row,size(AFstart,2)*size(AFmax,2)*size(AFdelta,2));
 for j=AFstart
     for jj=AFmax
         for jjj=AFdelta
             [Sar,turnPoints,longSort]= parabolicSAR(HIGH,LOW,j,jj,jjj);
             Sar=Sar';
             turnPoints=turnPoints';
             longSort=longSort';
             for i=2:row
                 if CLOSE(i-1,1) < Sar(i-1,1) && CLOSE(i,1) > Sar(i,1)
                     SignalSAR(i,idx)=1;
                 elseif CLOSE(i-1,1) > Sar(i-1,1) && CLOSE(i,1) < Sar(i,1)
                     SignalSAR(i,idx)=0;
                 else
                     SignalSAR(i,idx)=SignalSAR(i-1,idx);
                 end
             end
             SAR(:,idx)=Sar;
             idx=idx+1;
         end
     end
 end
 
 [~, num_cols] = size(SAR);
 NameIndicator=arrayfun(@(x) sprintf('SAR%d',x),1:num_cols,'UniformOutput',false);
 SAR=array2table(SAR,"VariableNames",NameIndicator);

 NameIndicator=arrayfun(@(x) sprintf('SignalSAR%d',x),1:num_cols,'UniformOutput',false);
 SignalSAR=array2table(SignalSAR,"VariableNames",NameIndicator);

end

%%
function[SignalThreeMovExponential]=ThreeMovavgExponentialindicator(CLOSE)

%% three exponential simple
[row,~]=size(CLOSE);
short_periods = [7, 14];    % اعداد فیبوناچی کوتاه‌مدت
medium_periods = [21, 34]; % اعداد فیبوناچی میان‌مدت
long_periods =[55 89]; % اعداد فیبوناچی بلندمدت
SignalThreeMovExponential=zeros(row,size(short_periods,2)*size(medium_periods,2)*size(long_periods,2));
idx=1;
for i=short_periods
    Mov1=movavg(CLOSE,'exponential',i);
    for j=medium_periods
        Mov2=movavg(CLOSE,'exponential',j);
        for h=long_periods
            Mov3=movavg(CLOSE,'exponential',h);
            SignalThreeMovExponential(:,idx)=1*(Mov1>Mov2 & Mov2>Mov3);
        end
        idx=idx+1;
    end
end
[~, num_cols] = size(SignalThreeMovExponential);
NameIndicator=arrayfun(@(x) sprintf('SignalThreeMovExponential%d',x),1:num_cols,'UniformOutput',false);
SignalThreeMovExponential=array2table(SignalThreeMovExponential,"VariableNames",NameIndicator);
end

%%
function[SignalCci,cci]=CciIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

% Model
Period=[7 13 21 34 63];
cci=zeros(row,size(Period,2));
idx=1;
for i=Period
    cci(:,idx)= indicators([HIGH,LOW,CLOSE],'cci',i,i,0.015);
    idx=idx+1;
end
% cci=fillmissing(cci,"linear","EndValues","nearest");
col=size(cci,2);
SignalCci=zeros(row,size(Period,2));
for j=1:col
    for i=2:row
        if cci(i-1,j)<-100 && cci(i,j)>-100
            SignalCci(i,j)=1;
        elseif cci(i-1,j)>100 && cci(i,j)<100
            SignalCci(i,j)=0;
        else
            SignalCci(i,j)=SignalCci(i-1,j);
        end
    end
end
[~, num_cols] = size(cci);
NameIndicator=arrayfun(@(x) sprintf('cci%d',x),1:num_cols,'UniformOutput',false);
cci=array2table(cci,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SignalCci%d',x),1:num_cols,'UniformOutput',false);
SignalCci=array2table(SignalCci,"VariableNames",NameIndicator);
end
%%
function[SignalCci2]=Cci2Indicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);
[~,cci]=CciIndicator(HIGH,LOW,CLOSE);
cci=table2array(cci);
[~,col_cci]=size(cci);

% Model
SignalCci2=zeros(row,col_cci);
for j=1:col_cci
    SignalCci2(1,j)=0;
    for i=2:row
        if cci(i-1,j)<0 && cci(i,j)>0
            SignalCci2(i,j)=1;
        elseif cci(i-1,j)>0 && cci(i,j)<0
            SignalCci2(i,j)=0;
        else
            SignalCci2(i,j)=SignalCci2(i-1,j);
        end
    end
end
[~, num_cols] = size(SignalCci2);
NameIndicator=arrayfun(@(x) sprintf('SignalCci2%d',x),1:num_cols,'UniformOutput',false);
SignalCci2=array2table(SignalCci2,"VariableNames",NameIndicator);
end

%%
function [SignalIchimoku, Tenkan, Kijun, SenkouA, SenkouB] = EnhancedIchimokuSimple(HIGH, LOW, CLOSE)
    % تنظیمات پیشفرض استاندارد ایچیموکو
    tenkan_period = 9;
    kijun_period = 26;
    senkouB_period = 52;
    senkou_shift = 26;
    
    [n, ~] = size(CLOSE);
    
    % پیش‌تخصیص حافظه برای بهبود سرعت
    Tenkan = zeros(n,1);
    Kijun = zeros(n,1);
    SenkouA = nan(n,1);
    SenkouB = nan(n,1);
    SignalIchimoku = zeros(n,1);
    
    % محاسبات برداری برای بهبود کارایی
    % تبدیل خط تنکان
    for i = tenkan_period:n
        Tenkan(i) = (max(HIGH(i-tenkan_period+1:i)) + min(LOW(i-tenkan_period+1:i))) / 2;
    end
    
    % تبدیل خط کیجون
    for i = kijun_period:n
        Kijun(i) = (max(HIGH(i-kijun_period+1:i)) + min(LOW(i-kijun_period+1:i))) / 2;
    end
    
    % ابر سنکو (A و B)
    SenkouA(senkou_shift+1:end) = (Tenkan(1:end-senkou_shift) + Kijun(1:end-senkou_shift)) / 2;
    
    for i = senkouB_period+senkou_shift:n
        SenkouB(i) = (max(HIGH(i-senkouB_period-senkou_shift+1:i-senkou_shift)) + ...
                     min(LOW(i-senkouB_period-senkou_shift+1:i-senkou_shift))) / 2;
    end
    
    % منطق سیگنالدهی بهبود یافته
    for i = senkouB_period+senkou_shift+1:n
        % شرایط اصلی خرید
        cond1 = CLOSE(i) > max(SenkouA(i), SenkouB(i)); % قیمت بالای ابر
        cond2 = Tenkan(i) > Kijun(i); % خط تنکان بالای خط کیجون
        cond3 = SenkouA(i) > SenkouB(i); % ابر سبز رنگ
        
        % شرایط اصلی فروش
        cond4 = CLOSE(i) < min(SenkouA(i), SenkouB(i)); % قیمت زیر ابر
        cond5 = Tenkan(i) < Kijun(i); % خط تنکان زیر خط کیجون
        cond6 = SenkouA(i) < SenkouB(i); % ابر قرمز رنگ
        
        % تولید سیگنال با فیلتر ساده
        if cond1 && cond2 && cond3
            SignalIchimoku(i) = 1; % سیگنال خرید قوی
        elseif cond4 && cond5 && cond6
            SignalIchimoku(i) =0; % سیگنال فروش قوی
        else
            SignalIchimoku(i) = SignalIchimoku(i-1); % حفظ وضعیت قبلی
        end
    end
    
    % حذف سیگنال‌های اولیه نامعتبر
    SignalIchimoku(1:senkouB_period+senkou_shift) = 0;

    [~, num_cols] = size(SignalIchimoku);
    NameIndicator=arrayfun(@(x) sprintf('SignalIchimoku%d',x),1:num_cols,'UniformOutput',false);
    SignalIchimoku=array2table(SignalIchimoku,"VariableNames",NameIndicator);

    NameIndicator=arrayfun(@(x) sprintf('Tenkan%d',x),1:num_cols,'UniformOutput',false);
    Tenkan=array2table(Tenkan,"VariableNames",NameIndicator);

    NameIndicator=arrayfun(@(x) sprintf('Kijun%d',x),1:num_cols,'UniformOutput',false);
    Kijun=array2table(Kijun,"VariableNames",NameIndicator);

    NameIndicator=arrayfun(@(x) sprintf('SenkouA%d',x),1:num_cols,'UniformOutput',false);
    SenkouA=array2table(SenkouA,"VariableNames",NameIndicator);

    NameIndicator=arrayfun(@(x) sprintf('SenkouB%d',x),1:num_cols,'UniformOutput',false);
    SenkouB=array2table(SenkouB,"VariableNames",NameIndicator);

end


%%
function[SignalMACD,MACD,histogram,LineSignal]=Macdindicator(CLOSE)

[row,~]=size(CLOSE);

SignalMACD=[];
WindowSizeSignal=9;
WindowSizeFast=12;
WindowSizeSlow=26;
movfast=movavg(CLOSE,'exponential',WindowSizeFast);
movslow=movavg(CLOSE,'exponential',WindowSizeSlow);
MACD=movfast-movslow;
LineSignal=movavg(MACD,'exponential',WindowSizeSignal);
histogram=MACD-LineSignal;
for q=2:row
    SignalMACD(1,1)=0;
    if MACD(q,1)<=0 && MACD(q,1)>LineSignal(q,1)
        SignalMACD(q,1)=1;
    elseif MACD(q,1)>0 && MACD(q,1)<LineSignal(q,1)
        SignalMACD(q,1)=0;
    else
        SignalMACD(q,1)=SignalMACD(q-1,1);
    end
end
[~, num_cols] = size(SignalMACD);
NameIndicator=arrayfun(@(x) sprintf('SignalMACD%d',x),1:num_cols,'UniformOutput',false);
SignalMACD=array2table(SignalMACD,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('MACD%d',x),1:num_cols,'UniformOutput',false);
MACD=array2table(MACD,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('histogram%d',x),1:num_cols,'UniformOutput',false);
histogram=array2table(histogram,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('LineSignal%d',x),1:num_cols,'UniformOutput',false);
LineSignal=array2table(LineSignal,"VariableNames",NameIndicator);

end

%%

function[SignalADX,SignalADXTrend,MDI_ADX,PDI_ADX,ADX]=AdxIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

% Model

[ ADX, ADXR, PDI_ADX, MDI_ADX] = calcDMI(HIGH,LOW,CLOSE);
ADX=ADX';
ADXR=ADXR';
PDI_ADX=PDI_ADX';   %+DM
MDI_ADX=MDI_ADX';   %-DM 
SignalADX=[];
 for i=2:row
     SignalADX(1,1)=0;
     if MDI_ADX(i-1,1) > PDI_ADX(i-1,1) && MDI_ADX(i,1) < PDI_ADX(i,1) && ADX(i,1)>20 
         SignalADX(i,1)=1;
     elseif MDI_ADX(i-1,1) < PDI_ADX(i-1,1) && MDI_ADX(i,1) > PDI_ADX(i,1)
         SignalADX(i,1)=0;
     else
         SignalADX(i,1)=SignalADX(i-1,1);
     end
 end
 SignalADXTrend=1*(MDI_ADX < PDI_ADX & ADX>20 );

 [~, num_cols] = size(SignalADX);
 NameIndicator=arrayfun(@(x) sprintf('SignalADX%d',x),1:num_cols,'UniformOutput',false);
 SignalADX=array2table(SignalADX,"VariableNames",NameIndicator);

 NameIndicator=arrayfun(@(x) sprintf('SignalADXTrend%d',x),1:num_cols,'UniformOutput',false);
 SignalADXTrend=array2table(SignalADXTrend,"VariableNames",NameIndicator);

 NameIndicator=arrayfun(@(x) sprintf('MDI_ADX%d',x),1:num_cols,'UniformOutput',false);
 MDI_ADX=array2table(MDI_ADX,"VariableNames",NameIndicator);

 NameIndicator=arrayfun(@(x) sprintf('PDI_ADX%d',x),1:num_cols,'UniformOutput',false);
 PDI_ADX=array2table(PDI_ADX,"VariableNames",NameIndicator);

 NameIndicator=arrayfun(@(x) sprintf('ADX%d',x),1:num_cols,'UniformOutput',false);
 ADX=array2table(ADX,"VariableNames",NameIndicator);

end

%%
function[SignalUltimateOscillator,UltimateOscillator]=UltimateOscillatorIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);
% Model

BP=zeros(row,1);
TR=zeros(row,1);
for i=2:row
    BP(i,1)=CLOSE(i,1)-min(LOW(i,1),CLOSE(i-1,1));
    TR(i,1)=max(HIGH(i,1),CLOSE(i-1,1))-min(LOW(i,1),CLOSE(i-1,1));
end
 
AVG7=zeros(row,1);
AVG14=zeros(row,1);
AVG28=zeros(row,1);
for i=28:row
    AVG7(i,1)=(sum(BP(i-6:i,1))/sum(TR(i-6:i,1)));
    AVG14(i,1)=(sum(BP(i-13:i,1))/sum(TR(i-13:i,1)));
    AVG28(i,1)=(sum(BP(i-27:i,1))/sum(TR(i-27:i,1)));
end

UltimateOscillator=zeros(row,1);
for i=28:row
    UltimateOscillator(i,1)=(((AVG7(i,1)*4)+(AVG14(i,1)*2)+AVG28(i,1))/(4+2+1))*100;
end

SignalUltimateOscillator=zeros(row,1);
for i=29:row
    if UltimateOscillator(i-1,1) <=30 && UltimateOscillator(i,1) >30 
        SignalUltimateOscillator(i,1)=1;
    elseif UltimateOscillator(i-1,1) >=70 && UltimateOscillator(i,1) <70
        SignalUltimateOscillator(i,1)=0;
    else
        if SignalUltimateOscillator(i-1,1)==1
            SignalUltimateOscillator(i,1)=1;
        end
    end
end
[~, num_cols] = size(SignalUltimateOscillator);
NameIndicator=arrayfun(@(x) sprintf('SignalUltimateOscillator%d',x),1:num_cols,'UniformOutput',false);
SignalUltimateOscillator=array2table(SignalUltimateOscillator,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('UltimateOscillator%d',x),1:num_cols,'UniformOutput',false);
UltimateOscillator=array2table(UltimateOscillator,"VariableNames",NameIndicator);

end

%%
function[SignalUltimateOscillatorTrend]=UltimateOscillatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,UltimateOscillator]=UltimateOscillatorIndicator(HIGH,LOW,CLOSE);
UltimateOscillator=table2array(UltimateOscillator);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);


% Model

SignalUltimateOscillatorTrend=zeros(row,col_MovRavand);

for h=1:col_MovRavand
    for i=29:row
        if UltimateOscillator(i-1,1) <=30 && UltimateOscillator(i,1) >30 && (MovTrend(i,h)==1)
            SignalUltimateOscillatorTrend(i,h)=1;
        elseif UltimateOscillator(i-1,1) >=70 && UltimateOscillator(i,1) <70
            SignalUltimateOscillatorTrend(i,h)=0;
        else
            if SignalUltimateOscillatorTrend(i-1,h)==1
                SignalUltimateOscillatorTrend(i,h)=1;
            end
        end
    end
end
[~, num_cols] = size(SignalUltimateOscillatorTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalUltimateOscillatorTrend%d',x),1:num_cols,'UniformOutput',false);
SignalUltimateOscillatorTrend=array2table(SignalUltimateOscillatorTrend,"VariableNames",NameIndicator);

end


%%
function[SignalUltimateOscillator_AC]=UltimateOscillatorWithACIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);
[~,UltimateOscillator]=UltimateOscillatorIndicator(HIGH,LOW,CLOSE);
UltimateOscillator=table2array(UltimateOscillator);

[~, AC] = ACIndicator(HIGH, LOW, CLOSE);
AC=table2array(AC);
col_AC=size(AC,2);

% Model

SignalUltimateOscillator_AC=zeros(row,col_AC);
idx=1;
for j=1:col_AC
    for i=29:row
        if UltimateOscillator(i-1,1) <=30 && UltimateOscillator(i,1) >30
            SignalUltimateOscillator_AC(i,idx)=1;
        elseif AC(i-1,j)>0 && AC(i,j)<0
            SignalUltimateOscillator_AC(i,idx)=0;
        else
            if SignalUltimateOscillator_AC(i-1,idx)==1
                SignalUltimateOscillator_AC(i,idx)=1;
            end
        end
    end
    idx=idx+1;
end
[~, num_cols] = size(SignalUltimateOscillator_AC);
NameIndicator=arrayfun(@(x) sprintf('SignalUltimateOscillator_AC%d',x),1:num_cols,'UniformOutput',false);
SignalUltimateOscillator_AC=array2table(SignalUltimateOscillator_AC,"VariableNames",NameIndicator);

end

%%
function[SignalACTrend]=AcWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~, AC] = ACIndicator(HIGH, LOW, CLOSE);
AC=table2array(AC);
col_AC=size(AC,2);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);

% Model

SignalACTrend=zeros(row,col_AC*col_MovRavand);
idx=1;
for j=1:col_MovRavand
    for jj=1:col_AC
        for i=2:row
            if MovTrend(i,j)==1 && AC(i-1,jj)<0 && AC(i,jj)>0
                SignalACTrend(i,idx)=1;
            elseif AC(i-1,jj)>0 && AC(i,jj)<0
                SignalACTrend(i,idx)=0;
            else
                SignalACTrend(i,idx)=SignalACTrend(i-1,idx);
            end
        end
        idx=idx+1;
    end
end
[~, num_cols] = size(SignalACTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalACTrend%d',x),1:num_cols,'UniformOutput',false);
SignalACTrend=array2table(SignalACTrend,"VariableNames",NameIndicator);

end

%%
function[SignalAlligator,Jaw_blue_Alligator,Teeth_Red_Alligator,Lips_Green_Alligator]=AlligatorIndicator(HIGH,LOW)

% Model

Jaw_blue_Alligator=movavg((HIGH+LOW),'simple',13)/2;
Teeth_Red_Alligator=movavg((HIGH+LOW),'simple',8)/2;
Lips_Green_Alligator=movavg((HIGH+LOW),'simple',5)/2;
[row_Jaw_blue,~]=size(Jaw_blue_Alligator);
SignalAlligator=zeros(row_Jaw_blue,1);
for i=14:row_Jaw_blue
    if Lips_Green_Alligator(i-3,1) > Lips_Green_Alligator(i-5,1) && Lips_Green_Alligator(i-5,1) > Jaw_blue_Alligator(i-13,1)
        SignalAlligator(i,1)=1;
    elseif Lips_Green_Alligator(i-3,1) < Lips_Green_Alligator(i-5,1) && Lips_Green_Alligator(i-5,1) < Jaw_blue_Alligator(i-13,1)
        SignalAlligator(i,1)=0;
    else
        SignalAlligator(i,1)=SignalAlligator(i-1,1);
    end
end
[~, num_cols] = size(SignalAlligator);
NameIndicator=arrayfun(@(x) sprintf('SignalAlligator%d',x),1:num_cols,'UniformOutput',false);
SignalAlligator=array2table(SignalAlligator,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('Jaw_blue_Alligator%d',x),1:num_cols,'UniformOutput',false);
Jaw_blue_Alligator=array2table(Jaw_blue_Alligator,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('Teeth_Red_Alligator%d',x),1:num_cols,'UniformOutput',false);
Teeth_Red_Alligator=array2table(Teeth_Red_Alligator,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('Lips_Green_Alligator%d',x),1:num_cols,'UniformOutput',false);
Lips_Green_Alligator=array2table(Lips_Green_Alligator,"VariableNames",NameIndicator);

end

%%
function[SignalAlligatorTrend]=AlligatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);

[~,Jaw_blue,~,Lips_Green]=AlligatorIndicator(HIGH,LOW);
Jaw_blue=table2array(Jaw_blue);
Lips_Green=table2array(Lips_Green);

[row_Jaw_blue,col_Jaw_blue]=size(Jaw_blue);

% Model

SignalAlligatorTrend=zeros(row,col_MovRavand*col_Jaw_blue);
for h=1:col_MovRavand
    for i=14:row_Jaw_blue
        if Lips_Green(i-3,1) > Lips_Green(i-5,1) && Lips_Green(i-5,1) > Jaw_blue(i-13,1) && MovTrend(i,h)==1
            SignalAlligatorTrend(i,h)=1;
        elseif Lips_Green(i-3,1) < Lips_Green(i-5,1) && Lips_Green(i-5,1) < Jaw_blue(i-13,1)
            SignalAlligatorTrend(i,h)=0;
        else
            SignalAlligatorTrend(i,h)=SignalAlligatorTrend(i-1,h);
        end
    end
end
[~, num_cols] = size(SignalAlligatorTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalAlligatorTrend%d',x),1:num_cols,'UniformOutput',false);
SignalAlligatorTrend=array2table(SignalAlligatorTrend,"VariableNames",NameIndicator);

end

%%
function[SignalRsi,Rsi]=RsiIndicator(CLOSE)

[row,~]=size(CLOSE);
% Model

period=[7 14 21 34 55];
Rsi=zeros(row,size(period,2));
idx=1;
for j=period
    Rsi(:,idx)=rsindex(CLOSE,j);
    idx=idx+1;
end
% Rsi=fillmissing(Rsi,"linear",1,"EndValues","nearest");
col_Rsi=size(Rsi,2);
SignalRsi=zeros(row,col_Rsi);

for j=1:col_Rsi
    for i=2:row
        % SignalRsi(1,j)=0;
        if Rsi(i-1,j)<30 && Rsi(i,j)>30 
            SignalRsi(i,j)=1;
        elseif Rsi(i-1,j)>70 && Rsi(i,j)<70
            SignalRsi(i,j)=0;
        else
            SignalRsi(i,j)=SignalRsi(i-1,j);
        end
    end
end
[~, num_cols] = size(SignalRsi);
NameIndicator=arrayfun(@(x) sprintf('SignalRsi%d',x),1:num_cols,'UniformOutput',false);
SignalRsi=array2table(SignalRsi,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('Rsi%d',x),1:num_cols,'UniformOutput',false);
Rsi=array2table(Rsi,"VariableNames",NameIndicator);

end

%%
function[SignalRsi2]=Rsi2Indicator(CLOSE)

[row,~]=size(CLOSE);

% Model
period=[7 14 21 34 55];
rsi_str2=zeros(row,size(period,2));
idx=1;
for j=period
    rsi_str2(:,idx)=rsindex(CLOSE,j);
    idx=idx+1;
end
[~,Col_Rsi]=size(rsi_str2);
SignalRsi2=zeros(row,Col_Rsi);
for j=1:Col_Rsi
    for i=2:row
        if rsi_str2(i-1,j)<50 && rsi_str2(i,j)>50 
            SignalRsi2(i,j)=1;
        elseif rsi_str2(i-1,j)>50 && rsi_str2(i,j)<50
            SignalRsi2(i,j)=0;
        else
            SignalRsi2(i,j)=SignalRsi2(i-1,j);
        end
    end
end
[~, num_cols] = size(SignalRsi2);
NameIndicator=arrayfun(@(x) sprintf('SignalRsiTwo%d',x),1:num_cols,'UniformOutput',false);
SignalRsi2=array2table(SignalRsi2,"VariableNames",NameIndicator);

end

%%
function[SignalMFI,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL)

row=size(CLOSE,1);
% Model
MFI=[];
for j=[8 14 21 34 55]
    [~,col_MFI]=size(MFI);
    MFI(:,col_MFI+1)=indicators([HIGH,LOW,CLOSE,VOL],'mfi',j);
end
% MFI=fillmissing(MFI,"linear",1,"EndValues","nearest");
[~,col_MFI]=size(MFI);
SignalMFI = zeros(row,col_MFI);

for j=1:col_MFI
    for i=2:row
        if MFI(i,j)>20 && MFI(i-1,j)<20
            SignalMFI(i,j)=1;
        elseif MFI(i,j)<80 && MFI(i-1,j)>80
            SignalMFI(i,j)=0;
        else
            SignalMFI(i,j)=SignalMFI(i-1,j);
        end
    end
end
[~, num_cols] = size(SignalMFI);
NameIndicator=arrayfun(@(x) sprintf('SignalMFI%d',x),1:num_cols,'UniformOutput',false);
SignalMFI=array2table(SignalMFI,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('MFI%d',x),1:num_cols,'UniformOutput',false);
MFI=array2table(MFI,"VariableNames",NameIndicator);

end

%%
function [SignalMfiTrend] = MfiWithTrendIndicator(DATE, OPEN, HIGH, LOW, CLOSE, VOL)
[row, ~] = size(CLOSE);

% محاسبه MFI و روند
[~, MFI] = MfiIndicator(HIGH, LOW, CLOSE, VOL);
MFI=table2array(MFI);
MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

% مدیریت NaN
% MFI = fillmissing(MFI, 'linear', 'EndValues', 'nearest');

% پیش‌تخصیص ماتریس خروجی
[~, col_MFI] = size(MFI);
[~, col_MovTrend] = size(MovTrend);
SignalMfiTrend = zeros(row, col_MFI * col_MovTrend);

% تولید سیگنال
signal_idx = 1;
for h = 1:col_MovTrend
    for j = 1:col_MFI
        for i = 2:row
            % شرایط خرید
            if MovTrend(i, h) == 1 && MFI(i, j) > 20 && MFI(i-1, j) <= 20
                SignalMfiTrend(i, signal_idx) = 1;
                % شرایط فروش
            elseif MFI(i,j)<80 && MFI(i-1,j)>80
                SignalMfiTrend(i, signal_idx) =0;
            else
                SignalMfiTrend(i, signal_idx) = SignalMfiTrend(i-1, signal_idx);
            end
        end
        signal_idx = signal_idx + 1;
    end
end
[~, num_cols] = size(SignalMfiTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalMfiTrend%d',x),1:num_cols,'UniformOutput',false);
SignalMfiTrend=array2table(SignalMfiTrend,"VariableNames",NameIndicator);

end

%%
function[SignalAOTrend]=AwesomeOscillatorWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~, AO] = AwesomeOscillatorIndicator(HIGH, LOW, CLOSE);
AO=table2array(AO);
colAo=size(AO,2);
MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);

% Model
SignalAOTrend=zeros(row,colAo*col_MovRavand);
signal_idx=1;
for h=1:col_MovRavand
    for A=1:colAo
        for i=2:row
            if (MovTrend(i,h)==1) && AO(i,A)<0 && (AO(i-1,A) < AO(i,A))
                SignalAOTrend(i,signal_idx)=1;
            elseif AO(i,A)>0 && (AO(i-1,A) > AO(i,A))
                SignalAOTrend(i,signal_idx)=0;
            else
                SignalAOTrend(i,signal_idx)=SignalAOTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalAOTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalAOTrend%d',x),1:num_cols,'UniformOutput',false);
SignalAOTrend=array2table(SignalAOTrend,"VariableNames",NameIndicator);

end

%%
function[SignalADXTrend2]=AdxWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,~,MDI,PDI,ADX]=AdxIndicator(HIGH,LOW,CLOSE);
MDI=table2array(MDI);
PDI=table2array(PDI);
ADX=table2array(ADX);
Col_ADX=size(ADX,2);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);
 % Model

 SignalADXTrend2=zeros(row,Col_ADX*col_MovRavand); 
 for h=1:col_MovRavand
     for i=2:row
         if MDI(i-1,1) > PDI(i-1,1) && MDI(i,1) < PDI(i,1) && ADX(i,1)>20 && MovTrend(i,h)==1
             SignalADXTrend2(i,h)=1;
         elseif MDI(i-1,1) < PDI(i-1,1) && MDI(i,1) > PDI(i,1)
             SignalADXTrend2(i,h)=0;
         else
             SignalADXTrend2(i,h)=SignalADXTrend2(i-1,h);
         end
     end
 end
 [~, num_cols] = size(SignalADXTrend2);
 NameIndicator=arrayfun(@(x) sprintf('SignalADXTrendTwo%d',x),1:num_cols,'UniformOutput',false);
 SignalADXTrend2=array2table(SignalADXTrend2,"VariableNames",NameIndicator);
 
end

%%
function[SignalWilliams,williams]=RWilliamsIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

% Model

Periods=[8 14 21 34 55];
williams=zeros(row,size(Periods,2));
Cycle=1;
for i=Periods
    williams(:,Cycle)=indicators([HIGH,LOW,CLOSE],'william',i);
    Cycle=Cycle+1;
end

SignalWilliams=zeros(row,size(williams,2));
col_williams=size(williams,2);
for j=1:col_williams
    for i=2:row
        if williams(i-1,j)<-80 && williams(i,j)>-80
            SignalWilliams(i,j)=1;
        elseif williams(i-1,j)>-20 && williams(i,j)<-20
            SignalWilliams(i,j)=0;
        else
            SignalWilliams(i,j)=SignalWilliams(i-1,j);
        end
    end
end
[~, num_cols] = size(SignalWilliams);
NameIndicator=arrayfun(@(x) sprintf('SignalWilliams%d',x),1:num_cols,'UniformOutput',false);
SignalWilliams=array2table(SignalWilliams,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('williams%d',x),1:num_cols,'UniformOutput',false);
williams=array2table(williams,"VariableNames",NameIndicator);

end

%%
function[SignalWilliamsWithTrend]=RWilliamsWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,williams]=RWilliamsIndicator(HIGH,LOW,CLOSE);
williams=table2array(williams);
[~,col_williams]=size(williams);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);

% Model
SignalWilliamsWithTrend=zeros(row,col_MovRavand*col_williams);
signal_idx=1;
for h=1:col_MovRavand
    for j=1:col_williams
        for i=2:row
            if MovTrend(i,h)==1 && williams(i-1,j)<-80 && williams(i,j)>-80
                SignalWilliamsWithTrend(i,signal_idx)=1;
            elseif williams(i-1,j)>-20 && williams(i,j)<-20
                SignalWilliamsWithTrend(i,signal_idx)=0;
            else
                SignalWilliamsWithTrend(i,signal_idx)=SignalWilliamsWithTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalWilliamsWithTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalWilliamsWithTrend%d',x),1:num_cols,'UniformOutput',false);
SignalWilliamsWithTrend=array2table(SignalWilliamsWithTrend,"VariableNames",NameIndicator);

end

%%
function[SignalSarTrend]=SarWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,SAR]=SarIndicator(HIGH,LOW,CLOSE);
SAR=table2array(SAR);
MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);

 % Model
col_SAR=size(SAR,2);
SignalSarTrend=zeros(row,col_MovRavand*col_SAR);
signal_idx=1;
for h=1:col_MovRavand
    for s=1:col_SAR
        for i=2:row
            if MovTrend(i,h)==1 && CLOSE(i-1,1) < SAR(i-1,s) && CLOSE(i,1) > SAR(i,s)
                SignalSarTrend(i,signal_idx)=1;
            elseif CLOSE(i-1,1) > SAR(i-1,s) && CLOSE(i,1) < SAR(i,s)
                SignalSarTrend(i,signal_idx)=0;
            else
                SignalSarTrend(i,signal_idx)=SignalSarTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalSarTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalSarTrend%d',x),1:num_cols,'UniformOutput',false);
SignalSarTrend=array2table(SignalSarTrend,"VariableNames",NameIndicator);

end

%%
function[SignalBollinger,Bollinger]=BollingerBondIndicator(CLOSE)

[row,~]=size(CLOSE);

% Model
%[middle,upper,lower] = indicators(price           ,'boll'   ,period,weight,nstd)

Period=[14 21 34];
SignalBollinger=zeros(row,size(Period,2));
Bollinger=[];
signal_idx=1;
for j=Period    
    Bollinger1=indicators(CLOSE,'boll',j,0,2);
    % Bollinger1=fillmissing(Bollinger1,"linear","EndValues","nearest");
    Bollinger=[Bollinger1 Bollinger];
    for i=2:row
        if CLOSE(i-1,1)<Bollinger1(i-1,2) && CLOSE(i,1)>Bollinger1(i,2)
            SignalBollinger(i,signal_idx)=1;
        elseif CLOSE(i-1,1)>Bollinger1(i-1,3) && CLOSE(i,1)<Bollinger1(i,3)
            SignalBollinger(i,signal_idx)=0;
        else
            SignalBollinger(i,signal_idx)=SignalBollinger(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalBollinger);
NameIndicator=arrayfun(@(x) sprintf('SignalBollinger%d',x),1:num_cols,'UniformOutput',false);
SignalBollinger=array2table(SignalBollinger,"VariableNames",NameIndicator);

[~, num_cols] = size(Bollinger);
NameIndicator=arrayfun(@(x) sprintf('Bollinger%d',x),1:num_cols,'UniformOutput',false);
Bollinger=array2table(Bollinger,"VariableNames",NameIndicator);

end

%%
function[SignalBollingerlower]=BollingerBandlowerIndicator(CLOSE)

[row,~]=size(CLOSE);
% Model

%[middle,upper,lower] = indicators(price           ,'boll'   ,period,weight,nstd)

Period=[8 13 21 34];
SignalBollingerlower=zeros(row,size(Period,2));
signal_idx=1;
for j=Period
    Bollinger=indicators(CLOSE,'boll',j,0,2);
    % Bollinger=fillmissing(Bollinger,"linear","EndValues","nearest");
    for i=2:row
        if CLOSE(i-1,1)<Bollinger(i-1,3) && CLOSE(i,1)>Bollinger(i,3)
            SignalBollingerlower(i,signal_idx)=1;
        elseif CLOSE(i-1,1)>Bollinger(i-1,3) && CLOSE(i,1)<Bollinger(i,3)
            SignalBollingerlower(i,signal_idx)=0;
        else
            SignalBollingerlower(i,signal_idx)=SignalBollingerlower(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalBollingerlower);
NameIndicator=arrayfun(@(x) sprintf('SignalBollingerlower%d',x),1:num_cols,'UniformOutput',false);
SignalBollingerlower=array2table(SignalBollingerlower,"VariableNames",NameIndicator);

end

%%
function[SignalBollingerMacd]=BollingerBandWithMacdIndicator(CLOSE)

[row,~]=size(CLOSE);
[~,MACD,~,LineSignal]=Macdindicator(CLOSE);
MACD=table2array(MACD);
LineSignal=table2array(LineSignal);
% Model
Period=[7 13 21 34];
SignalBollingerMacd=zeros(row,size(Period,2));

signal_idx=1;
for j=Period
    Bollinger=indicators(CLOSE,'boll',j,0,2);
    % Bollinger=fillmissing(Bollinger,"linear","EndValues","nearest");
    for i=2:row
        if (MACD(i,1)> 0 && MACD(i,1)>LineSignal(i,1)) && (CLOSE(i,1)>Bollinger(i,2))
            SignalBollingerMacd(i,signal_idx)=1;
        elseif MACD(i,1) < LineSignal(i,1) && CLOSE(i,1)<Bollinger(i,3)
            SignalBollingerMacd(i,signal_idx)=0;
        else
            SignalBollingerMacd(i,signal_idx)=SignalBollingerMacd(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalBollingerMacd);
NameIndicator=arrayfun(@(x) sprintf('SignalBollingerMacd%d',x),1:num_cols,'UniformOutput',false);
SignalBollingerMacd=array2table(SignalBollingerMacd,"VariableNames",NameIndicator);

end

%%
function[SignalBollingerMiddle]=BollingerBondMiddleIndicator(CLOSE)

[row,~]=size(CLOSE);
% Model

%[middle,upper,lower] = indicators(price           ,'boll'   ,period,weight,nstd)

Period=[7 13 21 34];
SignalBollingerMiddle=zeros(row,size(Period,2));
signal_idx=1;
for j=Period
    Bollinger=indicators(CLOSE,'boll',j,0,2);
    % Bollinger=fillmissing(Bollinger,"linear","EndValues","nearest");
    for i=2:row
        if CLOSE(i-1,1)<Bollinger(i-1,1) && CLOSE(i,1)>Bollinger(i,1)
            SignalBollingerMiddle(i,signal_idx)=1;
        elseif CLOSE(i-1,1)>Bollinger(i-1,1) && CLOSE(i,1)<Bollinger(i,1)
            SignalBollingerMiddle(i,signal_idx)=0;
        else
            SignalBollingerMiddle(i,signal_idx)=SignalBollingerMiddle(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalBollingerMiddle);
NameIndicator=arrayfun(@(x) sprintf('SignalBollingerMiddle%d',x),1:num_cols,'UniformOutput',false);
SignalBollingerMiddle=array2table(SignalBollingerMiddle,"VariableNames",NameIndicator);

end

%%
function[SignalCciWithTernd]=CciWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,cci]=CciIndicator(HIGH,LOW,CLOSE);
cci=table2array(cci);
MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
% Model
[~,col_cci]=size(cci);
[~,col_MovRavand]=size(MovTrend);
SignalCciWithTernd=zeros(row,col_cci*col_MovRavand);
signal_idx=1;
for h=1:col_MovRavand
    for j=1:col_cci
        for i=2:row
            if MovTrend(i,h)==1 && cci(i-1,j)<-100 && cci(i,j)>-100 
                SignalCciWithTernd(i,signal_idx)=1;
            elseif cci(i-1,j)>100 && cci(i,j)<100
                SignalCciWithTernd(i,signal_idx)=0;
            else
                SignalCciWithTernd(i,signal_idx)=SignalCciWithTernd(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalCciWithTernd);
NameIndicator=arrayfun(@(x) sprintf('SignalCciWithTernd%d',x),1:num_cols,'UniformOutput',false);
SignalCciWithTernd=array2table(SignalCciWithTernd,"VariableNames",NameIndicator);

end

%%
function[SignalSlowStochosc,FastPercentK,FastPercentD,SlowPercentK,SlowPercentD]=SlowStochasticIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

% Model
idx=1;
for NumPeriodsD=[3 9 13]
    for NumPeriodsK=[21 34 55]
        sto= stochosc([HIGH,LOW,CLOSE],'NumPeriodsD',NumPeriodsD,'NumPeriodsK',NumPeriodsK,'Type','exponential');
        FastPercentK(:,idx)=sto(:,1);
        FastPercentD(:,idx)=sto(:,2);
        SlowPercentK(:,idx)=sto(:,3);
        SlowPercentD(:,idx)=sto(:,4);
        idx=idx+1;
    end
end

col_FastPercentK=size(SlowPercentK,2);
SignalSlowStochosc=zeros(row,col_FastPercentK);
signal_idx=1;
for j=1:col_FastPercentK
    for i=2:row
        if SlowPercentK(i,j)<20 && SlowPercentD(i,j)<20 && SlowPercentK(i-1,j)<SlowPercentD(i-1,j) && SlowPercentK(i,j)>SlowPercentD(i,j)
            SignalSlowStochosc(i,signal_idx)=1;
        elseif SlowPercentK(i,j)>80 && SlowPercentD(i,j)>80 && SlowPercentK(i-1,j)>SlowPercentD(i-1,j) && SlowPercentK(i,j)< SlowPercentD(i,j)
            SignalSlowStochosc(i,signal_idx)=0;
        else
            SignalSlowStochosc(i,signal_idx)=SignalSlowStochosc(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalSlowStochosc);
NameIndicator=arrayfun(@(x) sprintf('SignalSlowStochosc%d',x),1:num_cols,'UniformOutput',false);
SignalSlowStochosc=array2table(SignalSlowStochosc,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('FastPercentK%d',x),1:num_cols,'UniformOutput',false);
FastPercentK=array2table(FastPercentK,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('FastPercentD%d',x),1:num_cols,'UniformOutput',false);
FastPercentD=array2table(FastPercentD,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SlowPercentK%d',x),1:num_cols,'UniformOutput',false);
SlowPercentK=array2table(SlowPercentK,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SlowPercentD%d',x),1:num_cols,'UniformOutput',false);
SlowPercentD=array2table(SlowPercentD,"VariableNames",NameIndicator);

end

%%
function[SignalFastStochosc,FastPercentK,FastPercentD,SlowPercentK,SlowPercentD]=FastStochasticIndicator(HIGH,LOW,CLOSE)

[row,~]=size(CLOSE);

% Model

FastPercentK=[];
idx=1;
for NumPeriodsD=[6 9 13]
    for NumPeriodsK=[21 34 55]
        sto= stochosc([HIGH,LOW,CLOSE],'NumPeriodsD',NumPeriodsD,'NumPeriodsK',NumPeriodsK,'Type','exponential');
        FastPercentK(:,idx)=sto(:,1);
        FastPercentD(:,idx)=sto(:,2);
        SlowPercentK(:,idx)=sto(:,3);
        SlowPercentD(:,idx)=sto(:,4);
        idx=idx+1;
    end
end

col_FastPercentK=size(FastPercentK,2);
SignalFastStochosc=zeros(row,col_FastPercentK);
signal_idx=1;
for j=1:col_FastPercentK
    for i=2:row
        if FastPercentK(i,j)<20 && FastPercentD(i,j)<20 && FastPercentK(i-1,j)<FastPercentD(i-1,j) && FastPercentK(i,j)>FastPercentD(i,j)
            SignalFastStochosc(i,signal_idx)=1;
        elseif FastPercentK(i,j)>80 && FastPercentD(i,j)>80 && FastPercentK(i-1,j)>FastPercentD(i-1,j) && FastPercentK(i,j)< FastPercentD(i,j)
            SignalFastStochosc(i,signal_idx)=0;
        else
            SignalFastStochosc(i,signal_idx)=SignalFastStochosc(i-1,signal_idx);
        end
    end
    signal_idx=signal_idx+1;
end
[~, num_cols] = size(SignalFastStochosc);
NameIndicator=arrayfun(@(x) sprintf('SignalFastStochosc%d',x),1:num_cols,'UniformOutput',false);
SignalFastStochosc=array2table(SignalFastStochosc,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('FastPercentK%d',x),1:num_cols,'UniformOutput',false);
FastPercentK=array2table(FastPercentK,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('FastPercentD%d',x),1:num_cols,'UniformOutput',false);
FastPercentD=array2table(FastPercentD,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SlowPercentK%d',x),1:num_cols,'UniformOutput',false);
SlowPercentK=array2table(SlowPercentK,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('SlowPercentD%d',x),1:num_cols,'UniformOutput',false);
SlowPercentD=array2table(SlowPercentD,"VariableNames",NameIndicator);

end

%%
function[SignalSlow_StoTrend]=SlowStochasticWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);

[~,~,~,SlowPercentK,SlowPercentD]=SlowStochasticIndicator(HIGH,LOW,CLOSE);
SlowPercentK=table2array(SlowPercentK);
SlowPercentD=table2array(SlowPercentD);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);
[~,col_SlowPercentK]=size(SlowPercentK);
% Model

SignalSlow_StoTrend=zeros(row,col_MovRavand*col_SlowPercentK);
signal_idx=1;
for h=1:col_MovRavand
    for j=1:col_SlowPercentK
        for i=2:row
            if MovTrend(i,h)==1 && SlowPercentK(i,j)<20 && SlowPercentD(i,j)<20 && SlowPercentK(i-1,j)<SlowPercentD(i-1,j) && SlowPercentK(i,j)>SlowPercentD(i,j)
                SignalSlow_StoTrend(i,signal_idx)=1;
            elseif SlowPercentK(i,j)>80 && SlowPercentD(i,j)>80 && SlowPercentK(i-1,j)>SlowPercentD(i-1,j) && SlowPercentK(i,j)< SlowPercentD(i,j)
                SignalSlow_StoTrend(i,signal_idx)=0;
            else
                SignalSlow_StoTrend(i,signal_idx)=SignalSlow_StoTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalSlow_StoTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalSlow_StoTrend%d',x),1:num_cols,'UniformOutput',false);
SignalSlow_StoTrend=array2table(SignalSlow_StoTrend,"VariableNames",NameIndicator);

end

%%
function[SignalFastStoTrend]=FastStochasticWithTrendIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);

[~,FastPercentK,FastPercentD,~,~]=FastStochasticIndicator(HIGH,LOW,CLOSE);
FastPercentK=table2array(FastPercentK);
FastPercentD=table2array(FastPercentD);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);
[~,col_MovRavand]=size(MovTrend);
[~,col_FastPercentK]=size(FastPercentK);
% Model

SignalFastStoTrend=zeros(row,col_MovRavand*col_FastPercentK);
signal_idx=1;
for h=1:col_MovRavand
    for j=1:col_FastPercentK
        for i=2:row
            if MovTrend(i,h)==1 && FastPercentK(i,j)<20 && FastPercentD(i,j)<20 && FastPercentK(i-1,j)<FastPercentD(i-1,j) && FastPercentK(i,j)>FastPercentD(i,j)
                SignalFastStoTrend(i,signal_idx)=1;
            elseif FastPercentK(i,j)>80 && FastPercentD(i,j)>80 && FastPercentK(i-1,j)>FastPercentD(i-1,j) && FastPercentK(i,j)< FastPercentD(i,j)
                SignalFastStoTrend(i,signal_idx)=0;
            else
                SignalFastStoTrend(i,signal_idx)=SignalFastStoTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalFastStoTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalFastStoTrend%d',x),1:num_cols,'UniformOutput',false);
SignalFastStoTrend=array2table(SignalFastStoTrend,"VariableNames",NameIndicator);

end

%%
function[SignalSlowStoMFI]=SlowStochasticWithMfiIndicator(HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);

[~,~,~,SlowPercentK,SlowPercentD]=SlowStochasticIndicator(HIGH,LOW,CLOSE);
SlowPercentK=table2array(SlowPercentK);
SlowPercentD=table2array(SlowPercentD);

[~,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL);
MFI=table2array(MFI);
% Model
[~,col_SlowPercentK]=size(SlowPercentK);
[~,col_MFI]=size(MFI);
SignalSlowStoMFI=zeros(row,col_MFI*col_SlowPercentK);
signal_idx=1;
for h=1:col_MFI
    for j=1:col_SlowPercentK
        for i=2:row
            if MFI(i,h)>MFI(i-1,h) && SlowPercentK(i,j)<20 && SlowPercentD(i,j)<20 && SlowPercentK(i-1,j)<SlowPercentD(i-1,j) && SlowPercentK(i,j)>SlowPercentD(i,j)
                SignalSlowStoMFI(i,signal_idx)=1;
            elseif SlowPercentK(i,j)>80 && SlowPercentD(i,j)>80 && SlowPercentK(i-1,j)>SlowPercentD(i-1,j) && SlowPercentK(i,j)< SlowPercentD(i,j) && MFI(i,h)>MFI(i-1,h)
                SignalSlowStoMFI(i,signal_idx)=0;
            else
                SignalSlowStoMFI(i,signal_idx)=SignalSlowStoMFI(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalSlowStoMFI);
NameIndicator=arrayfun(@(x) sprintf('SignalSlowStoMFI%d',x),1:num_cols,'UniformOutput',false);
SignalSlowStoMFI=array2table(SignalSlowStoMFI,"VariableNames",NameIndicator);

end

%%
function[SignalFastStoMFI]=FastStochasticWithMfiIndicator(HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);

[~,FastPercentK,FastPercentD,~,~]=FastStochasticIndicator(HIGH,LOW,CLOSE);
FastPercentK=table2array(FastPercentK);
FastPercentD=table2array(FastPercentD);

[~,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL);
MFI=table2array(MFI);
% Model
[~,col_FastPercentK]=size(FastPercentK);
[~,col_MFI]=size(MFI);
SignalFastStoMFI=zeros(row,col_MFI*col_FastPercentK);
signal_idx=1;
for h=1:col_MFI
    for j=1:col_FastPercentK
        for i=2:row
            if MFI(i,h)>MFI(i-1,h) && FastPercentK(i,j)<20 && FastPercentD(i,j)<20 && FastPercentK(i-1,j)<FastPercentD(i-1,j) && FastPercentK(i,j)>FastPercentD(i,j)
                SignalFastStoMFI(i,signal_idx)=1;
            elseif FastPercentK(i,j)>80 && FastPercentD(i,j)>80 && FastPercentK(i-1,j)>FastPercentD(i-1,j) && FastPercentK(i,j)< FastPercentD(i,j) && MFI(i,h)>MFI(i-1,h)
                SignalFastStoMFI(i,signal_idx)=0;
            else
                SignalFastStoMFI(i,signal_idx)=SignalFastStoMFI(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalFastStoMFI);
NameIndicator=arrayfun(@(x) sprintf('SignalFastStoMFI%d',x),1:num_cols,'UniformOutput',false);
SignalFastStoMFI=array2table(SignalFastStoMFI,"VariableNames",NameIndicator);

end

%%
function[SignalTwoMovavgExponentialRsi]=TwoMovavgExponentialWithRsiIndicator(CLOSE)

[row,~]=size(CLOSE);
[SignaltwoMovExponential]=TwoMovavgExponentialindicator(CLOSE);
SignaltwoMovExponential=table2array(SignaltwoMovExponential);

[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);

% Model
[~,col_Rsi]=size(Rsi);
[~,col_buy_signal_two_mo_exponential]=size(SignaltwoMovExponential);
SignalTwoMovavgExponentialRsi=zeros(row,col_Rsi*col_buy_signal_two_mo_exponential);
signal_idx=1;
for j=1:col_buy_signal_two_mo_exponential
    for h=1:col_Rsi
        for i=2:row
            if SignaltwoMovExponential(i,j)==1 && Rsi(i-1,h) < 50 && Rsi(i,h) > 50
                SignalTwoMovavgExponentialRsi(i,signal_idx)=1;
            elseif SignaltwoMovExponential(i,j)==0 && Rsi(i-1,h) > 50 && Rsi(i,h) < 50
                SignalTwoMovavgExponentialRsi(i,signal_idx)=0;
            else
                SignalTwoMovavgExponentialRsi(i,signal_idx)=SignalTwoMovavgExponentialRsi(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalTwoMovavgExponentialRsi);
NameIndicator=arrayfun(@(x) sprintf('SignalTwoMovavgExponentialRsi%d',x),1:num_cols,'UniformOutput',false);
SignalTwoMovavgExponentialRsi=array2table(SignalTwoMovavgExponentialRsi,"VariableNames",NameIndicator);

end

%%
function[Signal2Rsi]=TwoRsiIndicator(CLOSE)

[row,~]=size(CLOSE);
[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);
% Model

Signal2Rsi=[];
[~,col_Rsi]=size(Rsi);
signal_idx=1;
for i=1:(col_Rsi-1)
    for j=i+1:col_Rsi
        Signal2Rsi(1,signal_idx)=0;
        for g=2:row
            if Rsi(g,j) >50 && Rsi(g-1,i) <30 && Rsi(g,i)>30
                Signal2Rsi(g,signal_idx)=1;
            elseif Rsi(g-1,i)>70 && Rsi(g,i)<70
                Signal2Rsi(g,signal_idx)=0;
            else 
                Signal2Rsi(g,signal_idx)=Signal2Rsi(g-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(Signal2Rsi);
NameIndicator=arrayfun(@(x) sprintf('Signal2Rsi%d',x),1:num_cols,'UniformOutput',false);
Signal2Rsi=array2table(Signal2Rsi,"VariableNames",NameIndicator);

end

%%
function[SignalMacdRsiEma]=RsiMacdEmaIndicator(CLOSE)

[row,~]=size(CLOSE);
[~,MACD,~,LineSignal]=Macdindicator(CLOSE);
MACD=table2array(MACD);
LineSignal=table2array(LineSignal);

[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);
% Model

mo_exponential13=movavg(CLOSE,'exponential',13);
mo_exponential21=movavg(CLOSE,'exponential',21);

[~,col_MACD]=size(MACD);
[~,col_Rsi]=size(Rsi);
SignalMacdRsiEma=zeros(row,col_MACD*col_Rsi);

signal_idx=1;
for j=1:col_MACD
    for g=1:col_Rsi
        for i=2:row
            if MACD(i,j)<=0 && MACD(i,j)>LineSignal(i,j) && Rsi(i,g)<30
                SignalMacdRsiEma(i,signal_idx)=1;
            elseif mo_exponential13(i,1)<mo_exponential21(i,1) && Rsi(i,g)>60
                SignalMacdRsiEma(i,signal_idx)=0;
            else
                SignalMacdRsiEma(i,signal_idx)=SignalMacdRsiEma(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalMacdRsiEma);
NameIndicator=arrayfun(@(x) sprintf('SignalMacdRsiEma%d',x),1:num_cols,'UniformOutput',false);
SignalMacdRsiEma=array2table(SignalMacdRsiEma,"VariableNames",NameIndicator);

end

%%
function[SignalMfiWiliams]=MfiWithRWilliamsIndicator(HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL);
MFI=table2array(MFI);

[~,williams]=RWilliamsIndicator(HIGH,LOW,CLOSE);
williams=table2array(williams);
% Model

[~,col_MFI]=size(MFI);
[~,col_williams]=size(williams);
SignalMfiWiliams=zeros(row,col_MFI*col_williams);
signal_idx=1;
for j=1:col_MFI
    for i=1:col_williams
        for g=2:row
            if williams(g-1,i)<-80 && williams(g,i)>-80 && MFI(g,j)>MFI(g-1,j)
                SignalMfiWiliams(g,signal_idx)=1;
            elseif williams(g-1,i)>-20 && williams(g,i)<-20 && MFI(g,j)>MFI(g-1,j)
                SignalMfiWiliams(g,signal_idx)=0;
            else
                SignalMfiWiliams(g,signal_idx)=SignalMfiWiliams(g-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalMfiWiliams);
NameIndicator=arrayfun(@(x) sprintf('SignalMfiWiliams%d',x),1:num_cols,'UniformOutput',false);
SignalMfiWiliams=array2table(SignalMfiWiliams,"VariableNames",NameIndicator);

end

%%
function[SignalMfiRsi]=MfiWithRsiIndicator(HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,MFI]=MfiIndicator(HIGH,LOW,CLOSE,VOL);
MFI=table2array(MFI);

[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);
% Model

[~,col_MFI]=size(MFI);
[~,col_Rsi]=size(Rsi);
SignalMfiRsi=zeros(row,col_MFI*col_Rsi);
signal_idx=1;
for j=1:col_MFI
    for i=1:col_Rsi
        for g=2:row
            if Rsi(g-1,i)<20 && Rsi(g,i)>20 && MFI(g,j)>MFI(g-1,j)
                SignalMfiRsi(g,signal_idx)=1;
            elseif Rsi(g-1,i)>80 && Rsi(g,i)<80 && MFI(g,j)>MFI(g-1,j)
                SignalMfiRsi(g,signal_idx)=0;
            else
                SignalMfiRsi(g,signal_idx)=SignalMfiRsi(g-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalMfiRsi);
NameIndicator=arrayfun(@(x) sprintf('SignalMfiRsi%d',x),1:num_cols,'UniformOutput',false);
SignalMfiRsi=array2table(SignalMfiRsi,"VariableNames",NameIndicator);

end

%%
function[SignalMacdRsi]=RsiMacdIndicator(CLOSE)

[row,~]=size(CLOSE);
[~,MACD,~,LineSignal]=Macdindicator(CLOSE);
MACD=table2array(MACD);
LineSignal=table2array(LineSignal);

[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);
% Model

[~,Col_MACD]=size(MACD);
[~,Col_Rsi]=size(Rsi);
SignalMacdRsi=zeros(row,Col_MACD*Col_Rsi);
signal_idx=1;
for j=1:Col_MACD
    for g=1:Col_Rsi
        for i=2:row
            if MACD(i,j)<=0 && MACD(i,j)>LineSignal(i,j) && Rsi(i,g)<30
                SignalMacdRsi(i,signal_idx)=1;
            elseif MACD(i,j)>0 && MACD(i,j) < LineSignal(i,j)
                SignalMacdRsi(i,signal_idx)=0;
            else
                SignalMacdRsi(i,signal_idx)=SignalMacdRsi(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalMacdRsi);
NameIndicator=arrayfun(@(x) sprintf('SignalMacdRsi%d',x),1:num_cols,'UniformOutput',false);
SignalMacdRsi=array2table(SignalMacdRsi,"VariableNames",NameIndicator);

end

%%
function[SignalRsiTrend]=RsiWithTerndIndIndicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
[~,Rsi]=RsiIndicator(CLOSE);
Rsi=table2array(Rsi);

MovTrend=Terndindicator(DATE,OPEN,HIGH,LOW,CLOSE,VOL);

% Model
[~,Col_Rsi]=size(Rsi);
[~,col_MovRavand]=size(MovTrend);
SignalRsiTrend=zeros(row,Col_Rsi*col_MovRavand);
signal_idx=1;
for h=1:col_MovRavand
    for j=1:Col_Rsi
        for i=2:row
            if Rsi(i-1,j)<30 && Rsi(i,j)>30 && MovTrend(i,h)==1
                SignalRsiTrend(i,signal_idx)=1;
            elseif Rsi(i-1,j)>70 && Rsi(i,j)<70
                SignalRsiTrend(i,signal_idx)=0;
            else
                SignalRsiTrend(i,signal_idx)=SignalRsiTrend(i-1,signal_idx);
            end
        end
        signal_idx=signal_idx+1;
    end
end
[~, num_cols] = size(SignalRsiTrend);
NameIndicator=arrayfun(@(x) sprintf('SignalRsiTrend%d',x),1:num_cols,'UniformOutput',false);
SignalRsiTrend=array2table(SignalRsiTrend,"VariableNames",NameIndicator);

end

%%
function[OtherTechnicalIndicators]=OtherIndicators(DATE,OPEN,HIGH,LOW,CLOSE,VOL)

[row,~]=size(CLOSE);
% roc= indicators(price,'roc',period);
period=[8 14 21 34 55 63 89];
ROC=[];
idx=1;
for i=period
    ROC(:,idx)= indicators(CLOSE,'roc',i);
    idx=idx+1;
end
% ROC=fillmissing(ROC,"linear",1,"EndValues","nearest");
[~, num_cols] = size(ROC);
NameIndicator=arrayfun(@(x) sprintf('ROC%d',x),1:num_cols,'UniformOutput',false);
ROC=array2table(ROC,"VariableNames",NameIndicator);

%% KDJ
% kd=[14 3;20 5;21 8;50 14];
% % [fpctk,fpctd,jline]  = indicators([HIGH,LOW,CLOSE]      ,'kdj'    ,k,d);
% fpctk_kdi=zeros(row,size(kd,1));
% fpctd_kdi=zeros(row,size(kd,1));
% jline_kdi=zeros(row,size(kd,1));
% idx=1;
% for i=1:size(kd,1)
%     kdj=indicators([HIGH,LOW,CLOSE],'kdj',kd(i,1),kd(i,2));
%     fpctk_kdi(:,idx)=kdj(:,1);
%     fpctd_kdi(:,idx)=kdj(:,2);
%     jline_kdi(:,idx)=kdj(:,3);
%     idx=idx+1;
% end
% 
% [~, num_cols] = size(fpctk_kdi);
% NameIndicator=arrayfun(@(x) sprintf('fpctk_kdi%d',x),1:num_cols,'UniformOutput',false);
% fpctk_kdi=array2table(fpctk_kdi,"VariableNames",NameIndicator);
% 
% NameIndicator=arrayfun(@(x) sprintf('fpctd_kdi%d',x),1:num_cols,'UniformOutput',false);
% fpctd_kdi=array2table(fpctd_kdi,"VariableNames",NameIndicator);
% 
% NameIndicator=arrayfun(@(x) sprintf('jline_kdi%d',x),1:num_cols,'UniformOutput',false);
% jline_kdi=array2table(jline_kdi,"VariableNames",NameIndicator);

%% aroon
% [dn,up,os]=indicators([hi,lo],'aroon',period)
period=[7 14 21 30 63 100 34 55 89];
p=size(period,2);
dn_aroon=zeros(row,p);
up_aroon=zeros(row,p);
os_aroon=zeros(row,p);
idx=1;
for i=period
    aroon=indicators([HIGH,LOW],'aroon',i);
    dn_aroon(:,idx)=aroon(:,1);
    up_aroon(:,idx)=aroon(:,2);
    os_aroon(:,idx)=aroon(:,3);
    idx=idx+1;
end
% dn_aroon=fillmissing(dn_aroon,"linear",1,"EndValues","nearest");
% up_aroon=fillmissing(up_aroon,"linear",1,"EndValues","nearest");
% os_aroon=fillmissing(os_aroon,"linear",1,"EndValues","nearest");
[~, num_cols] = size(dn_aroon);
NameIndicator=arrayfun(@(x) sprintf('dn_aroon%d',x),1:num_cols,'UniformOutput',false);
dn_aroon=array2table(dn_aroon,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('up_aroon%d',x),1:num_cols,'UniformOutput',false);
up_aroon=array2table(up_aroon,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('os_aroon%d',x),1:num_cols,'UniformOutput',false);
os_aroon=array2table(os_aroon,"VariableNames",NameIndicator);

%% 

period = [25,13;  
          12,6; 
          20,10; 
          34,13; 
          50,25; 
          8,3;
          30,15;   
          40,20; 
          18,9; 
          60,30; 
          10,5; 
          15,7];

p=size(period,1);
tsi=zeros(row,p);
for i=1:p
    r=period(i,1);
    s=period(i,2);
    ts=indicators(CLOSE,'tsi',r,s);
    tsi(:,i)=ts;
end
% tsi=fillmissing(tsi,"linear",1,"EndValues","nearest");
[~, num_cols] = size(tsi);
NameIndicator=arrayfun(@(x) sprintf('tsi%d',x),1:num_cols,'UniformOutput',false);
tsi=array2table(tsi,"VariableNames",NameIndicator);

%%
% [pdi,mdi,adx]= indicators([hi,lo,cl],'adx',period)
% HIGH,LOW,CLOSE
period=[7, 9, 13, 14, 21, 25, 28, 34, 50];
p=size(period,2);
pdi_adx=zeros(row,p);
mdi_adx=zeros(row,p);
adx=zeros(row,p);
idx=1;
for i=period
    adx1=indicators([HIGH,LOW,CLOSE],'adx',i);
    pdi_adx(:,idx)=adx1(:,1);
    mdi_adx(:,idx)=adx1(:,2);
    adx(:,idx)=adx1(:,3);
    idx=idx+1;
end

[~, num_cols] = size(pdi_adx);
NameIndicator=arrayfun(@(x) sprintf('pdi_adx%d',x),1:num_cols,'UniformOutput',false);
pdi_adx=array2table(pdi_adx,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('mdi_adx%d',x),1:num_cols,'UniformOutput',false);
mdi_adx=array2table(mdi_adx,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('adx%d',x),1:num_cols,'UniformOutput',false);
adx=array2table(adx,"VariableNames",NameIndicator);

%%
% t3= indicators(price,'t3',period,volfact)
period1=[5 0.7;8 0.6;10 0.5;14 0.7;5 0.618;8 0.5;13 0.7;3 0.4;7 0.6];


p=size(period1,1);
t3=zeros(row,p);
idx=1;
for i=1:p
    period=period1(i,1);
    volfact=period1(i,2);
    t= indicators(CLOSE,'t3',period,volfact);
    t3(:,i)=t;
    idx=idx+1;
end
% t3=fillmissing(t3,"linear",1,"EndValues","nearest");
[~, num_cols] = size(t3);
NameIndicator=arrayfun(@(x) sprintf('t3%d',x),1:num_cols,'UniformOutput',false);
t3=array2table(t3,"VariableNames",NameIndicator);

%% OBV
obv= indicators([CLOSE,VOL],'obv');

[~, num_cols] = size(obv);
NameIndicator=arrayfun(@(x) sprintf('obv%d',x),1:num_cols,'UniformOutput',false);
obv=array2table(obv,"VariableNames",NameIndicator);

%%
% cmf= indicators([hi,lo,cl,vo],'cmf',period)
period1=[21 14 63 8 34 55 89];
p=size(period1,2);
idx=1;
for i=period1
    period=i;
    cmf(:,idx)= indicators([HIGH,LOW,CLOSE,VOL],'cmf',period);
    idx=idx+1;
end

[~, num_cols] = size(cmf);
NameIndicator=arrayfun(@(x) sprintf('cmf%d',x),1:num_cols,'UniformOutput',false);
cmf=array2table(cmf,"VariableNames",NameIndicator);

%%
% force= indicators([cl,vo],'force',period)
period1=[21 14 63 8 34 55 89];
p=size(period1,2);
idx=1;
for i=period1
    period=i;
    force(:,idx)= indicators([CLOSE,VOL],'force',period);
    idx=idx+1;
end

[~, num_cols] = size(force);
NameIndicator=arrayfun(@(x) sprintf('force%d',x),1:num_cols,'UniformOutput',false);
force=array2table(force,"VariableNames",NameIndicator);


%%
% [middle,upper,lower]=indicators([hi,lo,cl],'keltner',emaper,atrmul,atrper);
period1 = [
    20, 2.0, 20;   % تنظیمات استاندارد (EMA20, ATR20, ضریب ۲)
    10, 1.5, 10;   % کوتاه‌مدت (واکنش سریع)
    50, 2.5, 20;   % بلندمدت (فیلتر نویز قوی)
    21, 2.0, 21;   % فیبوناچی (EMA21, ATR21)
    34, 1.618, 34; % نسبت طلایی فیبوناچی
    13, 3.0, 14    % ترکیب پرنوسان (اسکالپینگ)
];

p=size(period1,1);
for i=1:p
    emaper=period1(i,1);
    atrmul=period1(i,2);
    atrper=period1(i,3);
    keltner=indicators([HIGH,LOW,CLOSE],'keltner',emaper,atrmul,atrper);
    middle_keltner(:,i)=keltner(:,1);
    upper_keltner(:,i)=keltner(:,2);
    lower_keltner(:,i)=keltner(:,3);
end

[~, num_cols] = size(middle_keltner);
NameIndicator=arrayfun(@(x) sprintf('middle_keltner%d',x),1:num_cols,'UniformOutput',false);
middle_keltner=array2table(middle_keltner,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('upper_keltner%d',x),1:num_cols,'UniformOutput',false);
upper_keltner=array2table(upper_keltner,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('lower_keltner%d',x),1:num_cols,'UniformOutput',false);
lower_keltner=array2table(lower_keltner,"VariableNames",NameIndicator);

%%
% atr= indicators([hi,lo,cl],'atr',period)
period1=[21 14 63 8 34 55 89];
% p=size(period1,2);
idx=1;
for i=period1
    period=i;
    atr(:,idx)= indicators([HIGH,LOW,CLOSE],'atr',period);
    idx=idx+1;
end

[~, num_cols] = size(atr);
NameIndicator=arrayfun(@(x) sprintf('atr%d',x),1:num_cols,'UniformOutput',false);
atr=array2table(atr,"VariableNames",NameIndicator);

%%
% vr= indicators([hi,lo,cl],'vr',period)
period1=[21 14 63 8 34 55 89];
idx=1;
for i=period1
    period=i;
    VolatilityRatio(:,idx)= indicators([HIGH,LOW,CLOSE],'vr',period);
    idx=idx+1;
end
[~, num_cols] = size(atr);
NameIndicator=arrayfun(@(x) sprintf('VolatilityRatio%d',x),1:num_cols,'UniformOutput',false);
VolatilityRatio=array2table(VolatilityRatio,"VariableNames",NameIndicator);

%% zigzag
% [index,value]= indicators(price,'zigzag' ,moveper)
% moveper1=[0.005:0.005:0.07];
% idx=1;
% for i=moveper1
%     moveper=i;
%     zigzag= indicators(CLOSE,'zigzag',moveper);
%     index_zigzag(:,idx)=zigzag(:,1);
%     value_zigzag(:,idx)=zigzag(:,2);
%     idx=idx+1;
% end

[PP_pivot, R1_pivot, R2_pivot, R3_pivot, S1_pivot, S2_pivot, S3_pivot] = pivot_points(HIGH, LOW, CLOSE);

[~, num_cols] = size(PP_pivot);
NameIndicator=arrayfun(@(x) sprintf('PP_pivot%d',x),1:num_cols,'UniformOutput',false);
PP_pivot=array2table(PP_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('R1_pivot%d',x),1:num_cols,'UniformOutput',false);
R1_pivot=array2table(R1_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('R2_pivot%d',x),1:num_cols,'UniformOutput',false);
R2_pivot=array2table(R2_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('R3_pivot%d',x),1:num_cols,'UniformOutput',false);
R3_pivot=array2table(R3_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('S1_pivot%d',x),1:num_cols,'UniformOutput',false);
S1_pivot=array2table(S1_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('S2_pivot%d',x),1:num_cols,'UniformOutput',false);
S2_pivot=array2table(S2_pivot,"VariableNames",NameIndicator);

NameIndicator=arrayfun(@(x) sprintf('S3_pivot%d',x),1:num_cols,'UniformOutput',false);
S3_pivot=array2table(S3_pivot,"VariableNames",NameIndicator);

OtherTechnicalIndicators=[ROC,dn_aroon,up_aroon,os_aroon,tsi,pdi_adx,mdi_adx,adx,...
    t3,obv,cmf,force,middle_keltner,upper_keltner,lower_keltner,atr,VolatilityRatio,PP_pivot,R1_pivot,R2_pivot,R3_pivot,S1_pivot,S2_pivot,S3_pivot];

% OtherTechnicalIndicators=[ROC,fpctk_kdi,fpctd_kdi,jline_kdi,dn_aroon,up_aroon,os_aroon,tsi,pdi_adx,mdi_adx,adx,...
%     t3,obv,cmf,force,middle_keltner,upper_keltner,lower_keltner,atr,VolatilityRatio,PP_pivot,R1_pivot,R2_pivot,R3_pivot,S1_pivot,S2_pivot,S3_pivot];
end

%%
function [PP, R1, R2, R3, S1, S2, S3] = pivot_points(high, low, close)
% محاسبه پیوت پوینت و سطوح حمایت/مقاومت
% ورودی:
%   high: بردار قیمتهای High
%   low: بردار قیمتهای Low
%   close: بردار قیمتهای Close
% خروجی:
%   PP: بردار پیوت پوینت
%   R1, R2, R3: سطوح مقاومت
%   S1, S2, S3: سطوح حمایت

n = length(high);
PP = zeros(n, 1);
R1 = zeros(n, 1);
R2 = zeros(n, 1);
R3 = zeros(n, 1);
S1 = zeros(n, 1);
S2 = zeros(n, 1);
S3 = zeros(n, 1);

for i = 2:n
    % محاسبه پیوت پوینت برای روز جاری بر اساس روز قبل
    PP(i) = (high(i-1) + low(i-1) + close(i-1)) / 3;

    % محاسبه سطوح مقاومت (R1, R2, R3)
    R1(i) = 2 * PP(i) - low(i-1);
    R2(i) = PP(i) + (high(i-1) - low(i-1));
    R3(i) = high(i-1) + 2 * (PP(i) - low(i-1));

    % محاسبه سطوح حمایت (S1, S2, S3)
    S1(i) = 2 * PP(i) - high(i-1);
    S2(i) = PP(i) - (high(i-1) - low(i-1));
    S3(i) = low(i-1) - 2 * (high(i-1) - PP(i));
end

% مقداردهی اولیه روز اول به NaN (چون داده قبلی وجود ندارد)
PP(1) = NaN;
R1(1) = NaN;
R2(1) = NaN;
R3(1) = NaN;
S1(1) = NaN;
S2(1) = NaN;
S3(1) = NaN;
end



















%%






