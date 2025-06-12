function results=BayesianOptimizationNeuralNetwork(CrossValidationData,ForTheModel)

numLayers=optimizableVariable('numLayers',[1 6],'Type','integer');
numHiddenUnits1=optimizableVariable('numHiddenUnits1',[8 512],'Type','integer');
numHiddenUnits2=optimizableVariable('numHiddenUnits2',[8 512],'Type','integer');
numHiddenUnits3=optimizableVariable('numHiddenUnits3',[8 512],'Type','integer');
numHiddenUnits4=optimizableVariable('numHiddenUnits4',[8 512],'Type','integer');
numHiddenUnits5=optimizableVariable('numHiddenUnits5',[8 512],'Type','integer');

numHiddenUnitsFC1=optimizableVariable('numHiddenUnitsFC1',[8 512],'Type','integer');

probability1=optimizableVariable('probability1',[0.1 0.5],'Type','real');
probability2=optimizableVariable('probability2',[0.1 0.5],'Type','real');
probability3=optimizableVariable('probability3',[0.1 0.5],'Type','real');
probability4=optimizableVariable('probability4',[0.1 0.5],'Type','real');
probability5=optimizableVariable('probability5',[0.1 0.5],'Type','real');
probability6=optimizableVariable('probability6',[0.1 0.5],'Type','real');

L2Value=optimizableVariable('L2Value',[0.0001 0.1]);
LearnRate=optimizableVariable('LearnRate',[0.001 0.01]);
numEpochs=optimizableVariable('numEpochs',[20 300],'Type','integer');
MiniBatchSize=optimizableVariable('MiniBatchSize',[16 256],'Type','integer');
gradientThreshold=optimizableVariable('gradientThreshold',[0.1 10],'Type','real');
%%

vars=[numLayers,numHiddenUnits1,numHiddenUnits2,numHiddenUnits3,numHiddenUnits4,numHiddenUnits5,...
    probability1,probability2,probability3,probability4,probability5,probability6,...
    numHiddenUnitsFC1,numEpochs,MiniBatchSize,gradientThreshold,L2Value,LearnRate];

fun=@(params)TrainNeuralNetwork(CrossValidationData,ForTheModel,params.numLayers,params.numHiddenUnits1,...
    params.numHiddenUnits2,params.numHiddenUnits3,params.numHiddenUnits4,params.numHiddenUnits5,...
    params.probability1,params.probability2,params.probability3,params.probability4,params.probability5,params.probability6,...
    params.numHiddenUnitsFC1,params.numEpochs,params.MiniBatchSize,params.gradientThreshold,params.L2Value,params.LearnRate);

results=bayesopt(fun,vars,"UseParallel",false,"Verbose",1,"MaxObjectiveEvaluations",3,"AcquisitionFunctionName","expected-improvement-per-second-plus");

end

