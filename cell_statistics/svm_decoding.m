function output=svm_decoding(xtrain,ytrain,xtest)

options=statset('MaxIter',1e6);
%svmStruct = svmtrain(xtrain,ytrain,'kernel_function','mlp','options',options);
svmStruct = svmtrain(xtrain,ytrain,'options',options);
output = svmclassify(svmStruct,xtest);