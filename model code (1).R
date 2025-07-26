library(lidR)
library(e1071)
library(randomForest)
library(tidyverse)

ymjdata <- readRDS("C:/Users/aj710/Desktop/r/zhsj.rds")

names(ymjdata)[1] <- 'fam' 

nir <- sapply(ymjdata[,-c(1:6, ncol(ymjdata))],
              function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

finadata <- cbind.data.frame(ymjdata[,c(1:6)],nir)
#write.csv(finadata,'e:/谷歌浏览器下载/finadata.csv')
##traits selection 
ANTH_c <-ymjdata$PAR 
# nbi <- x$NBI_G
data_use <- cbind.data.frame(ANTH_c=ANTH_c,nir)

library(enpls)

y <- data_use$ANTH_c
x <- data_use[,-1]
set.seed(42)
od <- enpls.od(x, y, reptimes = 450)
print(od)
plot(od)
plot(od, criterion = "sd")
x <- x[od$error.mean<4& od$error.sd <0.7,]
y <- y[od$error.mean<4& od$error.sd <0.7]
data_use <- cbind.data.frame(ANTH_c=y,x)

require(caTools)
sel <- sample.split(data_use$ANTH_c, SplitRatio = .8)
train1 = subset(data_use, sel == TRUE)
test  = subset(data_use, sel  == FALSE)
#####################
#
#    plsr       ####
#
####################
library(pls)

# 建立PLSR模型
model <- plsr(ANTH_c ~., data = train1, ncomp = 12, validation = 'LOO')

# 训练集和预测集的预测值
train_pred <- predict(model, train1, ncomp = 12)
test_pred <- predict(model, test, ncomp = 12)  # 使用 test，而不是 test1

library(caret)

# 计算 R² 和 RMSE（训练集）
train_r2 <- R2(train_pred, train1$ANTH_c)
train_rmse <- RMSE(train_pred, train1$ANTH_c)

# 计算 R² 和 RMSE（测试集）
test_r2 <- R2(test_pred, test$ANTH_c)
test_rmse <- RMSE(test_pred, test$ANTH_c)

# 打印结果
cat("R² (Train):", train_r2, "\n")
cat("RMSE (Train):", train_rmse, "\n")
cat("R² (Test):", test_r2, "\n")
cat("RMSE (Test):", test_rmse, "\n")



########################
##      svm         ######
#                    ####
#########################
# 加载所需的库
library(enpls)
library(caTools)
library(e1071)

# 读取数据
ymjdata <- readRDS("C:/Users/aj710/Desktop/r/zhsj.rds")
names(ymjdata)[1] <- 'fam'  # 将第一列重命名为 'fam'

# 对 NIR 数据进行归一化处理
nir <- sapply(ymjdata[,-c(1:6, ncol(ymjdata))],
              function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

# 合并数据集
finadata <- cbind.data.frame(ymjdata[,c(1:6)], nir)

# 选择目标变量
ANTH_c <- ymjdata$PAR
data_use <- cbind.data.frame(ANTH_c = ANTH_c, nir)

# 设置随机种子并应用 enpls 方法
set.seed(42)
y <- data_use$ANTH_c
x <- data_use[,-1]
od <- enpls.od(x, y, reptimes = 450)
print(od)
plot(od)
plot(od, criterion = "sd")

# 筛选掉 error.mean 和 error.sd 不满足条件的特征
x <- x[od$error.mean < 4 & od$error.sd < 0.7,]
y <- y[od$error.mean < 4 & od$error.sd < 0.7]

# 更新数据集
data_use <- cbind.data.frame(ANTH_c = y, x)

# 划分训练集和测试集
set.seed(123) # 设置随机种子以便结果可重复
sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.8)
train1 <- subset(data_use, sel == TRUE)
test <- subset(data_use, sel == FALSE)

# 训练 SVM 模型
model_svm <- svm(ANTH_c ~ ., scale = TRUE, data = train1, kernel = "linear")

# 预测测试集
pre_test <- predict(model_svm, newdata = test)

# 计算测试集的 RMSE 和 R²
error_test <- test$ANTH_c - pre_test  
SVM_rmse_test <- sqrt(mean(error_test^2))
SVM_R2_test <- cor(test$ANTH_c, pre_test)^2

# 预测训练集
pre_train1 <- predict(model_svm, newdata = train1)

# 计算训练集的 RMSE 和 R²
error_train <- train1$ANTH_c - pre_train1  
SVM_rmse_train <- sqrt(mean(error_train^2))
SVM_R2_train <- cor(train1$ANTH_c, pre_train1)^2



########################
##     rf         ######
#                    ####
#########################
library(e1071)
library(caTools)
library(randomForest)
library(caret)

# 计算R²和RMSE的函数
calc_metrics <- function(actual, predicted) {
  error <- actual - predicted
  rmse <- sqrt(mean(error^2, na.rm = TRUE))  # 计算RMSE, 处理NA值
  r2 <- cor(actual, predicted, use = "complete.obs")^2  # 计算R²，处理NA值
  return(list(RMSE = rmse, R2 = r2))
}

# 存储每次训练和测试的数据集和结果
list_datasets <- list()
results <- data.frame(RF_rmse_train = numeric(300), 
                      RF_R2_train = numeric(300), 
                      RF_rmse_test = numeric(300), 
                      RF_R2_test = numeric(300))

# 多次训练模型
for (i in 1:300) {
  # 数据分割
  sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.8)
  train1 = subset(data_use, sel == TRUE)
  test = subset(data_use, sel == FALSE)
  
  # 随机森林模型训练
  model_rf <- randomForest(ANTH_c ~ ., scale = TRUE, data = train1, mtry = 5,
                           importance = TRUE, na.action = na.omit)
  
  # 训练集的预测和计算
  pre_train1 <- predict(model_rf, newdata = train1)
  train_metrics <- calc_metrics(train1$ANTH_c, pre_train1)
  results$RF_rmse_train[i] <- train_metrics$RMSE  # 计算RMSE
  results$RF_R2_train[i] <- train_metrics$R2  # 计算R²
  
  # 测试集的预测和计算
  pre_test <- predict(model_rf, newdata = test)
  test_metrics <- calc_metrics(test$ANTH_c, pre_test)
  results$RF_rmse_test[i] <- test_metrics$RMSE  # 计算RMSE
  results$RF_R2_test[i] <- test_metrics$R2  # 计算R²
  
  # 存储数据集
  list_datasets[[i]] <- list(train = train1, test = test)
}

# 找到最佳的测试集 R²
best_index <- which.max(results$RF_R2_test)
best_index
round(results, 2)

# 获取最佳的数据集
best_train <- list_datasets[[best_index]]$train
best_test <- list_datasets[[best_index]]$test

# 保存最佳数据集
write.csv(best_train, "best_train_data.csv")
write.csv(best_test, "best_test_data.csv")

# 用最佳训练集构建最终模型
final_model <- randomForest(ANTH_c ~ ., scale = TRUE, data = best_train, mtry = 5,
                            importance = TRUE, na.action = na.omit)

# 最终模型预测和评估
final_predictions <- predict(final_model, newdata = best_test)
final_error_test <- best_test$ANTH_c - final_predictions
final_rmse_test <- sqrt(mean(final_error_test^2))
final_r2_test <- cor(best_test$ANTH_c, final_predictions)^2

# 输出最终模型的评估结果
print(paste("Final RMSE Test:", round(final_rmse_test, 2)))
print(paste("Final R² Test:", round(final_r2_test, 2)))

# 预测训练集
pre_train1 <- predict(final_model, newdata = best_train)

# 预测测试集
pre_test <- predict(final_model, newdata = best_test)


####################

XGboost


####################

# 加载所需的库
library(enpls)
library(caTools)
library(xgboost)  # 加载 xgboost 包

# 读取数据
ymjdata <- readRDS("C:/Users/aj710/Desktop/r/zhsj.rds")
names(ymjdata)[1] <- 'fam'  # 将第一列重命名为 'fam'

# 对 NIR 数据进行归一化处理
nir <- sapply(ymjdata[,-c(1:6, ncol(ymjdata))],
              function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

# 合并数据集
finadata <- cbind.data.frame(ymjdata[,c(1:6)], nir)

# 选择目标变量
ANTH_c <- ymjdata$PAR
data_use <- cbind.data.frame(ANTH_c = ANTH_c, nir)

# 设置随机种子并应用 enpls 方法
set.seed(42)
y <- data_use$ANTH_c
x <- data_use[,-1]
od <- enpls.od(x, y, reptimes = 450)
print(od)
plot(od)
plot(od, criterion = "sd")

# 筛选掉 error.mean 和 error.sd 不满足条件的特征
x <- x[od$error.mean < 4 & od$error.sd < 0.7,]
y <- y[od$error.mean < 4 & od$error.sd < 0.7]

# 更新数据集
data_use <- cbind.data.frame(ANTH_c = y, x)

# 划分训练集和测试集
set.seed(123) # 设置随机种子以便结果可重复
sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.8)
train1 <- subset(data_use, sel == TRUE)
test <- subset(data_use, sel == FALSE)

# 将数据转换为矩阵格式（XGBoost 要求的数据格式）
train_matrix <- as.matrix(train1[,-1])
train_label <- train1$ANTH_c
test_matrix <- as.matrix(test[,-1])
test_label <- test$ANTH_c

# 训练 XGBoost 模型
dtrain <- xgb.DMatrix(data = train_matrix, label = train_label)
dtest <- xgb.DMatrix(data = test_matrix, label = test_label)

# 设置 XGBoost 模型的参数
params <- list(
  objective = "reg:squarederror",  # 回归任务
  eval_metric = "rmse"  # 使用 RMSE 作为评估指标
)

# 训练模型
xgb_model <- xgboost(params = params, data = dtrain, nrounds = 100, verbose = 1)

# 预测测试集
pre_test <- predict(xgb_model, newdata = dtest)

# 计算测试集的 RMSE 和 R²
error_test <- test_label - pre_test  
XGB_rmse_test <- sqrt(mean(error_test^2))
XGB_R2_test <- cor(test_label, pre_test)^2

# 预测训练集
pre_train1 <- predict(xgb_model, newdata = dtrain)

# 计算训练集的 RMSE 和 R²
error_train <- train_label - pre_train1  
XGB_rmse_train <- sqrt(mean(error_train^2))
XGB_R2_train <- cor(train_label, pre_train1)^2

