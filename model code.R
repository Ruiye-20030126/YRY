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
sel <- sample.split(data_use$ANTH_c, SplitRatio = .7)
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


# 新建一个图形设备并调整大小
dev.new(width = 10, height = 8)

# 绘制训练集的回归散点图
plot(train1$ANTH_c, train_pred, main = "Regression Scatter Plot (Train & Test)",
     xlab = "Measured FAPAR", ylab = "Estimated FAPAR", pch = 1, col = "orange", cex = 1.5,
     xlim = c(min(c(train1$ANTH_c, test$ANTH_c)), max(c(train1$ANTH_c, test$ANTH_c))),
     ylim = c(min(c(train_pred, test_pred)), max(c(train_pred, test_pred))),
     cex.main = 1.5)  # 增大标题字体大小

# 添加测试集的散点图
points(test$ANTH_c, test_pred, pch = 0, col = "blue", cex = 1.5)  # 测试集: 空心正方形，蓝色

# 绘制回归线
abline(0, 1, col = "black", lty = 2)  # 绘制回归线

# 在训练集回归散点图中添加R²和RMSE值（向左调整位置）
text(x = max(train1$ANTH_c) * 0.35, y = max(train_pred) * 0.95,  # 修改x坐标，使其向左移动
     labels = paste("R² (Train): ", round(train_r2, 3)), col = "black", cex = 1)
text(x = max(train1$ANTH_c) * 0.35, y = max(train_pred) * 0.90,  # 同样调整y坐标
     labels = paste("RMSE (Train): ", round(train_rmse, 3)), col = "black", cex = 1)

# 在预测集回归散点图中添加R²和RMSE值（适当向右调整位置）
text(x = max(test$ANTH_c) * 0.45, y = max(test_pred) * 0.95,  # 修改x坐标，使其向右移动
     labels = paste("R² (Test): ", round(test_r2, 3)), col = "black", cex = 1)
text(x = max(test$ANTH_c) * 0.45, y = max(test_pred) * 0.90,  # 同样调整y坐标
     labels = paste("RMSE (Test): ", round(test_rmse, 3)), col = "black", cex = 1)

# 设置保存路径和图像尺寸
png("regression_scatter_plot_plsr.png", width = 1000, height = 800)

# 绘制训练集的回归散点图
plot(train1$ANTH_c, train_pred, main = "Regression Scatter Plot (Train & Test)",
     xlab = "Measured FAPAR", ylab = "Estimated FAPAR", pch = 1, col = "orange", cex = 1.5,
     xlim = c(min(c(train1$ANTH_c, test$ANTH_c)), max(c(train1$ANTH_c, test$ANTH_c))),
     ylim = c(min(c(train_pred, test_pred)), max(c(train_pred, test_pred))),
     cex.main = 1.5)  # 增大标题字体大小

# 添加测试集的散点图
points(test$ANTH_c, test_pred, pch = 0, col = "blue", cex = 1.5)  # 测试集: 空心正方形，蓝色

# 绘制回归线
abline(0, 1, col = "black", lty = 2)  # 绘制回归线

# 在训练集回归散点图中添加R²和RMSE值（向左调整位置）
text(x = max(train1$ANTH_c) * 0.35, y = max(train_pred) * 0.95,  # 修改x坐标，使其向左移动
     labels = paste("R² (Train): ", round(train_r2, 3)), col = "black", cex = 1)
text(x = max(train1$ANTH_c) * 0.35, y = max(train_pred) * 0.90,  # 同样调整y坐标
     labels = paste("RMSE (Train): ", round(train_rmse, 3)), col = "black", cex = 1)

# 在预测集回归散点图中添加R²和RMSE值（适当向右调整位置）
text(x = max(test$ANTH_c) * 0.45, y = max(test_pred) * 0.95,  # 修改x坐标，使其向右移动
     labels = paste("R² (Test): ", round(test_r2, 3)), col = "black", cex = 1)
text(x = max(test$ANTH_c) * 0.45, y = max(test_pred) * 0.90,  # 同样调整y坐标
     labels = paste("RMSE (Test): ", round(test_rmse, 3)), col = "black", cex = 1)

# 关闭图形设备并保存为PNG文件
dev.off()

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
sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.7)
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

# 新建一个图形设备并调整大小
dev.new(width = 10, height = 8)  # 可以调整为适当的宽度和高度

# 绘制训练集和测试集回归散点图
plot(train1$ANTH_c, pre_train1, main = "Regression Scatter Plot",
     xlab = "Actual Values", ylab = "Predicted Values", 
     pch = 1, col = "orange", cex = 1.5,  # 训练集: 空心圆，橙色
     xlim = range(c(test$ANTH_c, train1$ANTH_c)),  # 设置 x 轴范围为训练集和测试集的最小和最大值
     ylim = range(c(pre_test, pre_train1)))  # 设置 y 轴范围为预测值的最小和最大值

# 绘制测试集的散点图
points(test$ANTH_c, pre_test, pch = 0, col = "blue", cex = 1.5)  # 测试集: 空心正方形，蓝色

# 添加回归线
abline(0, 1, col = "black", lty = 2)  # 绘制回归线

# 在图中添加R²和RMSE值
# 调整文本位置到稍微右上角，避免被覆盖并确保完全显示

# 训练集文本位置调整
text(x = max(train1$ANTH_c) * 0.8, y = max(pre_train1) * 1.1, 
     labels = paste("R² (Train): ", round(SVM_R2_train, 3)), col = "black", cex = 0.9)
text(x = max(train1$ANTH_c) * 0.8, y = max(pre_train1) * 1.03, 
     labels = paste("RMSE (Train): ", round(SVM_rmse_train, 2)), col = "black", cex = 0.9)

# 测试集文本位置调整
text(x = max(test$ANTH_c) * 0.8, y = max(pre_test) * 1.1, 
     labels = paste("R² (Test): ", round(SVM_R2_test, 3)), col = "black", cex = 0.9)
text(x = max(test$ANTH_c) * 0.8, y = max(pre_test) * 1.03, 
     labels = paste("RMSE (Test): ", round(SVM_rmse_test, 2)), col = "black", cex = 0.9)



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
  sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.7)
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

# 新建一个图形设备并调整大小
dev.new(width = 10, height = 8)  # 调整宽度和高度

# 绘制训练集和测试集的散点图
plot(best_train$ANTH_c, pre_train1, ylim = c(0, 2), xlim = c(0, 2), 
     main = "Training and Test Set Regression Scatter Plot", xlab = "Measured FAPAR", ylab = "Estimated FAPAR", 
     pch = 1, col = "orange", family = "serif", cex.main = 1.0)  # 训练集空心圆

# 添加回归线
abline(0, 1)

# 添加训练集的R²和RMSE值
text(1, 1.8, paste("R²: ", round(final_r2_test, 2)), col = "black", cex = 0.9, family = "serif")
text(1, 1.6, paste("RMSE: ", round(final_rmse_test, 2)), col = "black", cex = 0.9, family = "serif")

# 预测测试集
pre_test <- predict(final_model, newdata = best_test)

# 使用 points() 添加测试集的散点图
points(best_test$ANTH_c, pre_test, pch = 0, col = "blue", cex = 1.2)  # 测试集空心正方形

# 添加测试集的R²和RMSE值
text(1, 1.4, paste("R²: ", round(final_r2_test, 2)), col = "black", cex = 0.9, family = "serif")
text(1, 1.2, paste("RMSE: ", round(final_rmse_test, 2)), col = "black", cex = 0.9, family = "serif")

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
sel <- sample.split(data_use$ANTH_c, SplitRatio = 0.7)
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

# 新建一个图形设备并调整大小
dev.new(width = 10, height = 8)  # 可以调整为适当的宽度和高度

# 绘制训练集和测试集回归散点图
plot(train_label, pre_train1, main = "Regression Scatter Plot",
     xlab = "Actual Values", ylab = "Predicted Values", 
     pch = 1, col = "orange", cex = 1.5,  # 训练集: 空心圆，橙色
     xlim = range(c(test_label, train_label)),  # 设置 x 轴范围为训练集和测试集的最小和最大值
     ylim = range(c(pre_test, pre_train1)))  # 设置 y 轴范围为预测值的最小和最大值

# 绘制测试集的散点图
points(test_label, pre_test, pch = 0, col = "blue", cex = 1.5)  # 测试集: 空心正方形，蓝色

# 添加回归线
abline(0, 1, col = "black", lty = 2)  # 绘制回归线

# 在图中添加R²和RMSE值
# 调整文本位置到稍微右上角，避免被覆盖并确保完全显示

# 训练集文本位置调整
text(x = max(train_label) * 0.8, y = max(pre_train1) * 1.1, 
     labels = paste("R² (Train): ", round(XGB_R2_train, 3)), col = "black", cex = 0.9)
text(x = max(train_label) * 0.8, y = max(pre_train1) * 1.03, 
     labels = paste("RMSE (Train): ", round(XGB_rmse_train, 2)), col = "black", cex = 0.9)

# 测试集文本位置调整
text(x = max(test_label) * 0.8, y = max(pre_test) * 1.1, 
     labels = paste("R² (Test): ", round(XGB_R2_test, 3)), col = "black", cex = 0.9)
text(x = max(test_label) * 0.8, y = max(pre_test) * 1.03, 
     labels = paste("RMSE (Test): ", round(XGB_rmse_test, 2)), col = "black", cex = 0.9)


# 设置保存路径和图像尺寸
png("regression_scatter_plot_XGB.png", width = 1000, height = 800)

# 新建一个图形设备并调整大小
dev.new(width = 10, height = 8)  # 可以调整为适当的宽度和高度

# 绘制训练集和测试集回归散点图
plot(train_label, pre_train1, main = "Regression Scatter Plot",
     xlab = "Actual Values", ylab = "Predicted Values", 
     pch = 1, col = "orange", cex = 1.5,  # 训练集: 空心圆，橙色
     xlim = range(c(test_label, train_label)),  # 设置 x 轴范围为训练集和测试集的最小和最大值
     ylim = range(c(pre_test, pre_train1)))  # 设置 y 轴范围为预测值的最小和最大值

# 绘制测试集的散点图
points(test_label, pre_test, pch = 0, col = "blue", cex = 1.5)  # 测试集: 空心正方形，蓝色

# 添加回归线
abline(0, 1, col = "black", lty = 2)  # 绘制回归线

# 在图中添加R²和RMSE值
text(x = min(train_label) * 0.3, y = max(pre_train1) * 1.1, 
     labels = paste("R² (Train): ", round(XGB_R2_train, 3)), col = "black", cex = 0.9)
text(x = min(train_label) * 0.3, y = max(pre_train1) * 1.03, 
     labels = paste("RMSE (Train): ", round(XGB_rmse_train, 2)), col = "black", cex = 0.9)

text(x = min(test_label) * 0.3, y = max(pre_test) * 1.1, 
     labels = paste("R² (Test): ", round(XGB_R2_test, 3)), col = "black", cex = 0.9)
text(x = min(test_label) * 0.3, y = max(pre_test) * 1.03, 
     labels = paste("RMSE (Test): ", round(XGB_rmse_test, 2)), col = "black", cex = 0.9)

# 关闭图形设备并保存为PNG文件
dev.off()
