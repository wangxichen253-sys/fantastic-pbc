#PBC_Predictive_Model_Code#
#user:DongZhan,XichenWang
#Last_modified_date:2026-6-29

####1.Load_the_necessary_R_packages####
library(survminer) 
library(survival) 
library(corrplot) 
library(glmnet) 
library(caret) 
library(CBCgrps)
library(nortest)
library(tidyverse)
library(ggpubr)
library(rms)
library(pROC)
library(mice)
library(zoo)
library(VIM)
library(ggplot2)
library(stringr)
library(dplyr)

####2_Data_organization_and_cleaning####
###2.1####
#data <- read.csv()
data <- na.omit(data)
###2.2####
missing_percentage <- colMeans(is.na(data)) * 100
selected_columns <- names(missing_percentage[missing_percentage <=20])
cleaned_data <- data[selected_columns]
###2.3####
data <- mice(cleaned_data,m=5,maxit=50,meth='rf')
summary(data)
data <- complete(data,1)
aggr(data)
write.csv(data,"data.csv")
###2.4####
two_value_cols <- sapply(data, function(x) length(unique(x))) == c(2,3,4,5,6,7,8,9,10,11,12)
data[, two_value_cols] <- lapply(data[, two_value_cols], as.factor)
data <- data[order(substr(rownames(data), 1, 1)), ]

####3_Descriptive_Statistical_Analysis####
###3.1_Draw_a_baseline_data_table_PFS####
data <- as.data.frame(data)#有时候需要
tab1 <-twogrps(data, gvar = "pfs", sim = TRUE)
write.csv(tab1$Table,file="tab1.csv")
###3.2_Draw_a_heat_map(Univariate_Correlation_Analysis, numeric)####
numdata <- data[,-216]
numdata[,1:217] <- lapply(data[,1:217],as.numeric)
numdata = data%>%
  select(where(is.numeric))
numdata <- numdata[,-c(1,2)]
M<-cor(numdata)
M
testRes<-cor.mtest(numdata,conf.level = 0.95)
testRes
par(mfrow=c(2,3))
corrplot(M,method ='circle')
dev.off()#again and again
corrplot(M,method ='circle')
corrplot(
  M,
  method='color',
  type = 'upper',
  add = T ,
  tl.pos = "n",
  cl.pos = "n",
  diag = F,
  p.mat = testRes$p,
  sig.level = c(0.001,0.01,0.05),
  pch.cex = 1.5,
  insig = 'label_sig'
)
###3.3_Draw_a_heat_map(Data_Presentation)####
data <- data[order(substr(rownames(data), 1, 1)), ]
##3.3.1_Function-normalize_clinical_data####
library(dplyr)
library(tidyr)
library(purrr)
normalize_clinical_data <- function(data,
                                    max_cat_levels = 5,
                                    log_vars = NULL,
                                    binary_ref = "min0") {
  data <- as.data.frame(data)
  type_map <- list()
  norm_list <- list()
  for (colname in names(data)) {
    x <- data[[colname]]
    # 字符/逻辑 -> 因子
    if (is.character(x) || is.logical(x)) {
      x <- as.factor(x)
    }
    # 缺失率过高警告
    if (mean(is.na(x)) > 0.9) {
      warning(sprintf("列 '%s' 缺失值超过90%%", colname))
    }
    # --- 因子处理 ---
    if (is.factor(x)) {
      n_lvl <- length(levels(x))
      if (is.ordered(x)) {
        type_map[[colname]] <- "ordered"
        num_x <- as.numeric(x) - 1
        norm_list[[colname]] <- num_x / max(num_x, na.rm = TRUE)
      } else if (n_lvl == 2) {
        type_map[[colname]] <- "binary"
        if (binary_ref == "min0") {
          norm_list[[colname]] <- ifelse(as.numeric(x) == 1, 0, 1)
        } else {
          norm_list[[colname]] <- ifelse(as.numeric(x) == 2, 1, 0)
        }
      } else {
        # 多分类 -> one-hot
        type_map[[colname]] <- "categorical"
        dummies <- model.matrix(~ . - 1, data = data.frame(x))
        colnames(dummies) <- paste0(colname, "_", levels(x))
        for (dcol in colnames(dummies)) {
          norm_list[[dcol]] <- as.numeric(dummies[, dcol])
          type_map[[dcol]] <- "binary"
        }
      }
      # --- 数值处理 ---
    } else if (is.numeric(x)) {
      uniq <- setdiff(unique(x), NA)
      n_uniq <- length(uniq)
      # 仅含0/1 -> 二分类
      if (n_uniq <= 2 && all(uniq %in% c(0, 1))) {
        type_map[[colname]] <- "binary"
        norm_list[[colname]] <- ifelse(x == 0, 0, 1)
        # 少量整数 -> 有序
      } else if (n_uniq <= max_cat_levels &&
                 all(uniq == floor(uniq), na.rm = TRUE)) {
        type_map[[colname]] <- "ordered"
        rng <- range(x, na.rm = TRUE)
        norm_list[[colname]] <- (x - rng[1]) / (rng[2] - rng[1])
        # 其余 -> 连续
      } else {
        type_map[[colname]] <- "continuous"
        if (colname %in% log_vars) {
          x <- log1p(x)
        }
        rng <- range(x, na.rm = TRUE)
        if (rng[2] == rng[1]) {
          norm_list[[colname]] <- rep(0.5, length(x))
        } else {
          norm_list[[colname]] <- (x - rng[1]) / (rng[2] - rng[1])
        }
      }
    } else {
      warning(sprintf("列 '%s' 类型无法识别，跳过", colname))
      next
    }
  }
  norm_df <- as.data.frame(norm_list, stringsAsFactors = FALSE)
  rownames(norm_df) <- rownames(data)
  type_df <- stack(type_map) %>%
    select(variable = ind, type = values) %>%
    as.data.frame()
  list(norm_df = norm_df, type_map = type_df)
}
##3.3.2_normalization####
raw_data<- data
drop_cols <- c("os")#no feature
feature_df <- raw_data %>% select(-any_of(drop_cols))
library(e1071)
num_cols <- names(feature_df)[sapply(feature_df, is.numeric)]
skew_vals <- sapply(feature_df[num_cols], function(x) {
  x_clean <- x[!is.na(x)]
  if(length(unique(x_clean)) < 5) return(NA)  # 唯一值太少，非连续
  skewness(x_clean)
})
high_skew_vars <- names(which(skew_vals > 1))
log_vars <- high_skew_vars

result <- normalize_clinical_data(
  data = feature_df,
  max_cat_levels = 5,
  log_vars = log_vars,
  binary_ref = "min0"
)
data_standardized <- result$norm_df
data1 <- data_standardized
##3.3.3_Remove_dummy_variables####
original_cat_vars <- result$type_map %>% filter(type == "categorical") %>% pull(variable)
cat_dummy_cols <- names(result$norm_df)[
  grepl(paste0("^(", paste(original_cat_vars, collapse = "|"), ")_"), names(result$norm_df))
]
data_heat <- result$norm_df %>% select(-all_of(cat_dummy_cols))
data <- data_heat
##3.3.4_Draw_a_heat_map####
library(pheatmap)
data1 <- as.data.frame(data1)
rownames(data1) <- rownames(data)
data1 <- data1[,-1]
heatmap1 <- pheatmap(data1,
                     color = colorRampPalette(c("blue", "white", "red"))(100), 
                     cluster_rows = FALSE,     
                     cluster_cols = FALSE,     
                     show_rownames = TRUE,    
                     show_colnames = TRUE,     
                     fontsize = 8,             
                     border_color = NA,angle_col = 90)
dev.off()
col_names <- colnames(data1)
col_names[seq(2, length(col_names), by = 2)] <- paste0(col_names[seq(2, length(col_names), by = 2)], "—————")
colnames(data1) <- col_names
pheatmap(data1,
         color = colorRampPalette(c("#94c6c9","#94c6c9", "#e2eed2", "#c82935"))(100), 
         cluster_rows = F,     
         cluster_cols = F,     
         show_rownames = F, 
         show_colnames = T,     
         fontsize = 6,     
         border_color = "grey",angle_col = 90)
####4_Lasso_Regression####
###4.1_Normalization####
min_max_scale = function(x){
  (x-min(x))/(max(x)-min(x))
}
data2 = data%>%
  mutate_if(.predicate = is.numeric,
            .funs = min_max_scale)%>%
  as.data.frame()
###4.2_Variable_Selection####
set.seed(123) #random number generator
x <- data.matrix(data2[, -1])
x <- x[, -1]
y <- data.matrix(Surv(data$time,data$pfs))
###4.3####
lasso <- glmnet(x, y, family = "cox",nlambda = 1000, alpha = 1)
print(lasso)
plot(lasso, xvar = "lambda", label = TRUE)
###4.4_Cross_Validation####
lasso.cv = cv.glmnet(x, y,alpha = 1,nfolds =20,family="cox")
plot(lasso.cv)
lasso.cv$lambda.min #minimum
lasso.cv$lambda.1se #one standard error away
###4.5_output_coefficients####
coef_min <- coef(lasso.cv, s = "lambda.min")
coef_df <- data.frame(
  variable = rownames(as.matrix(coef_min)),
  coefficient = as.matrix(coef_min)[, 1],
  row.names = NULL)
coef_df_nonzero <- coef_df[coef_df$coefficient != 0, ]
#coef_df_nonzero <- coef_df_nonzero[-15, ]
write.csv(coef_df, file = "lasso_coefficients_all.csv", row.names = FALSE)
write.csv(coef_df_nonzero, file = "lasso_coefficients_nonzero.csv", row.names = FALSE)
###4.6_Plot_LASSO_regression_coefficient_graph####
#Extract_Variable
Coefficients <- coef(lasso.cv, s = lasso.cv$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
lassonames <- row.names(Coefficients)[Active.Index]
reformulate(lassonames)
# plot
library(MetBrewer)
df_lasso <- data.frame(Variable = lassonames, Coefficient = Active.Coefficients) %>% 
  filter(Variable != "Intercept") %>% arrange(Coefficient)
df_lasso$Variable <- factor(df_lasso$Variable, levels = df_lasso$Variable)
clabs <- sprintf("%.3f", df_lasso$Coefficient)
p_lpop <- ggplot(df_lasso, aes(Coefficient, Variable)) +
  geom_segment(aes(yend = Variable, xend = 0), linewidth = 0.8) +
  geom_point(aes(size = abs(Coefficient), color = Coefficient)) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_color_gradientn(colors = met.brewer("VanGogh2")) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  scale_y_discrete(sec.axis = dup_axis(labels = clabs)) +
  labs(x = "LASSO Coefficient", y = "NULL") +
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 7)) +
  theme(plot.margin = unit(rep(0.5, 4), "cm"), panel.background = element_blank(),
        panel.grid = element_line(color = "grey40", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks.y.right = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x.right = element_text(size = 10, hjust = 1, color = ifelse(df_lasso$Coefficient > 0, "red", "#004c6d")),
        axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 13, face = "bold"))
p_lpop
###4.7.1_Lasso_Correlation_Analysis1####
numdata <- data[, c("pfs", as.character(coef_df_nonzero$variable)), drop = FALSE]
numdata <- numdata[,-18]#delete character
numdata[,1:17] <- lapply(numdata[,1:17],as.numeric)
M<-cor(numdata)
M
testRes<-cor.mtest(numdata,conf.level = 0.95)
testRes
par(mfrow=c(2,3))
corrplot(M,method ='circle')
dev.off()
corrplot(M,method ='circle')
corrplot(
  M,
  method='color',
  type = 'upper',
  add = T ,
  tl.pos = "n",
  cl.pos = "n",
  diag = F,
  p.mat = testRes$p,
  sig.level = c(0.001,0.01,0.05),
  pch.cex = 1.5,
  insig = 'label_sig'
)
###4.7.2_Lasso_Correlation_Analysis2####

####5_Split_the_dataset####
set.seed(1)
train_id = sample(1:nrow(data),0.85*nrow(data))#Set Parameters
train=data[train_id,]
test=data[-train_id,]
save(data, train, test, file = "model_data.rda")

####6_Cox####
###6.1_Integrate_Data & set_Outcome_Variable####
#coef_df_nonzero <- read.csv(..)#input lasso
coef_df_nonzero <- coef_df_nonzero[-17,]
attach(train)
dd<-datadist(train) 
options(datadist='dd') #packed
train$status <- as.numeric(train$pfs)#outcome_variable
###6.2_multi_factor1####
#nr <- "y~a+b+c+d"
#reformulate(colnames(train),"y")
#coxm <- cph(Surv(train$time,train$pfs)~ ultrasound.cirrhosis +
#              ultrasound.liver.injury + pathologic.diagnosis,
#            x=T,y=T,data=train,surv=T) 
#coxm
###6.3_multi_factor2####
nr <- "y~a+b+c+d"
reformulate(colnames(train),"y")
vars <- coef_df_nonzero$variable
formula <- as.formula(
  paste("Surv(time, status) ~", paste(vars, collapse = " + "))
)
coxm <- cph(formula,
            x = TRUE, y = TRUE,
            data = train,
            surv = TRUE) 
coxm
###6.4_save####
model_stats <- data.frame(
  LR_chi2 = coxm$stats["Model L.R."],
  df = coxm$stats["d.f."],
  P_value = coxm$stats["P"],
  R2 = coxm$stats["R2"]
)
###6.5.1_Nomogram_of_Median_Survival_Time####
med <- Quantile(coxm)
surv <- Survival(coxm)
nom <- nomogram(coxm,fun=function(x) med(lp=x),
                funlabel = "Median Survival Time")
plot(nom)
###6.5.2_dynamic_nomogram####
library(regplot)
coxm <- cph(formula,
            x = TRUE, y = TRUE,
            data = train,
            surv = TRUE)
regplot(coxm,
        observation = train[4, ],# appointed patient
        points = T,                     
        plots = c("density", "no plot"), 
        failtime = c(12, 12 * 3, 12 * 5, 12 * 8),# time
        leftlabel = T,                
        prfail = TRUE,# if,Cox,TRUE; else
        showP = T,                     
        droplines = T,                 
        #colors = "orange",# color
        rank = "range",                 
        interval = "confidence",
        title = "Nomogram"
        # plots = c("violin", "boxes")  
        # clickable = T                
)
###6.5.3_Counting_socre####
options(datadist="dd")
library(rms)  # needed for nomogram
library(nomogramFormula)
results <- formula_lp(nom)  # 基于线性预测值计算
results$formula
points <- points_cal(formula = results$formula, lp = coxm$linear.predictors)
head(points)
results <- formula_rd(nomogram = nom)  # 基于自变量值计算
str(data)
results$formula
points = points_cal(formula = results$formula, rd = data)
head(data)
data$points <- points
###6.6_KM####
library(ggplot2)
library(survminer)
library(survival)
data$points2 <- ifelse(data$points > median(data$points), 'High-Score', 'Low-Score')
data$pfs <- as.numeric(data$pfs) 
fit <- survfit(Surv(time, pfs) ~ points2, data = data)
p <- ggsurvplot(fit, 
                data = data,               
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                legend.labs = c('High-Score', 'Low-Score'),
                legend.title = 'Risk Score',
                palette = c('dodgerblue2', 'orchid2'),
                risk.table.height = 0.3)
p
###6.7_Survival_Probability_Nomogram####
surv1 <- function(x)surv(1*25,lp=x) 
surv2 <- function(x)surv(1*50,lp=x)
surv3 <- function(x)surv(1*75,lp=x)
surv4 <- function(x)surv(1*100,lp=x)
plot(nomogram(coxm,
              fun=list(surv1,surv2,surv3,surv4),
              lp= F,
              funlabel=c('25 months Survival','50 months Survival','75 months Survival','100 months Survival'),
              maxscale=100,
              fun.at=c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','0.95')))
###6.8_C_index####
rcorrcens(Surv(time,pfs)~predict(coxm),data=train)
###6.9_Calibration_Curve(again)####
train$status <- as.numeric(train$pfs)
vars <- coef_df_nonzero$variable
formula <- as.formula(
  paste("Surv(time, pfs) ~", paste(vars, collapse = " + "))
)
coxm1 <- cph(formula,
            x = TRUE, y = TRUE,
            data = train,
            surv = TRUE,time.inc = 75) #time
cal<- calibrate(coxm1, cmethod = 'KM', method = 'boot', 
                u = 25, m = 50, B = 100)#m_sample,u_time
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 75 months",#
     ylab="Actual 75-months live proportion", #
     col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,
      col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
###6.10_timeROC####
library(timeROC)
library(survivalROC)
train$status <- as.factor(train$pfs)
time_roc <- timeROC(
  T=train$time,
  delta=train$pfs,
  marker=coxm$linear.predictors,
  cause = 1,
  weighting = "marginal",
  times = c(25,50,75),#time
  ROC=TRUE,
  iid=TRUE)
time_ROC_df <- data.frame(
  TP_25_months = time_roc$TP[, 1],#time1
  FP_25_months = time_roc$FP_1[, 1],#time1
  TP_50_months = time_roc$TP[, 2],#time2
  FP_50_months = time_roc$FP_2[, 2],#time2
  TP_75_months = time_roc$TP[, 3],#time3
  FP_75_months = time_roc$FP_2[, 3])#time3
#25months
{auc_25  <- time_roc[["AUC_1"]][["t=25"]]
sd_25   <- time_roc[["inference"]][["vect_sd_1"]][["t=25"]]
ci95low <- sprintf("%.3f", auc_25 - 1.96 * sd_25)
ci95up  <- sprintf("%.3f", auc_25 + 1.96 * sd_25)
ci9525months <- paste0("(", ci95low, "-", ci95up, ")")
#50months
auc_50  <- time_roc[["AUC_1"]][["t=50"]]
sd_50   <- time_roc[["inference"]][["vect_sd_1"]][["t=50"]]
ci95low <- sprintf("%.3f", auc_50 - 1.96 * sd_50)
ci95up  <- sprintf("%.3f", auc_50 + 1.96 * sd_50)
ci9550months <- paste0("(", ci95low, "-", ci95up, ")")
#75months
auc_75  <- time_roc[["AUC_1"]][["t=75"]]
sd_75   <- time_roc[["inference"]][["vect_sd_1"]][["t=75"]]
ci95low <- sprintf("%.3f", auc_75 - 1.96 * sd_75)
ci95up  <- sprintf("%.3f", auc_75 + 1.96 * sd_75)
ci9575months <- paste0("(", ci95low, "-", ci95up, ")")}
ggplot(data = time_ROC_df) +
  geom_line(aes(y = FP_25_months, x = TP_25_months), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(y = FP_50_months, x = TP_50_months), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(y = FP_75_months, x = TP_75_months), linewidth = 1, color = "#6bc72b") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 25 months = ", sprintf("%.3f", 1-time_roc$AUC_1[[1]]),ci9525months), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 50 months = ", sprintf("%.3f", 1-time_roc$AUC_1[[2]]),ci9550months), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 75 months = ", sprintf("%.3f", 1-time_roc$AUC_1[[3]]),ci9575months), color = "#6bc72b"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(t = 0, r = 15, b = 0, l = 0))
  )

####7_COX_test####
###7.2_timeROC####
library(timeROC)
library(survivalROC)
train$status <- as.factor(train$pfs)
time_roc_test <- timeROC(
  T = test$time,
  delta = test$pfs,
  marker = predict(coxm, newdata = test, type = "lp"),
  cause = 1,
  weighting = "marginal",
  times = c(25,50,75),#time
  ROC = TRUE,
  iid = TRUE
)
time_ROC_test_df <- data.frame(
  TP_25_months = time_roc_test$TP[, 1],#time1
  FP_25_months = time_roc_test$FP_1[, 1],#time1
  TP_50_months = time_roc_test$TP[, 2],#time2
  FP_50_months = time_roc_test$FP_2[, 2],#time2
  TP_75_months = time_roc_test$TP[, 3],#time3
  FP_75_months = time_roc_test$FP_2[, 3])#time3
#25months
{auc_25  <- time_roc_test[["AUC_1"]][["t=25"]]
sd_25   <- time_roc_test[["inference"]][["vect_sd_1"]][["t=25"]]
ci95low <- sprintf("%.3f", auc_25 - 1.96 * sd_25)
ci95up  <- sprintf("%.3f", auc_25 + 1.96 * sd_25)
ci9525months <- paste0("(", ci95low, "-", ci95up, ")")
#50months
auc_50  <- time_roc_test[["AUC_1"]][["t=50"]]
sd_50   <- time_roc_test[["inference"]][["vect_sd_1"]][["t=50"]]
ci95low <- sprintf("%.3f", auc_50 - 1.96 * sd_50)
ci95up  <- sprintf("%.3f", auc_50 + 1.96 * sd_50)
ci9550months <- paste0("(", ci95low, "-", ci95up, ")")
#75months
auc_75  <- time_roc_test[["AUC_1"]][["t=75"]]
sd_75   <- time_roc_test[["inference"]][["vect_sd_1"]][["t=75"]]
ci95low <- sprintf("%.3f", auc_75 - 1.96 * sd_75)
ci95up  <- sprintf("%.3f", auc_75 + 1.96 * sd_75)
ci9575months <- paste0("(", ci95low, "-", ci95up, ")")}
ggplot(data = time_ROC_test_df) +
  geom_line(aes(y = FP_25_months, x = TP_25_months), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(y = FP_50_months, x = TP_50_months), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(y = FP_75_months, x = TP_75_months), linewidth = 1, color = "#6bc72b") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 25 months = ", sprintf("%.3f", 1-time_roc_test$AUC_1[[1]]),ci9525months), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 50 months = ", sprintf("%.3f", 1-time_roc_test$AUC_1[[2]]),ci9550months), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 75 months = ", sprintf("%.3f", 1-time_roc_test$AUC_1[[3]]),ci9575months), color = "#6bc72b"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(t = 0, r = 15, b = 0, l = 0))
  )
###7.1_Calibration_Curve(again)####
attach(test)
dd<-datadist(test) 
options(datadist='dd') 
test$status <- as.numeric(test$pfs)
nr <- "y~a+b+c+d"
reformulate(colnames(test),"y")
vars <- coef_df_nonzero$variable
formula <- as.formula(
  paste("Surv(time, status) ~", paste(vars, collapse = " + "))
)
coxm <- cph(formula,
            x = TRUE, y = TRUE,
            data = test,
            surv = TRUE) 
coxm
#again
test$status <- as.numeric(test$pfs)
vars <- coef_df_nonzero$variable
formula <- as.formula(
  paste("Surv(time, status) ~", paste(vars, collapse = " + "))
)
coxm1 <- cph(formula,
             x = TRUE, y = TRUE,
             data = test,
             surv = TRUE,time.inc = 50) #time
cal<- calibrate(coxm1, cmethod = 'KM', method = 'boot', 
                u = 50, m = 30, B = 100)#m_sample,u_time
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 50 months",#
     ylab="Actual 50-months live proportion", #
     col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,
      col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
###7.3_baseline_data_table_Cox_train_test####
train$grp <- 1
test$grp <- 0
datatt <- rbind(train,test)
datatt <- as.data.frame(datatt)
library(CBCgrps)
tab1 <-twogrps(datatt, gvar = "grp", sim = TRUE)
write.csv(tab1$Table,file="tab1.csv")