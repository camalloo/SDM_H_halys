test <- predsCV |> 
  dplyr::select(1:4, MXT_fold_3_pred, GLM_fold_3_pred, GAM_fold_3_pred) |> 
  mutate(MXT_fold_3_pred = log(MXT_fold_3_pred / (1-MXT_fold_3_pred)), 
         GLM_fold_3_pred = log(GLM_fold_3_pred / (1-GLM_fold_3_pred)),
         GAM_fold_3_pred = log(GAM_fold_3_pred / (1-GAM_fold_3_pred)))
## serie perfetta
low <- c(mean = 0, sd = 0.001)
high <- c(mean = 1, sd = 0.001)

n <- nrow(test)
test$quasiperfectFit <- case_when(
  test$presence == 0 ~ abs(0.001+rnorm(n, mean = low["mean"], sd = low["sd"])),
  test$presence == 1 ~ pmin(rnorm(n, mean = high["mean"], sd = high["sd"]), 0.99))

test <- test |> 
  mutate(perf = log(quasiperfectFit / (1 - quasiperfectFit)))


## maxent
form <- reformulate(termlabels = 'MXT_fold_3_pred', response = 'presence')
logistic <- glm(formula = form, family = binomial, data = test)
summary(logistic)
test$predMXT <- ifelse(predict(logistic, test, type = 'response') <= 0.5, 0, 1)
table(test$predMXT, test$presence)

## glm
form1 <- reformulate(termlabels = 'GLM_fold_3_pred', response = 'presence')
logistic1 <- glm(formula = form1, family = binomial, data = test)
summary(logistic1)
test$predGLM <- ifelse(predict(logistic1, test, type = 'response') <= 0.5, 0, 1)
table(test$predGLM, test$presence)

## synthetic
form3 <- reformulate(termlabels = 'perf', response = 'presence')
logistic3 <- glm(formula = form3, family = binomial, data = test)
summary(logistic3)
test$perfPred <- ifelse(predict(logistic3, test, type = 'response') <= 0.5, 0, 1)
table(test$perfPred, test$presence)

## gam
form2 <- reformulate(termlabels = 'GAM_fold_3_pred', response = 'presence')
logistic2 <- glm(formula = form2, family = binomial, data = test)
summary(logistic2)



## plot prediction
ggplot(data = test, aes(x = presence)) +
  geom_smooth(aes(y = perfPred),
              method = 'glm',
              method.args = list(family = binomial)) +
  geom_smooth(aes(y = predGLM),
              method = 'glm',
              method.args = list(family = binomial), col = 'black') +
  geom_smooth(aes(y = predMXT),
              method = 'glm',
              method.args = list(family = binomial), col = 'red')

## plot density
ggplot(data = test) +
  # geom_density(aes(x = MXT_fold_3_pred), col = 'red') +
  # geom_density(aes(x = GLM_fold_3_pred)) +
  # geom_density(aes(x = GAM_fold_3_pred), col = 'green') +
  geom_density(aes(x = perf), col = 'blue')

ggplot(data = test) +
  geom_smooth(aes(x = GLM_fold_3_pred, y = MXT_fold_3_pred), method = 'lm') +
  geom_abline()+
  coord_equal(ylim = c(0,1))

linear <- lm(GLM_fold_3_pred ~ MXT_fold_3_pred, data = test)



ggplot(data = test) +
  geom_density(aes(x = perfPred))
