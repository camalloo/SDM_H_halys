logo <- rast(system.file("ex/logo.tif", package="terra"))
plot(logo[[1]])
p <- prcomp(logo)
summary(p)
x <- predict(logo, p)
plot(x)

test <- crop(preds[[8]], world[world$GID_0 == 'ARG'])
plot(test)
plot(world, add =T)

v='gddlgd10'
t <- rast("data/CHELSAbioclim/CHELSA_gddlgd10_1981-2010_V.2.1.tif")
plot(t)
plot(world, add = T)
c <- crop(t, world[world$GID_0 == 'NAM'])
plot(c, colNA = 'red')
plot(world, add = T)
d <- data.frame(c)
e <- crop(c, ext(22, 24, -28, -27))
plot(e, colNA = 'red')


NAs <- app(bio[[9:11]], fun = function(x) any(is.na(x)))
plot(NAs)

ff <- subst(bio[[21]], NA, 0)
plot(ff)



y = 2000
bartRCs <- embarcadero::partial(bartMod, trace = FALSE, smooth = 3)
bartRCs_df <- do.call(rbind, lapply(bartRCs, function(i) {
  df  <- i$data
  var <- i@labels$title   # taken from the object, not from vars
  
  data.frame(
    year     = y,
    variable = var,
    x        = df$x,
    fit      = df$med,
    ci_low   = df$q05,
    ci_high  = df$q95
  )
}))

# Inspect the structure before assuming anything
class(bartRCs)
class(bartRCs[[1]])
str(bartRCs[[1]])         # full structure of first element
names(bartRCs[[1]]$data)  # column names of the underlying data
head(bartRCs[[1]]$data)   # see actual values and column names


pdf(NULL)
gamRCs <- plot(gamMod, pages = 1, shade = FALSE)
dev.off()



rfRCs <- list()
for (var in vars) {
  thisRC <- do.call(randomForest::partialPlot,
                    list(x = rfMod, pred.data = data, x.var = var, plot = F))
  rfRCs[[var]] <- thisRC
}

library('pdp')
library('future.apply')


bartRCs_df <- future_lapply(vars, function(v) {
  embarcadero::partial(bartMod,
                       var = v,
                       trace = FALSE,
                       smooth = 5)
})

# set up parallel plan
plan(multisession, workers = parallel::detectCores() - 1)

rfRCs_list <- future_lapply(vars, function(var) {
  library(randomForest)
  thisRC <- pdp::partial(
    object = rfMod,
    pred.var = var,
    train = data,
    grid.resolution = 100,
    prob = TRUE
  )
  
  data.frame(
    year     = y,
    variable = var,
    x        = thisRC[[var]],
    fit      = thisRC$yhat
  )
}) |> bind_rows()

rfRCs_df <- rfRCs_list |> 
  bind_rows()

pkg <- as.data.frame(installed.packages())

library(dplyr)



gamRCs <- plot(gamMod, pages = 1, shade = F)
gamRCs_df <- do.call(rbind, lapply(gamRCs, function(i) {
  data.frame(
    year     = y,
    variable = i$xlab,
    x        = i$x,
    fit      = i$fit,
    se       = i$se,
    ci_low   = i$fit - 1.96 * i$se,
    ci_high  = i$fit + 1.96 * i$se
  )
}))
summary(gamMod)
gam.check(gamMod)


bartRCs_df <- do.call(rbind, future_lapply(vars, function(var) {
  library(embarcadero)
  options(mc.cores = 1)  # prevent BART internal threading conflicting with future
  plot_obj <- embarcadero::partial(bartMod,
                                   x.vars = var,
                                   equal  = TRUE,
                                   smooth = 5,
                                   trace  = FALSE,
                                   ci     = TRUE)
  df <- plot_obj$data
  data.frame(
    year     = y,
    variable = plot_obj@labels$title,
    x        = df$x,
    fit      = df$med,
    ci_low   = df$q05,
    ci_high  = df$q95
  )
}))

