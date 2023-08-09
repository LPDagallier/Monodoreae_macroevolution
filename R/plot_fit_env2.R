plot_fit_env2 <- function (fit.env, env_data, tot_time) 
{
  if (!inherits(fit.env, "fit.env")) 
    stop("object is not of class \"fit.env\"")
  t <- seq(0, tot_time, length.out = 100)
  df <- smooth.spline(env_data[, 1], env_data[, 2])$df
  spline_result <- pspline::sm.spline(env_data[, 1], env_data[, 2], 
                             df = df)
  env_func <- function(t) {
    predict(spline_result, t)
  }

  if ("f.mu" %in% attributes(fit.env)$names) {
    par(mfrow = c(3,2))
    plot(-t, fit.env$f.lamb(t), type = "l", xlab = "time", ylab = "speciation rate", 
         main = "Fitted speciation rate")
    plot(env_func(t), fit.env$f.lamb(t), xlab = "Environmental data", 
         ylab = "speciation rate", main = "Fitted speciation rate")
    plot(-t, fit.env$f.mu(t), type = "l", xlab = "time", 
         ylab = "extinction rate", main = "Fitted extinction rate")
    plot(env_func(t), fit.env$f.mu(t), xlab = "Environmental data", 
         ylab = "extinction rate", main = "Fitted extinction rate")
    r <- function(t) {
      fit.env$f.lamb(t) - fit.env$f.mu(t)
    }
    plot(-t, r(t), type = "l", xlab = "time", ylab = "net diversification rate", 
         main = "Fitted net diversification rate")
    plot(env_func(t), r(t), xlab = "Environmental data", 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
  }
  else {
    par(mfrow = c(1,2))
    # plot(-t, fit.env$f.lamb(t), type = "l", xlab = "time", ylab = "speciation rate", 
    #      main = "Fitted speciation rate")
    # plot(env_func(t), fit.env$f.lamb(t), xlab = "Environmental data", 
    #      ylab = "speciation rate", main = "Fitted speciation rate")
    plot(-t, fit.env$f.lamb(t), type = "l", xlab = "time", 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
    plot(env_func(t), fit.env$f.lamb(t), xlab = "Environmental data", 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
  }
}
