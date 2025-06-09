plot_predictive_positive <- function(label,
                                     predicted,
                                     xlim= NULL,
                                     ylim= NULL,
                                     plot= FALSE,
                                     ...) {
  # Create a data table with observed and predicted values
  dat <- data.table(obs = label,
                    pred = predicted)
  
  # Sort the data table by the predicted values in descending order
  setorderv(dat,
            cols = "pred",
            order = -1)
  
  # Calculate the cumulative percentage of actually positive values
  dat[, perc:= sapply(seq(obs), function(i) sum(obs[1:i]) / i * 100)]
  setorderv(dat,
            cols = "pred")
  
  # Fit a smooth spline and find peaks
  x <- dat[1:(.N-100), pred]
  y <- dat[1:(.N-100), perc]
  spline_fit <- smooth.spline(x, y)
  spline_derivative <- predict(spline_fit, deriv = 1)
  
  # Find peaks
  peaks <- which(diff(sign(spline_derivative$y)) == -2) + 1
  peak_x_values <- spline_derivative$x[peaks]
  peak_y_values <- predict(spline_fit, x = peak_x_values)$y
  peak_x_values <- c(peak_x_values, last(x)) # value at .N-100
  peak_y_values <- c(peak_y_values, last(y)) # value at .N-100 
  
  # Select ideal cutoff
  sel <- min(which(peak_y_values>max(0.95*peak_y_values)))
  x_cutoff <- peak_x_values[sel]
  y_cutoff <- peak_y_values[sel]
  
  # Plot PPV and cutoffs
  if(plot)
  {
    dat[, {
      # initiate plot
      if(is.null(xlim))
        xlim <- range(pred, na.rm= T)
      if(is.null(ylim))
        ylim <- range(perc, na.rm= T)
      plot(NA,
           type = "n",
           xlab = "Prediction score",
           ylab = "Positive pred. value (%)",
           xlim= xlim,
           ylim= ylim,
           ...)
      
      # Plot lines
      lines(pred[1:(.N-100)], perc[1:(.N-100)])
      lines(pred[(.N-100):.N], perc[(.N-100):.N], lty = "33")
    }]
    
    # Plot the smooth spline
    lines(spline_fit, col = "blue", lwd = 2)
    points(x_cutoff,
           y_cutoff,
           col = "red",
           pch = 19)
    segments(x_cutoff,
             0,
             x_cutoff,
             y_cutoff, lty= "33")
    text(x_cutoff,
         y_cutoff/2,
         round(x_cutoff, 2),
         pos= 4,
         col= "red",
         offset= 1)
    segments(0,
             y_cutoff,
             x_cutoff,
             y_cutoff, lty= "33")
    text(x_cutoff/2,
         y_cutoff,
         paste0(round(y_cutoff, 1), "%"),
         pos= 3,
         col= "red",
         offset= 1)
  }

  # Return cutoffs
  return(list(min_PPV= dat$perc[1],
              predict_cutoff= x_cutoff,
              PPV_at_cutoff= y_cutoff))
}