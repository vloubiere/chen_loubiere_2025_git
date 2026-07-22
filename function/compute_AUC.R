# Function to compute AUC ----
vl_ROC_AUC <- function(label, predicted, plot.line= F, col= "red", ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  
  # Make data table ----
  dat <- data.table(label= label,
                    predicted= predicted)
  
  # Order ----
  setorderv(dat, "predicted", -1)
  
  # Compute AUC ----
  dat[, pos:= sum(label)]
  dat[, neg:= .N-pos]
  dat[, TPR:= cumsum(label)/pos]
  dat[, FPR:= cumsum(!label)/neg]
  dat[, ROC_AUC:= c(0, sapply(seq(.N)[-1], function(i) (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1]) / 2))]

  # Plot ----
  if(plot.line)
    lines(dat$FPR,
          dat$TPR,
          col= col, 
          ...)
  
  # AUC function
  # AUC <- pROC::roc(AUC$label, AUC$score)$auc
  
  return(sum(dat$ROC_AUC))
}
