# Creates 10x1001bp tiles per region (5 + strand, 5 - strand)
augTile <- function(bed, width= 1001, shifts= c(-400, -200, 0, 200, 400))
{
  # Hard copy ----
  bed <- data.table::copy(bed)
  
  # Augment regions using tiling ----
  aug <- lapply(shifts, function(x) {
    .c <- data.table::copy(bed)
    .c[, start:= (start+end)/2] # Center
    .c$start <- round(.c$start-width/2) # Start
    .c$start <- .c$start+x # shift
    .c$end <- .c$start+(width-1) # Extend
    .c[, shift:= x]
    # Plus minus strand
    .c[, strand:= "+"]
    ns <- data.table::copy(.c)
    ns[, strand:= "-"]
    .c <- rbind(.c, ns)
    return(.c)
  })
  
  # Clean ----
  aug <- rbindlist(aug)
  return(aug)
}

# Creates tiles for each region using a sliding window
augSlideFunction <- function(bed, width= 1001, step= 50)
{
  # Hard copy ----
  bed <- data.table::copy(bed)
  # Augment using sliding window ----
  bed$shift <- lapply(bed[, end-start+1], function(x) seq(0, x-(width-step), step))
  bed <- bed[, .(shift= unique(unlist(shift))), setdiff(names(bed), "shift")]
  bed[, start:= start+shift]
  bed[, end:= start+width-1]
  bed[, strand:= "+"]
  
  # Plus and minus strands
  ns <- data.table::copy(bed)
  ns[, strand:= "-"]
  bed <- rbind(bed, ns)
  
  # Return
  return(bed)
}
