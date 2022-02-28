################################################################################
# Useful functions for working with motus data
# 09-Jan-2020
# adapted from Ariel Lenske
#
################################################################################

# set custom theme for all plots

theme_ajm <- function() {
  theme_classic() %+replace%
    theme(axis.title.x = element_text(size=10),
          axis.text.x  = element_text(size=8, colour = "black"),
          axis.title.y = element_text(size=10, angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.text.y = element_text(size=8, colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          strip.text.x = element_text(size=12, face = "bold"),
          legend.text = element_text(size=8),
          legend.key.height = unit(1, "mm"),
          plot.title = element_text(size = 12, hjust = 0, vjust = 1.5),
          panel.border = element_rect(size =0.5, fill = "transparent"),
          plot.margin = margin(10, 10, 10, 15))
}

################################################################################

# function for mapping tower detections - simplifies so not too much data to plot####
fun.getpath <- function(df) 
{
  df %>%
    filter(!is.na(recvLat) | !(recvLat == 0)) %>% # drops data without lon/lat
    group_by(motusTagID, runID, recvName, ambigID, 
             tagDepLon, tagDepLat, recvLat, recvLon, tagDeployStart, motusFilter) %>%
    #summarizing by runID to get max run length and mean time stamp:
    summarize(max.runLen = max(runLen), ts.h = mean(ts)) %>% 
    mutate(year = as.numeric(format(ts.h,'%Y')),
           month = as.numeric(format(ts.h,'%m'))) %>%
    arrange(motusTagID, ts.h) %>%
    ungroup()
} # end of function call

# function for plotting signal strength against time by recv
plotTagSig_mod <- function (df, motusTagID) 
{
  tag.id <- motusTagID
  df <- filter(df, motusTagID == tag.id)
  df <- select(df, motusTagID, sig, noise, ts, recvLat, 
               recvLon, recvName, runLen, motusFilter)
  df <- mutate(df, 
               recvName = paste(recvName, 
                                round(recvLat, digits = 1), sep = "\n"), 
               recvName = paste(recvName, 
                                round(recvLon, digits = 1), sep = ", "), 
               ts = as_datetime(ts, tz = "UTC")
               # antBearing.ch = if_else(is.na(antBearing),"none", as.character(antBearing))
  )
  df <- within(df, recvName <- reorder(recvName, (recvLat)))
  
  ggplot(df, aes(ts, sig, col = as.factor(motusFilter))) +
    geom_point() +
    geom_point(aes(ts, noise), col = "lightgrey") +
    theme_bw() + 
    labs(title = paste("Detection Time vs Signal Strength, coloured by motusFilter \n ID ", 
                       motusTagID),
         x = "Date", 
         y = "Signal Strength", 
         colour = "motusFilter") + 
    facet_grid(recvName ~ .) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} # end of function call

# calculate distance between two points
# used to be in motus package?
latLonDist = function(lat1, lon1, lat2, lon2) {
  a = 6378137
  b = 6356752.314245
  f = 1/298.257223563  ## WGS-84 ellipsoid params
  
  llmat = cbind(lat1, lon1, lat2, lon2) ## recycles coordinates to match
  
  s = rep(-1, nrow(llmat)) ## return values; -1 means not yet computed
  for (i in 1:nrow(llmat)) {  ## calculate distance between i'th pair of points
    if (!all(is.finite(llmat[i,]))) {
      s[i] = NA
      next
    }
    
    L = rad(llmat[i, 4]-llmat[i, 2])
    U1 = atan((1-f) * tan(rad(llmat[i, 1])))
    U2 = atan((1-f) * tan(rad(llmat[i, 3])))
    sinU1 = sin(U1)
    cosU1 = cos(U1)
    sinU2 = sin(U2)
    cosU2 = cos(U2)
    lambda = L
    iterLimit = 100
    repeat {
      sinLambda = sin(lambda)
      cosLambda = cos(lambda)
      sinSigma = sqrt((cosU2*sinLambda) * (cosU2*sinLambda) + 
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda))
      if (abs(sinSigma) < 1e-10) {
        s[i] = 0 ## co-incident points
        break
      }
      cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
      sigma = atan2(sinSigma, cosSigma)
      sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
      cosSqAlpha = 1 - sinAlpha*sinAlpha
      cos2SigmaM = cosSigma - 2*sinU1*sinU2 / cosSqAlpha
      if (is.nan(cos2SigmaM))
        cos2SigmaM = 0  ## equatorial line: cosSqAlpha=0 (§6)
      C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
      lambdaP = lambda
      lambda = L + (1-C) * f * sinAlpha *
        (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
      iterLimit = iterLimit - 1
      if (abs(lambda-lambdaP) <= 1e-12 || iterLimit == 0)
        break
    } 
    
    if (iterLimit==0) {
      s[i] = NaN  ## formula failed to converge
    } else if (s[i] < 0) {
      uSq = cosSqAlpha * (a*a - b*b) / (b*b)
      A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
      B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
      deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
                                                 B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))
      s[i] = b*A*(sigma-deltaSigma)
    }
  }
  s = round(s, 3)
  return (s)
}

# radians
rad = function(x) {
  return (x * (pi/180))
}

################################################################################
#'functions used to calculate site transition statistics

################################################################################
#'main function:
################################################################################

#'modifided siteTrans function (from the Motus package) - Summarize transitions between sites for each tag

site_transit_min <- function(data, latCoord = "recvLat", lonCoord = "recvLon"){
  tmp <- data
  data <- rename(tmp, lat = latCoord, lon = lonCoord)
  
  data <- data %>% select(ts, motusTagID, lat, lon, recvName, recvProjID, month, year) %>% #pull out relevent columns
    dplyr::group_by(motusTagID) %>% do(consec.fun(.)) %>% #re-format data step#1
    dplyr::group_by(motusTagID) %>% do(site.fun(.)) %>% #re-format data step#2
    mutate(tot_ts = round(difftime(ts.y, ts.x, units = "secs"))) %>% #time diff between detections in seconds
    mutate(tot_ts = ifelse(tot_ts == 0, 1, tot_ts), #deal with simultaneous detections by making them one sec apart
           dist = round(latLonDist(lat.x, lon.x, lat.y, lon.y)), #meters
           dist.min = dist - 50000) %>% #account for detection ranges of towers
    mutate(dist.min = ifelse(dist.min < 0, 0, dist.min),
           rate = round(dist.min/(as.numeric(tot_ts)), digits = 4)) %>% #rate of travel in m/s
    mutate(suspect.transit = ifelse(rate > 72, "suspect", "fine"))
  
  return(data)
}

site_transit_stopover <- function(data, latCoord = "recvLat", lonCoord = "recvLon"){
  tmp <- data
  data <- rename(tmp, lat = latCoord, lon = lonCoord)
  
  data <- data %>% 
    do(add_deploy_loc_rgroup(.)) %>%
    select(ts, motusTagID, lat, lon, recvName, recvProjID, recvGroup, month, year) %>% #pull out relevent columns
    dplyr::group_by(motusTagID) %>% do(consec.fun(.)) %>% #re-format data step#1
    dplyr::group_by(motusTagID) %>% do(site.fun(.)) %>% #re-format data step#2
    mutate(tot_ts = round(difftime(ts.y, ts.x, units = "secs"))) %>% #time diff between detections in seconds
    mutate(tot_ts = ifelse(tot_ts == 0, 1, tot_ts), #deal with simultaneous detections by making them one sec apart
           dist = round(latLonDist(lat.x, lon.x, lat.y, lon.y))) %>% #meters
    mutate(rate = round(dist/(as.numeric(tot_ts)), digits = 4)) %>% #rate of travel in m/s
    mutate(stopover = ifelse(rate < 5, "yes", "no"))
  
  return(data)
}




###############################################################################
#'functions that are used in main function:
###############################################################################

#'function to add row for tag deployment location/time to use for calculating transitions

add_deploy_loc <- function(df){
  
  df.r1 <- df %>% 
    group_by(motusTagID) %>%
    arrange(ts) %>%
    slice(1) %>%
    mutate(ts = tagDeployStart,
           recvName = "deploy location",
           recvLat = tagDepLat, 
           recvLon = tagDepLon,
           month = as.numeric(format(ts,'%m')),
           year = as.numeric(format(ts,'%Y'))) 
  
  df <- bind_rows(df.r1, df) %>% group_by(motusTagID) %>% arrange(ts) %>% as.data.frame()
  
}

add_deploy_loc_rgroup <- function(df){
  
  df.r1 <- df %>% 
    group_by(motusTagID) %>%
    arrange(ts) %>%
    slice(1) %>%
    mutate(ts = tagDeployStart,
           recvName = "deploy location",
           recvLat = tagDepLat, 
           recvLon = tagDepLon,
           recvGroup = "deploy location",
           month = as.numeric(format(ts,'%m')),
           year = as.numeric(format(ts,'%Y'))) 
  
  df <- bind_rows(df.r1, df) %>% group_by(motusTagID) %>% arrange(ts) %>% as.data.frame()
  
}



################################################################################
#' Create dataframe for siteTrans function
#'
#' @param data dataframe of Motus detection data containing at a minimum fullID, ts, lat, lon
#'
#' @author Zoe Crysler \email{zcrysler@@gmail.com}
#'

## site.fun and consec.fun adapted from "between.locs.R" script written by Phil

consec.fun <- function(df) {
  
  df <- df[order(df$ts),]
  a <- df$recvName[-length(df$recvName)]
  b <- df$recvName[-1]
  tmp <- c(0, 1 - (as.numeric(a==b)))
  run <- cumsum(tmp)
  transitions <- which(diff(run) != 0)
  transitions <- c(transitions, transitions+1, length(df$recvName))
  out.df <- df[transitions,]
  out.df <- out.df[order(out.df$ts),]
  return(out.df)
  
}

##################################################################################
#' Create dataframe for siteTrans function
#'
#' @param data dataframe of Motus detection data containing at a minimum fullID, ts, lat, lon
#'
#' @author Zoe Crysler \email{zcrysler@@gmail.com}
#'

## site.fun and consec.fun adapted from "between.locs.R" script written by Phil
site.fun <- function(df) {
  
  df <- subset(df, select = -c(motusTagID))
  df <- df[order(df$ts),] ## should already be in order, but just in case
  out.df.x <- df[1:(length(df$recvName)-1), ]
  names(out.df.x) <- paste(names(df), "x", sep=".")
  out.df.y <- df[2:length(df$recvName), ]
  names(out.df.y) <- paste(names(df), "y", sep=".")
  out.df <- cbind(out.df.x, out.df.y)
  out.df <- subset(out.df, ((recvName.x != recvName.y)))
  
  
  return(out.df)
}

####################################################################

# function to check deployment coordinates against coordinates for each receiver and select
# the minimum distance to a receiver and assign the region of that receiver to the deployment
getDepRegion <- function(birds){
  
  dep.regions <- list()
  
  for (i in 1:length(sesa.deps.bird)){
    
    bird.dep <- sesa.tag.deps %>% 
      filter(motusTagID == sesa.deps.bird[i])
    
    bird.dep.reg <- recv.locs %>% 
      mutate(motusTagID = bird.dep$motusTagID,
             tagDepLat = bird.dep$tagDepLat,
             tagDepLon = bird.dep$tagDepLon) %>% 
      mutate(dist = latLonDist(recvLat, recvLon, tagDepLat, tagDepLon)) %>% 
      filter(dist == min(dist))
    
    dep.regions[[i]] <- bird.dep.reg
    
  }
  
  return(dep.regions)
  
}

################################################################################
# modified correlation plot

# function body
ggcorrplot_mod <- function(corr,
                           method = c("square", "circle"),
                           type = c("full", "lower", "upper"),
                           ggtheme = ggplot2::theme_minimal,
                           title = "",
                           show.legend = TRUE,
                           legend.title = "Corr",
                           show.diag = FALSE,
                           colors = c("blue", "white", "red"),
                           outline.color = "gray",
                           hc.order = FALSE,
                           hc.method = "complete",
                           lab = FALSE,
                           lab_col = "black",
                           lab_size = 4,
                           p.mat = NULL,
                           sig.level = 0.05,
                           insig = c("pch", "blank"),
                           pch = 4,
                           pch.col = "black",
                           pch.cex = 5,
                           tl.cex = 12,
                           tl.col = "black",
                           tl.srt = 45,
                           digits = 2,
                           as.is = FALSE) {
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  
  if(inherits(corr, "cor_mat")){
    # cor_mat object from rstatix
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  
  corr <- base::round(x = corr, digits = digits)
  
  if (hc.order) {
    ord <- .hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  
  # Get lower or upper triangle
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  
  # Melt corr and pmat
  corr <- reshape2::melt(corr, na.rm = TRUE, as.is = as.is)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value > sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  
  
  corr$abs_corr <- abs(corr$value) * 10
  
  # heatmap
  p <-
    ggplot2::ggplot(
      data = corr,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")
    )
  
  # modification based on method
  if (method == "square") {
    p <- p +
      ggplot2::geom_tile(color = outline.color)
  } else if (method == "circle") {
    p <- p +
      ggplot2::geom_point(
        color = outline.color,
        shape = 21,
        ggplot2::aes_string(size = "abs_corr")
      ) +
      ggplot2::scale_size(range = c(4, 10)) +
      ggplot2::guides(size = FALSE)
  }
  
  # adding colors
  p <-
    p + ggplot2::scale_fill_gradient2(
      low = colors[1],
      high = colors[3],
      mid = colors[2],
      midpoint = 0.2,
      limit = c(0, 0.4),
      space = "Lab",
      name = legend.title
    )
  
  # depending on the class of the object, add the specified theme
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  } else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  
  
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = tl.srt,
        vjust = 1,
        size = tl.cex,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = tl.cex)
    ) +
    ggplot2::coord_fixed()
  
  label <- round(x = corr[, "value"], digits = digits)
  if(!is.null(p.mat) & insig == "blank"){
    ns <- corr$pvalue > sig.level
    if(sum(ns) > 0) label[ns] <- " "
  }
  
  # matrix cell labels
  if (lab) {
    p <- p +
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
        label = label,
        color = lab_col,
        size = lab_size
      )
  }
  
  # matrix cell glyphs
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(
      data = p.mat,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
      shape = pch,
      size = pch.cex,
      color = pch.col
    )
  }
  
  # add titles
  if (title != "") {
    p <- p +
      ggplot2::ggtitle(title)
  }
  
  # removing legend
  if (!show.legend) {
    p <- p +
      ggplot2::theme(legend.position = "none")
  }
  
  # removing panel
  p <- p +
    .no_panel()
  p
}



#' Compute the matrix of correlation p-values
#'
#' @param x numeric matrix or data frame
#' @param ... other arguments to be passed to the function cor.test.
#' @rdname ggcorrplot
#' @export

cor_pmat <- function(x, ...) {
  
  # initializing values
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # creating the p-value matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # name rows and columns of the p-value matrix
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  # return the final matrix
  p.mat
}



#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}


# Convert a tbl to matrix
.tibble_to_matrix <- function(x){
  x <-  as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}