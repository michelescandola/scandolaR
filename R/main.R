#' strInNum.
#'
#' Function that takes a factor that is numeric into a number.
#' The factors should only contain digits.
#'
#' @param input is the factor.
#' @returns the numeric version of the factor.
#' @export
#'
#' @examples
#' strInNum(as.factor(c("1", "2", "3"))) ## returns 1 2 3
strInNum <- function(input) as.numeric(as.character(input))

#' GimmeResult
#'
#' A function that takes the results of a t.test or other common tests and
#' returns a string formatted similarly to APA style
#'
#' @param test is the object resulting from a test.
#' @param digits a vector containing the number of required digits for the
#' degrees of freedom (second element in the vector), statistics and p-value
#' (first element of the vector). Default = c(3, 2)
#' @returns a string.
#' @export
#'
#' @examples
#' result <- t.test(1:10, 11:20)
#' GimmeResult(result)
GimmeResult <- function(test, digits = c(3, 2)) {
  dof <- ""
  if (!is.null(test$parameter)) dof <- paste("(", round(test$parameter, digits[2]), ")", sep = "")
  return(paste(test$method, " ", names(test$statistic), "$_{", dof, "}$", " = ",
    round(test$statistic, digits = digits[1]), "; p = ",
    round(test$p.value, digits = digits[1]), " ",
    sep = ""
  ))
}

#' GimmeResult2
#'
#' A function that takes the results of a chi squared test and
#' returns a string formatted similarly to APA style
#'
#' @param test is the object resulting from a test.
#' @param digits a vector containing the number of required digits for the
#' degrees of freedom (second element in the vector), statistics and p-value
#' (first element of the vector). Default = c(3, 2)
#' @returns a string.
#' @export
GimmeResult2 <- function(test, digits = c(3, 2)) {
  dof <- ""
  if (!is.null(test$parameter)) dof <- paste("(", round(test$parameter, digits[2]), ")", sep = "")
  return(paste(test$method, " OR= ", test$estimate, "; p = ",
    round(test$p.value, digits = digits[1]), " ",
    sep = ""
  ))
}

#' GimmeStar
#'
#' A function that reads the p-values and returns the asterisks
#'
#' @param p the p-value.
#' @param not_trend if FALSE returns "+" for p-values between 0.05 and 0.1
#' @returns a string.
#' @export
#'
#' @examples
#' GimmeStar(0.01) ## returns **
#' GimmeStar(0.010)
#' GimmeStar(0.07)
#' GimmeStar(0.07, not_trend = TRUE)
#' sapply(c(0.1, 0.04, 0.64, 0.02), GimmeStar) ## \"+\" \"*\" \"\"  \"*\"
GimmeStar <- function(p, not_trend = FALSE) {
  if (is.na(p)) {
    warning("presence of NAs")
    return("")
  }
  if (p > 0.10) {
    return("")
  }
  if (p > 0.05 & !not_trend) {
    return("+")
  }
  if (p > 0.01) {
    return("*")
  }
  if (p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

#' GimmeStarCor
#'
#' A function that reads correlations and returns the asterisks
#'
#' @param r the correlation value.
#' @param not_trend if FALSE returns "+" for rhos between 0.3 and 0.1
#' @returns a string.
#' @export
#'
#' @examples
#' GimmeStarCor(0.5) ## returns **
#' GimmeStarCor(0.1)
#' GimmeStarCor(0.2)
#' GimmeStarCor(0.2, not_trend = TRUE)
#' sapply(c(0.1, 0.04, 0.64, 0.02), GimmeStarCor)
GimmeStarCor <- function(r, not_trend = FALSE) {
  r2 <- abs(r)
  if (is.na(r)) {
    warning("presence of NAs")
    return("")
  }
  if (r2 > 0.50) {
    return("**")
  }
  if (r2 > 0.30) {
    return("*")
  }
  if (r2 > 0.10 & !not_trend) {
    return("+")
  }
  return("")
}

#' GimmeMSD
#'
#' A function that computes means and standard deviations (or standard error)
#' of a numeric vector and returns a string
#'
#' @param y the numerical vector.
#' @param digits the number of digits.
#' @param na.rm if they want to remove NAs from the computation
#' @param se if TRUE it will use the standard error instead of
#' the standard deviation
#' @returns a string.
#' @export
#'
#' @examples
#' GimmeMSD(c(1, 4, 8, 2, 3, 9, 10))
GimmeMSD <- function(y, digits = 3, na.rm = TRUE, se = FALSE) {
  if (!se) {
    return(paste0(
      round(mean(y, na.rm = na.rm), digits = digits),
      " (", round(sd(y, na.rm = na.rm), digits = digits), ")"
    ))
  } else {
    return(paste0(
      round(mean(y, na.rm = na.rm), digits = digits),
      " (", round(se(y, na.rm = na.rm), digits = digits), ")"
    ))
  }
}

#' Detect outliers within groups
#'
#' Identifies outlying observations in a dependent variable, either across
#' the entire dataset or separately within groups defined by one or more
#' independent variables.
#'
#' Several criteria can be used to identify outliers: the interquartile
#' range rule, percentile cut-offs, a specified number of standard deviations
#' from the mean, or a specified number of median absolute deviations from
#' the median.
#'
#' When requested, the function also creates a PDF file named
#' \code{outlier.pdf} containing diagnostic plots. If the number of groups is
#' larger than 20, groups are split across multiple plots.
#'
#' @param data A data frame containing the variables used in the analysis.
#'
#' @param IVs A character vector containing the names of the independent
#'   variables used to define groups. Outliers are detected separately within
#'   each combination of these variables. Use \code{"all"} to detect outliers
#'   across the entire distribution without grouping.
#'
#' @param DV A character string giving the name of the dependent variable in
#'   which outliers should be detected. Only one variable can be specified.
#'
#' @param criterion Character string specifying the criterion used to detect
#'   outliers. Available options are:
#'   \itemize{
#'     \item \code{"IQR"}: observations identified as outliers according to
#'       \code{\link[grDevices]{boxplot.stats}}.
#'     \item \code{"perc"}: observations below the 2.5th percentile or above
#'       the 97.5th percentile.
#'     \item \code{"norm-N"}: observations lying more than \code{N} standard
#'       deviations from the mean, for example \code{"norm-2"} or
#'       \code{"norm-3"}.
#'     \item \code{"MAD-N"}: observations lying more than \code{N} median
#'       absolute deviations from the median, for example \code{"MAD-2.5"}.
#'   }
#'   The default is \code{"IQR"}.
#'
#' @param graph Logical. If \code{TRUE}, diagnostic plots are written to a PDF
#'   file named \code{outlier.pdf} in the current working directory.
#'   The default is \code{TRUE}.
#'
#' @param out.csv Optional character string specifying the name or path of a
#'   CSV file in which the dataset used for the diagnostic plots should be
#'   saved. The exported dataset includes an \code{out} variable indicating
#'   whether each observation was classified as an outlier. If \code{NULL},
#'   no CSV file is written.
#'
#' @return An integer vector containing the row indices of the observations
#'   identified as outliers in \code{data}. If no outliers are detected,
#'   \code{NULL} is returned.
#'
#' @details
#' Groups are created using all combinations of the variables specified in
#' \code{IVs}. For example, if \code{IVs = c("Condition", "Sex")}, outliers
#' are detected independently within each Condition-by-Sex combination.
#'
#' If \code{IVs = "all"}, the grouping step is skipped and outliers are
#' detected using the complete distribution of \code{DV}.
#'
#' For the \code{"norm-N"} criterion, the lower and upper limits are:
#'
#' \deqn{mean(DV) \pm N \times SD(DV)}
#'
#' For the \code{"MAD-N"} criterion, the limits are:
#'
#' \deqn{median(DV) \pm N \times MAD(DV)}
#'
#' The function returns row numbers rather than the corresponding observations.
#' These indices can therefore be used directly to inspect or remove detected
#' observations.
#'
#' @note
#' The function requires the \pkg{ggplot2} package when
#' \code{graph = TRUE}.
#'
#' The criterion name for the MAD method is case-sensitive in the current
#' implementation and should contain \code{"MAD"} in uppercase.
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B"), each = 10),
#'   score = c(rnorm(9), 10, rnorm(9), -10)
#' )
#'
#' # Detect outliers separately within each group
#' rimozione.outlier(
#'   data = dat,
#'   IVs = "group",
#'   DV = "score",
#'   criterion = "IQR",
#'   graph = FALSE
#' )
#'
#' # Detect observations more than 2 SD from the group mean
#' rimozione.outlier(
#'   data = dat,
#'   IVs = "group",
#'   DV = "score",
#'   criterion = "norm-2",
#'   graph = FALSE
#' )
#'
#' # Detect outliers across the complete sample
#' rimozione.outlier(
#'   data = dat,
#'   IVs = "all",
#'   DV = "score",
#'   criterion = "perc",
#'   graph = FALSE
#' )
#'
#' @export
rimozione.outlier <- function(data,
                              IVs,
                              DV,
                              criterion = "IQR",
                              graph = TRUE,
                              out.csv = NULL) {
  # restituisce un vettore con gli indici degli outlier della
  # variabile dipendente
  # DV, divisi per le interazioni di IVs, con diversi criteri possibili:
  #   IQR   - Interquartile Rule
  #   perc  - rimuove i minore del 2.5o percentile e maggiori il 97.5 percentile
  #   norm-N- rimuove al di fuori di N*SD dalla media
  #   mad-N - rimuove al di fuori di N*MAD dalla mediana (Iglewicz & Hoaglin, 1993)
  # se IVs == "all" allora cerca gli outlier sull'intera distribuzione
  selezione.outliers <- NULL
  if (!(nrow(data) > 0 && ncol(data) > 0)) stop("Check your dataframe in data")
  if (!length(IVs) > 0 || !is.character(IVs)) stop("Please set in IVs a string or an array of strings")
  if (length(DV) > 1) stop("Please set in DV a single dependent varaible ")
  if (IVs[1] == "all") {
    tmp <- factor(rep("1", nrow(data)))
  } else {
    tmp <- interaction(data[, IVs])
  }

  for (ii in levels(tmp)) {
    sel <- which(tmp == ii)
    if (length(sel) > 0) {
      selout <- NULL
      if (criterion == "IQR") {
        outliers <- boxplot.stats(data[sel, DV])$out
        selout <- sel[data[sel, DV] %in% outliers]
      } else if (criterion == "perc") {
        lims <- quantile(data[sel, DV], probs = c(0.025, 0.975))
        selout <- sel[data[sel, DV] < lims[1] | data[sel, DV] > lims[2]]
      } else if (grepl("norm", criterion)) {
        N <- as.numeric(strsplit(criterion, "-")[[1]][2])
        lims <- c(
          mean(data[sel, DV], na.rm = TRUE) - N * sd(data[sel, DV], na.rm = TRUE),
          mean(data[sel, DV], na.rm = TRUE) + N * sd(data[sel, DV], na.rm = TRUE)
        )
        selout <- sel[data[sel, DV] < lims[1] | data[sel, DV] > lims[2]]
      } else if (grepl("MAD", criterion)) {
        N <- as.numeric(strsplit(criterion, "-")[[1]][2])
        lims <- c(
          median(data[sel, DV], na.rm = TRUE) - N * mad(data[sel, DV], na.rm = TRUE),
          median(data[sel, DV], na.rm = TRUE) + N * mad(data[sel, DV], na.rm = TRUE)
        )
        selout <- sel[data[sel, DV] < lims[1] | data[sel, DV] > lims[2]]
      } else {
        stop("criterion not recognized")
      }
      selezione.outliers <- c(selezione.outliers, selout)
    }
  }

  graphic <- data
  graphic$out <- "NO"
  graphic$out[selezione.outliers] <- "YES"
  graphic$out <- factor(graphic$out)
  graphic$IV <- tmp
  graphic$IV2 <- as.numeric(graphic$IV)
  graphic$y <- graphic[, DV]
  graphic$riga <- 1:nrow(graphic)
  mmax <- max(graphic$y)
  mmin <- min(graphic$y)
  if (graph) {
    nn <- (length(levels(tmp)) %/% 20) + 1
    graphic$div <- NA # factor(rep(1:nn,each=20,length.out=nrow(graphic)))
    for (inn in 1:nn) {
      sel <- which((graphic$IV2 < ((inn * 20) + 1)) & (graphic$IV2 >= (((inn - 1) * 20) + 1)))
      graphic$div[sel] <- inn
    }
    graphic$div <- factor(graphic$div)

    pdf("outlier.pdf")
    for (dd in levels(graphic$div)) {
      tttmp <- subset(graphic, div == dd)
      tttmp$IV <- droplevels(tttmp$IV)

      print(ggplot(tttmp, aes(y = y, x = IV, colour = out, label = riga)) +
        geom_text() +
        theme(axis.text.x = element_text(angle = 90)) +
        geom_boxplot(colour = "black", alpha = 0.5) +
        ggtitle(paste("Graphic number", dd)))
      # +
      #   coord_cartesian(ylim=c(mmin,mmax)))
      # scan(n=1)
      Sys.sleep(0.2)
    }
    dev.off()
  }
  if (graph) {
    warning(
      paste(
        "Printed",
        length(levels(graphic$div)),
        "graphics in", getwd(),
        " outliers.pdf"
      )
    )
  }
  if (!is.null(out.csv)) write.csv2(graphic, file = out.csv)
  return(selezione.outliers)
}

#' Visual inspection of a dependent variable by group
#'
#' Creates a separate boxplot for each group defined by one or more
#' independent variables. Each observation is labelled with its corresponding
#' row number in the original data frame, allowing potentially unusual or
#' outlying observations to be visually identified.
#'
#' @param data A data frame containing the variables used in the analysis.
#'
#' @param DV A character string specifying the name of the dependent variable
#'   to be plotted. Only one dependent variable can be specified.
#'
#' @param IVs A character vector containing the names of the independent
#'   variables used to define groups. When more than one variable is supplied,
#'   groups are defined by the interaction among all variables in \code{IVs}.
#'
#' @details
#' The function creates one plot for each combination of the levels of the
#' variables specified in \code{IVs}.
#'
#' Each plot contains a boxplot of \code{DV}, with every observation labelled
#' by its row number in the original \code{data} object. The same y-axis range
#' is used for all plots, making visual comparisons across groups easier.
#'
#' After each plot is displayed, the function pauses and waits for user input
#' before showing the next plot. Press Enter in the R console to continue.
#'
#' Groups are constructed internally using \code{\link[base]{interaction}}.
#'
#' @return
#' The function is used for its graphical output and does not explicitly
#' return an object.
#'
#' @note
#' The function requires the \pkg{ggplot2} package.
#'
#' Missing values in \code{DV} may cause problems when determining the common
#' y-axis limits because \code{min} and \code{max} are currently called
#' without \code{na.rm = TRUE}.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   condition = rep(c("A", "B"), each = 10),
#'   sex = rep(c("Female", "Male"), times = 10),
#'   score = rnorm(20)
#' )
#'
#' # Inspect score separately for each condition
#' visual.inspection(
#'   data = dat,
#'   DV = "score",
#'   IVs = "condition"
#' )
#'
#' # Inspect score for each Condition-by-Sex combination
#' visual.inspection(
#'   data = dat,
#'   DV = "score",
#'   IVs = c("condition", "sex")
#' )
#' }
#'
#' @export
visual.inspection <- function(data, DV, IVs) {
  # funzione che ti restituisce un grafico per ogni IVs dove ci sarà
  # il boxplot e gli indici di ogni data point

  if (!(nrow(data) > 0 && ncol(data) > 0)) stop("Check your dataframe in data")
  if (!length(IVs) > 0 || !is.character(IVs)) stop("Please set in IVs a string or an array of strings")
  if (length(DV) > 1) stop("Please set in DV a single dependent varaible ")
  tmp <- interaction(data[, IVs])
  tmp.data <- data
  tmp.data$riga <- 1:nrow(tmp.data)
  mmin <- min(data[, DV])
  mmax <- max(data[, DV])
  for (ii in levels(tmp)) {
    ttmp <- subset(tmp.data, tmp == ii) ## ogni singolo soggetto
    print(ggplot(ttmp, aes(y = ttmp[, DV], label = riga)) +
      geom_text(position = "dodge") +
      geom_boxplot(colour = "black", alpha = 0.5) +
      coord_cartesian(ylim = c(mmin, mmax)))
    scan(n = 1)
  }
}

## -------------------------------------------------------------------------
## UTILITY FUNCTIONS FOR REGRESSION
## -------------------------------------------------------------------------


#' Mean-center a numeric variable
#'
#' Centers a numeric vector by subtracting its mean while retaining the
#' original scale.
#'
#' @param x A numeric vector or matrix.
#'
#' @return A centered object with the same dimensions as \code{x}. For a
#'   vector, the returned object is a one-column matrix because the function
#'   relies on \code{\link[base]{scale}}.
#'
#' @details
#' This is a convenience wrapper around:
#'
#' \code{scale(x, scale = FALSE)}
#'
#' Each variable is centered around its mean but is not divided by its
#' standard deviation.
#'
#' @examples
#' x <- c(10, 12, 14, 16)
#' c.(x)
#'
#' mean(c.(x))
#'
#' @seealso
#' \code{\link[base]{scale}}, \code{\link{z.}}
#'
#' @export
c. <- function(x) scale(x, scale = FALSE)


#' Standardize a numeric variable
#'
#' Standardizes a numeric vector or matrix by subtracting the mean and
#' dividing by the standard deviation.
#'
#' @param x A numeric vector or matrix.
#'
#' @return A standardized object with the same dimensions as \code{x}. For a
#'   vector, the returned object is a one-column matrix because the function
#'   relies on \code{\link[base]{scale}}.
#'
#' @details
#' This is a convenience wrapper around \code{\link[base]{scale}} using its
#' default arguments. For a numeric variable \eqn{x}, the transformation is
#'
#' \deqn{z = (x - \bar{x}) / s}
#'
#' where \eqn{\bar{x}} is the sample mean and \eqn{s} is the sample standard
#' deviation.
#'
#' @examples
#' x <- c(10, 12, 14, 16)
#' z.(x)
#'
#' mean(z.(x))
#' sd(z.(x))
#'
#' @seealso
#' \code{\link[base]{scale}}, \code{\link{c.}}
#'
#' @export
z. <- function(x) scale(x)


#' Compute standardized residuals from a linear model
#'
#' Fits a linear model and returns its standardized residuals.
#'
#' @param formula A model formula passed to \code{\link[stats]{lm}}.
#'
#' @param ... Additional arguments passed to \code{\link[stats]{lm}}, such as
#'   \code{data}, \code{subset}, or \code{na.action}.
#'
#' @return A numeric vector containing the standardized residuals from the
#'   fitted linear model.
#'
#' @details
#' The function first fits a model using \code{\link[stats]{lm}} and then
#' applies \code{\link[stats]{rstandard}} to the fitted model.
#'
#' It is equivalent to:
#'
#' \code{rstandard(lm(formula, ...))}
#'
#' @examples
#' dat <- data.frame(
#'   y = c(2, 4, 5, 8, 10),
#'   x = c(1, 2, 3, 4, 5)
#' )
#'
#' r.(y ~ x, data = dat)
#'
#' @seealso
#' \code{\link[stats]{lm}}, \code{\link[stats]{rstandard}}
#'
#' @export
r. <- function(formula, ...) rstandard(lm(formula, ...))


#' Log-transform a numeric variable
#'
#' Applies the natural logarithm to a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector containing the natural logarithm of \code{x}.
#'
#' @details
#' This is a convenience wrapper around \code{\link[base]{log}}. Values must
#' satisfy the usual requirements of the logarithm: zero values return
#' \code{-Inf}, while negative values produce \code{NaN}.
#'
#' @examples
#' x <- c(1, 2, 5, 10)
#' l.(x)
#'
#' @seealso
#' \code{\link[base]{log}}
#'
#' @export
l. <- function(x) log(x)


#' Rescale a variable to the interval -1 to 1
#'
#' Linearly transforms a numeric variable so that its minimum becomes
#' \code{-1} and its maximum becomes \code{1}.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector rescaled to the interval \code{[-1, 1]}.
#'
#' @details
#' The transformation is:
#'
#' \deqn{
#' x^* =
#' \frac{2x - \max(x) - \min(x)}
#'      {\max(x) - \min(x)}
#' }
#'
#' This transformation has been suggested before computing polynomial terms
#' because rescaling the predictor can improve numerical properties.
#'
#' The function currently does not remove missing values when computing the
#' minimum and maximum. Therefore, the presence of \code{NA} values will
#' propagate through the transformation.
#'
#' If all values of \code{x} are identical, the denominator is zero and the
#' transformation is undefined.
#'
#' @references
#' Seber, G. A. F. (1977).
#' \emph{Linear Regression Analysis}.
#' Wiley.
#'
#' @examples
#' x <- seq(10, 20, by = 2)
#' s.(x)
#'
#' range(s.(x))
#'
#' @seealso
#' \code{\link{p.}}, \code{\link[base]{scale}}
#'
#' @export
s. <- function(x) {
  (2 * x - max(x) - min(x)) / (max(x) - min(x))
}


#' Generate orthogonal or raw polynomial terms
#'
#' Creates polynomial terms for a numeric predictor.
#'
#' @param x A numeric vector.
#'
#' @param ... Additional arguments passed to \code{\link[stats]{poly}}, such
#'   as \code{degree}, \code{raw}, or \code{simple}.
#'
#' @return The object returned by \code{\link[stats]{poly}}, usually a matrix
#'   containing polynomial terms.
#'
#' @details
#' This function is a convenience wrapper around \code{\link[stats]{poly}}.
#'
#' By default, \code{poly} returns orthogonal polynomial terms. To obtain raw
#' polynomial terms, use \code{raw = TRUE}.
#'
#' @examples
#' x <- 1:10
#'
#' # Orthogonal quadratic polynomial
#' p.(x, degree = 2)
#'
#' # Raw polynomial terms
#' p.(x, degree = 2, raw = TRUE)
#'
#' @seealso
#' \code{\link[stats]{poly}}, \code{\link{s.}}
#'
#' @export
p. <- function(x, ...) poly(x, ...)



## -------------------------------------------------------------------------
## UTILITY FUNCTIONS FOR LINEAR MIXED MODELS
## -------------------------------------------------------------------------


#' Variance inflation factors for a mixed-effects model
#'
#' Computes variance inflation factors for the fixed-effect coefficients of a
#' fitted mixed-effects model.
#'
#' @param fit A fitted mixed-effects model for which \code{vcov()} and
#'   \code{fixef()} methods are available.
#'
#' @return A named numeric vector containing the variance inflation factor
#'   associated with each fixed-effect coefficient, excluding the intercept.
#'
#' @details
#' The calculation is based on the covariance matrix of the fixed-effect
#' estimates. The intercept, when present, is excluded before the VIFs are
#' calculated.
#'
#' The covariance matrix is first converted to a correlation matrix. The VIF
#' for each coefficient is then obtained from the diagonal of the inverse
#' correlation matrix.
#'
#' This implementation was adapted from the approach used in
#' \code{rms::vif}.
#'
#' Large VIF values indicate that a coefficient is strongly linearly related
#' to other predictors in the model. Interpretation thresholds should not be
#' treated as absolute rules and should be considered together with the model
#' design.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#'
#' fit <- lmer(
#'   Reaction ~ Days + I(Days^2) + (1 | Subject),
#'   data = sleepstudy
#' )
#'
#' vif.mer(fit)
#' }
#'
#' @seealso
#' \code{\link{colldiag.mer}}, \code{\link{kappa.mer}},
#' \code{\link{maxcorr.mer}}
#'
#' @export
vif.mer <- function(fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))

  # exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))

  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }

  d <- diag(v)^0.5
  v <- diag(solve(v / (d %o% d)))

  names(v) <- nam

  v
}


#' Condition number for a mixed-effects model
#'
#' Computes the condition number of the fixed-effects design matrix of a
#' mixed-effects model.
#'
#' @param fit A fitted mixed-effects model containing a fixed-effects design
#'   matrix.
#'
#' @param scale Logical. If \code{TRUE}, columns of the design matrix are
#'   divided by their standard deviations before the condition number is
#'   calculated. The default is \code{TRUE}.
#'
#' @param center Logical. If \code{TRUE}, columns of the design matrix are
#'   centered before the calculation. The default is \code{FALSE}.
#'
#' @param add.intercept Logical. If \code{TRUE}, an intercept column is added
#'   after the original model intercept has been removed. The default is
#'   \code{TRUE}.
#'
#' @param exact Logical passed to \code{\link[base]{kappa}}. If \code{TRUE},
#'   the exact 2-norm condition number is computed. If \code{FALSE}, the
#'   default approximation used by \code{kappa} is used.
#'
#' @return A single numeric value containing the condition number of the
#'   fixed-effects design matrix.
#'
#' @details
#' The condition number is a diagnostic for collinearity in a design matrix.
#' Larger values indicate that the columns of the matrix are increasingly
#' linearly dependent.
#'
#' The model intercept is first removed from the original fixed-effects
#' design matrix. Depending on \code{add.intercept}, an intercept may then be
#' added back after scaling or centering.
#'
#' This function accesses the fixed-effects model matrix through
#' \code{fit@X}, which corresponds to an older internal representation of
#' mixed-effects model objects. For newer versions of \pkg{lme4},
#' \code{getME(fit, "X")} is generally preferable.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#'
#' fit <- lmer(
#'   Reaction ~ Days + I(Days^2) + (1 | Subject),
#'   data = sleepstudy
#' )
#'
#' kappa.mer(fit)
#' }
#'
#' @seealso
#' \code{\link[base]{kappa}}, \code{\link{vif.mer}},
#' \code{\link{colldiag.mer}}
#'
#' @export
kappa.mer <- function(fit,
                      scale = TRUE,
                      center = FALSE,
                      add.intercept = TRUE,
                      exact = FALSE) {
  X <- fit@X
  nam <- names(fixef(fit))

  # exclude intercepts
  nrp <- sum(1 * (nam == "(Intercept)"))

  if (nrp > 0) {
    X <- X[, -(1:nrp), drop = FALSE]
    nam <- nam[-(1:nrp)]
  }

  if (add.intercept) {
    X <- cbind(
      rep(1),
      scale(X, scale = scale, center = center)
    )

    kappa(X, exact = exact)
  } else {
    kappa(
      scale(X, scale = scale, center = scale),
      exact = exact
    )
  }
}


#' Collinearity diagnostics for mixed-effects models
#'
#' Computes condition indices and variance-decomposition proportions for the
#' fixed-effects design matrix of a mixed-effects model or for a supplied
#' numeric matrix.
#'
#' @param fit A fitted mixed-effects model, matrix, or data frame containing
#'   the predictors to be evaluated.
#'
#' @param scale Logical. If \code{TRUE}, columns of the design matrix are
#'   scaled before the singular value decomposition. The default is
#'   \code{TRUE}.
#'
#' @param center Logical. If \code{TRUE}, columns are centered before the
#'   analysis. When \code{TRUE}, \code{add.intercept} is automatically set to
#'   \code{FALSE}. The default is \code{FALSE}.
#'
#' @param add.intercept Logical. If \code{TRUE}, an intercept column is added
#'   to the design matrix when one is not already included. The default is
#'   \code{TRUE}.
#'
#' @return A data frame containing:
#'
#' \describe{
#'   \item{\code{cond.index}}{
#'     The condition index associated with each singular dimension.
#'   }
#'   \item{Remaining columns}{
#'     Variance-decomposition proportions for the corresponding model
#'     coefficients.
#'   }
#' }
#'
#' @details
#' The function implements a collinearity diagnostic based on the approach of
#' Belsley, Kuh, and Welsch (1980), using the singular value decomposition of
#' the design matrix.
#'
#' Condition indices are calculated as the ratio between the largest singular
#' value and each singular value:
#'
#' \deqn{
#' CI_j = d_{max} / d_j
#' }
#'
#' Collinearity is of particular concern when a large condition index is
#' accompanied by large variance-decomposition proportions for two or more
#' coefficients.
#'
#' A commonly used heuristic is to inspect condition indices above
#' approximately 30 together with the corresponding variance proportions.
#' This threshold should be interpreted as a diagnostic guideline rather than
#' a formal statistical test.
#'
#' The implementation was adapted from the approach used by
#' \code{perturb::colldiag}.
#'
#' @references
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980).
#' \emph{Regression Diagnostics: Identifying Influential Data and Sources of
#' Collinearity}. Wiley.
#'
#' @examples
#' x <- data.frame(
#'   x1 = 1:20,
#'   x2 = (1:20) + rnorm(20, sd = 0.5),
#'   x3 = rnorm(20)
#' )
#'
#' colldiag.mer(x)
#'
#' \dontrun{
#' library(lme4)
#'
#' fit <- lmer(
#'   Reaction ~ Days + I(Days^2) + (1 | Subject),
#'   data = sleepstudy
#' )
#'
#' colldiag.mer(fit)
#' }
#'
#' @seealso
#' \code{\link{vif.mer}}, \code{\link{kappa.mer}}
#'
#' @export
colldiag.mer <- function(fit,
                         scale = TRUE,
                         center = FALSE,
                         add.intercept = TRUE) {
  result <- NULL

  if (center) {
    add.intercept <- FALSE
  }

  if (is.matrix(fit) || is.data.frame(fit)) {
    X <- as.matrix(fit)
    nms <- colnames(fit)
  } else if (inherits(fit) == "mer") {
    nms <- names(fixef(fit))
    X <- fit@X

    if (any(grepl("(Intercept)", nms))) {
      add.intercept <- FALSE
    }
  }

  X <- X[!is.na(apply(X, 1, all)), ]

  if (add.intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "(Intercept)"
  }

  X <- scale(
    X,
    scale = scale,
    center = center
  )

  svdX <- svd(X)

  condindx <- max(svdX$d) / svdX$d
  dim(condindx) <- c(length(condindx), 1)

  Phi <- svdX$v %*% diag(1 / svdX$d)
  Phi <- t(Phi^2)

  pi <- prop.table(Phi, 2)

  colnames(condindx) <- "cond.index"

  if (!is.null(nms)) {
    rownames(condindx) <- nms
    colnames(pi) <- nms
    rownames(pi) <- nms
  } else {
    rownames(condindx) <- 1:length(condindx)
    colnames(pi) <- 1:ncol(pi)
    rownames(pi) <- 1:nrow(pi)
  }

  result <- data.frame(
    cbind(condindx, pi)
  )

  zapsmall(result)
}


#' Maximum absolute correlation among fixed-effect estimates
#'
#' Extracts the correlation matrix of the fixed-effect coefficient estimates
#' from a mixed-effects model and returns the correlation with the largest
#' absolute magnitude.
#'
#' @param fit A fitted mixed-effects model.
#'
#' @param exclude.intercept Logical. If \code{TRUE}, correlations involving
#'   the intercept are excluded. The default is \code{TRUE}.
#'
#' @return A single numeric value corresponding to the fixed-effect
#'   correlation with the largest absolute magnitude. The original sign of
#'   the correlation is retained.
#'
#' @details
#' Strong correlations among estimated fixed-effect coefficients may indicate
#' substantial dependence among predictors or model terms.
#'
#' The function considers only one triangle of the fixed-effect correlation
#' matrix to avoid duplicated correlations.
#'
#' This implementation accesses the correlation matrix using slots associated
#' with older \pkg{lme4} model objects:
#'
#' \code{summary(fit)@vcov@factors$correlation}
#'
#' This representation is not compatible with many current \pkg{lme4}
#' versions. In modern code, the correlation matrix can instead be derived
#' from \code{vcov(fit)}.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#'
#' fit <- lmer(
#'   Reaction ~ Days + I(Days^2) + (1 | Subject),
#'   data = sleepstudy
#' )
#'
#' maxcorr.mer(fit)
#' }
#'
#' @seealso
#' \code{\link{vif.mer}}, \code{\link{colldiag.mer}}
#'
#' @export
maxcorr.mer <- function(fit,
                        exclude.intercept = TRUE) {
  so <- summary(fit)
  corF <- so@vcov@factors$correlation
  nam <- names(fixef(fit))

  # exclude intercepts
  ns <- sum(
    1 * (nam == "Intercept" | nam == "(Intercept)")
  )

  if (ns > 0 & exclude.intercept) {
    corF <- corF[
      -(1:ns),
      -(1:ns),
      drop = FALSE
    ]

    nam <- nam[-(1:ns)]
  }

  corF[!lower.tri(corF)] <- 0

  maxCor <- max(corF)
  minCor <- min(corF)

  if (abs(maxCor) > abs(minCor)) {
    zapsmall(maxCor)
  } else {
    zapsmall(minCor)
  }
}


#' Add informative names to random-effects PCA results
#'
#' Performs a principal component analysis of the random-effects covariance
#' structure of a mixed-effects model and assigns informative column names to
#' the resulting PCA summaries.
#'
#' @param fit A fitted mixed-effects model supported by
#'   \code{\link[lme4]{rePCA}} and \code{\link[lme4]{VarCorr}}.
#'
#' @return A list containing the summarized random-effects PCA results.
#'   Column names of the \code{importance} matrices are modified to identify
#'   the corresponding grouping factor and random-effect term.
#'
#' @details
#' The function applies \code{\link[lme4]{rePCA}} to examine the dimensionality
#' of the estimated random-effects covariance structure.
#'
#' Names are constructed by combining the grouping-factor name and the
#' corresponding random-effect coefficient name, separated by \code{"*"}.
#'
#' For example, a random intercept and random slope for \code{Days} grouped by
#' \code{Subject} may produce names such as:
#'
#' \code{"Subject*(Intercept)"} and \code{"Subject*Days"}.
#'
#' Principal component analysis of the random-effects covariance matrix can be
#' useful when investigating overparameterized or singular mixed-effects
#' models.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#'
#' fit <- lmer(
#'   Reaction ~ Days + (Days | Subject),
#'   data = sleepstudy
#' )
#'
#' rePCA.names(fit)
#' }
#'
#' @seealso
#' \code{\link[lme4]{rePCA}},
#' \code{\link[lme4]{VarCorr}},
#' \code{\link[lme4]{isSingular}}
#'
#' @export
rePCA.names <- function(fit) {
  obj <- summary(rePCA(fit))
  model <- VarCorr(fit)

  if (length(obj) == length(model)) {
    obj <- Map(
      function(x, z) {
        colnames(x$importance) <- paste(
          z,
          unique(sapply(model, colnames)),
          sep = "*"
        )

        x
      },
      obj,
      names(obj)
    )
  } else if (length(obj) == 1) {
    colnames(obj[[1]]$importance) <- unlist(
      mapply(
        paste,
        names(model),
        sapply(model, colnames),
        MoreArgs = list(sep = "*")
      )
    )
  }

  return(obj)
}


#' Convert an F statistic to partial eta squared
#'
#' Converts an observed F statistic and its degrees of freedom to partial
#' eta squared.
#'
#' @param f A numeric value or vector containing F statistics.
#'
#' @param df A numeric value or vector containing the numerator degrees of
#'   freedom associated with each F statistic.
#'
#' @param df_error A numeric value or vector containing the denominator,
#'   or error, degrees of freedom associated with each F statistic.
#'
#' @return A numeric value or vector containing partial eta squared
#'   (\eqn{\eta_p^2}).
#'
#' @details
#' Partial eta squared is calculated from an F statistic as:
#'
#' \deqn{
#' \eta_p^2 =
#' \frac{F \times df}
#'      {F \times df + df_{error}}
#' }
#'
#' The function is vectorized, so multiple F statistics can be converted in
#' a single call.
#'
#' @examples
#' # F = 5.4, numerator df = 2, denominator df = 97
#' my_F_to_eta2(
#'   f = 5.4,
#'   df = 2,
#'   df_error = 97
#' )
#'
#' # Convert several F statistics
#' my_F_to_eta2(
#'   f = c(2.1, 5.4, 10.2),
#'   df = c(1, 2, 1),
#'   df_error = c(98, 97, 98)
#' )
#'
#' @export
my_F_to_eta2 <- function(f, df, df_error) {
  (f * df) / (f * df + df_error)
}

#' Extract the legend from a ggplot object
#'
#' Extracts the legend grob from a plot created with \pkg{ggplot2}.
#' This can be useful when arranging multiple plots while displaying
#' a shared legend only once.
#'
#' @param a.gplot A \code{ggplot} object containing a legend.
#'
#' @return A graphical object (grob) containing the legend of
#'   \code{a.gplot}.
#'
#' @details
#' The function builds the supplied ggplot using
#' \code{\link[ggplot2]{ggplot_build}}, converts it to a graphical table
#' using \code{\link[ggplot2]{ggplot_gtable}}, and extracts the grob named
#' \code{"guide-box"}, which contains the plot legend.
#'
#' The supplied plot must contain a legend. If no legend is present,
#' the function may return an indexing error because no
#' \code{"guide-box"} grob can be found.
#'
#' This function is particularly useful when combining several ggplot
#' objects and a common legend needs to be extracted and positioned
#' separately.
#'
#' @examples
#' library(ggplot2)
#'
#' p <- ggplot(
#'   iris,
#'   aes(
#'     x = Sepal.Length,
#'     y = Sepal.Width,
#'     colour = Species
#'   )
#' ) +
#'   geom_point()
#'
#' legend <- g_legend(p)
#' legend
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot_build}},
#' \code{\link[ggplot2]{ggplot_gtable}}
#'
#' @importFrom ggplot2 ggplot_build ggplot_gtable
#' @export
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(
    sapply(tmp$grobs, function(x) x$name) == "guide-box"
  )
  legend <- tmp$grobs[[leg]]

  return(legend)
}
