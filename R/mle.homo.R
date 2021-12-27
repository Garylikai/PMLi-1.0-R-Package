#' @title Ekbohm's MLE-Based Test under Homoscedasticity
#' @description \code{mle.homo} is used to perform Ekbohm's MLE-based test
#'   under homoscedasticity for partially matched samples, specified by giving
#'   a data frame or matrix, testing hypothesis, true mean and confidence level.
#' @param data a 2-by-n or n-by-2 partially matched pairs samples data frame or
#'   matrix.
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" or "less". You can specify
#'   just the initial letter.
#' @param mu a number indicating the true difference in means.
#' @param conf.level confidence level for the returned confidence interval.
#' @param ...	 further arguments to be passed to or from methods.
#' @details Whether the dataset input is partially matched samples will be
#'   checked. If not, appropriate hypothesis tests, such as one and two sample
#'   t-tests, will be performed instead.
#' @return A list with class "htest" containing the following components:
#'   \item{\code{statistic}}{the value of the test statistic.}
#'   \item{\code{parameter}}{the degrees of freedom of the test statistic or
#'   NA.} \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{conf.int}}{a confidence interval for the mean appropriate to
#'   the specified alternative hypothesis.} \item{\code{estimate}}{the estimated
#'   mean or difference in means.} \item{\code{null.value}}{the specified
#'   hypothesized value of the mean or mean difference.}
#'   \item{\code{stderr}}{the standard error of the mean (difference).}
#'   \item{\code{alternative}}{a character string describing the alternative
#'   hypothesis.} \item{\code{method}}{a character string indicating what type
#'   of test was performed.} \item{\code{data.name}}{a character string giving
#'   the name(s) of the data.}
#' @import stats
#' @author Kai Li \email{kai.li@stonybrook.edu}
#' @references Kuan P F, Huang B. A simple and robust method for partially
#'   matched samples using the p-values pooling approach. \emph{Statistics in
#'   medicine}. 2013; 32(19): 3247-3259.
#'
#'   Ekbohm G. On comparing means in the paired case with incomplete data on
#'   both responses. \emph{Biometrika}. 1976; 63(2):299-304.
#' @seealso \code{\link[stats]{t.test}} for the one and two-sample t-tests.
#'
#'   \code{\link[stats]{var.test}} for the F test to compare two variances.
#'
#'   \code{\link[stats]{shapiro.test}} for the Shapiro-Wilk test of normality.
#'
#'   \code{\link[stats]{wilcox.test}} for the one and two-sample Wilcoxon tests.
#'
#'   \code{\link{weighted.z}}, \code{\link{modified.t}},
#'   \code{\link{corrected.z}}, and \code{\link{mle.hetero}} for other
#'   statistical approaches for partially matched samples.
#' @examples
#' # pm is a sample dataset for the PMLi package
#'
#' # Ekbohm's MLE-based test under homoscedasticity formula interface
#' mle.homo(pm, "less", conf.level = 0.99)
#'
#' # p-value of Ekbohm's MLE-based test under homoscedasticity
#' p.value1 <- mle.homo(pm)$p.value
#' @export

mle.homo <- function(data, alternative = c("two.sided", "less", "greater"),
                     mu = 0, conf.level = 0.95, ...){
  
  DNAME <- paste(deparse(substitute(data)))
  alternative <- match.arg(alternative)
  if(length(data) == 0){
    stop("The data is empty. Please input a nonempty dataset to begin.\n")
  }else if(is.vector(data)){
    stop("The data is a vector. Please input a nonempty dataset to begin.\n")
  }else if(length(data[1, ]) == 2){
    data <- as.data.frame(t(data))
  }else if(length(data[, 1]) == 2){
    data <- as.data.frame(data)
  }else{
    stop("The data input is not partially matched.\n")
  }
  
  data1 <- data[which(is.na(data[1, ]) == FALSE & is.na(data[2, ]) == FALSE)]
  data2 <- data[which(is.na(data[1, ]) == FALSE & is.na(data[2, ]))]
  data3 <- data[which(is.na(data[1, ]) & is.na(data[2, ]) == FALSE)]
  
  n1 <- length(data1)
  n2 <- length(data2)
  n3 <- length(data3)
  
  if(n1 == 0 & n2 == 0 & n3 == 0){
    stop("The data only contains missing values. Please input a valid dataset
  to begin.\n")
  }else if(n1 == 0 & n2 > 0 & n3 == 0){
    if(shapiro.test(unlist(data2[1, ]))$p.value <= 0.05){
      test <- wilcox.test(unlist(data2[1, ]), alternative = alternative, 
                          mu = mu, paired = FALSE, conf.int = TRUE, 
                          conf.level = conf.level)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data 
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and only the first sample is given. 
  Because the Shapiro-Wilk test is significant, the p-value for the
  given sample will be outputted using the one-sample Wilcoxon 
  signed rank test method.\n")
    }else{
      test <- t.test(data2[1, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and only the first sample is given. 
  Because the Shapiro-Wilk test is not significant, the p-value for
  the given sample will be outputted using the one-sample t-test 
  method.\n")
    }
    STATISTIC <- test$statistic
    PARAMETER <- test$parameter
    PVAL <- test$p.value
    CINT <- c(test$conf.int[1], test$conf.int[2])
    ESTIMATE <- test$estimate
    STDERR <- test$stderr
    METHOD <- test$method
    DNAME <- "data2"
    names(STATISTIC) <- names(test$statistic)
    names(ESTIMATE) <- names(test$estimate)
    names(mu) <- names(test$null.value)
    attr(CINT, "conf.level") <- conf.level
  }else if(n1 == 0 & n2 == 0 & n3 > 0){
    if(shapiro.test(unlist(data3[2, ]))$p.value <= 0.05){
      test <- wilcox.test(unlist(data3[2, ]), alternative = alternative, 
                          mu = mu, paired = FALSE, conf.int = TRUE, 
                          conf.level = conf.level)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and only the second sample is given. 
  Because the Shapiro-Wilk test is significant, the p-value for the
  given sample will be outputted using the one-sample Wilcoxon 
  signed rank test method.\n")
    }else{
      test <- t.test(data3[2, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and only the second sample is given. 
  Because the Shapiro-Wilk test is not significant, the p-value for
  the given sample will be outputted using the one-sample t-test 
  method.\n")
    }
    STATISTIC <- test$statistic
    PARAMETER <- test$parameter
    PVAL <- test$p.value
    CINT <- c(test$conf.int[1], test$conf.int[2])
    ESTIMATE <- test$estimate
    STDERR <- test$stderr
    METHOD <- test$method
    DNAME <- "data3"
    names(STATISTIC) <- names(test$statistic)
    names(ESTIMATE) <- names(test$estimate)
    names(mu) <- names(test$null.value)
    attr(CINT, "conf.level") <- conf.level
  }else if(n1 == 0 & n2 > 0 & n3 > 0){
    if(shapiro.test(unlist(data2[1, ]))$p.value <= 0.05 | 
       shapiro.test(unlist(data3[2, ]))$p.value <= 0.05){
      test <- wilcox.test(unlist(data2[1, ]), unlist(data3[2, ]),
                          alternative = alternative, mu = mu, paired = FALSE, 
                          conf.int = TRUE, conf.level = conf.level)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data 
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and two independent samples are given. 
  Because the Shapiro-Wilk test is significant for at least for one
  of the samples, the p-value for the given independent samples 
  will be outputted using the two-sample Wilcoxon rank sum test 
  method.\n")
    }else{  
      if(var.test(unlist(data2[1, ]), unlist(data3[2, ]))$p.value <= 0.05){
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data does not 
  have any matched pairs, and two independent samples are given. 
  Because the Shapiro-Wilk test is not significant for both 
  independent samples, and the F test to compare the variances of
  the independent samples is significant, the p-value for the 
  given independent samples will be outputted using the 
  two-sample t-test method with unequal variance assumption.\n")
      }else{
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = TRUE)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data does not
  have any matched pairs, and two independent samples are given. 
  Because the Shapiro-Wilk test is not significant for both 
  independent samples, and the F test to compare the variances of
  the independent samples is not significant, the p-value for the
  given independent samples will be outputted using the 
  two-sample t-test method with equal variance assumption.\n")
      }
    }
    STATISTIC <- test$statistic
    PARAMETER <- test$parameter
    PVAL <- test$p.value
    CINT <- c(test$conf.int[1], test$conf.int[2])
    ESTIMATE <- test$estimate
    STDERR <- test$stderr
    METHOD <- test$method
    DNAME <- "data2 + data3"
    names(STATISTIC) <- names(test$statistic)
    names(ESTIMATE) <- names(test$estimate)
    names(mu) <- names(test$null.value)
    attr(CINT, "conf.level") <- conf.level
  }else if(n1 > 0 & n2 == 0 & n3 == 0){
    if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      test <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                          alternative = alternative, mu = mu, paired = TRUE,
                          conf.int = TRUE, conf.level = conf.level)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data is 
  completly matched. Because the Shapiro-Wilk test is significant
  for at least one of the samples, the p-value for the given
  data will be outputted using Wilcoxon signed rank test method
  for matched pairs.\n")
    }else{
      test <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                     alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = TRUE, var.equal = FALSE)
      warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data is 
  completly matched. Because the Shapiro-Wilk test is not 
  significant for both independent samples, the p-value for the 
  given data will be outputted using the matched t-test 
  method.\n")
    }
    STATISTIC <- test$statistic
    PARAMETER <- test$parameter
    PVAL <- test$p.value
    CINT <- c(test$conf.int[1], test$conf.int[2])
    ESTIMATE <- test$estimate
    STDERR <- test$stderr
    METHOD <- test$method
    DNAME <- "data1"
    names(STATISTIC) <- names(test$statistic)
    names(ESTIMATE) <- names(test$estimate)
    names(mu) <- names(test$null.value)
    attr(CINT, "conf.level") <- conf.level
  }else if(n1 > 0 & n2 > 0 & n3 == 0){
    if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      test1 <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                           alternative = alternative, mu = mu, paired = TRUE,
                           conf.int = TRUE, conf.level = conf.level)
      if(shapiro.test(unlist(cbind(data1, data2)[1, ]))$p.value <= 0.05 |
         shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
        test2 <- wilcox.test(unlist(cbind(data1, data2)[1, ]), 
                             unlist(data1[2, ]), alternative = alternative, 
                             mu = mu, paired = FALSE, conf.int = TRUE, 
                             conf.level = conf.level)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent samples
  using the two-sample Wilcoxon rank sum test method will also be
  outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the first sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is significant, the p-value for the complete
  matched pairs using the Wilcoxon signed rank test method
  for matched pairs will be outputted. Because the 
  Shapiro-Wilk test for the two independent samples is not 
  significant and the F test to compare the variances of the
  independent samples is significant, the two independent
  samples using the two-sample t-test method with unequal
  variance assumption will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the first sample
  is given. Because the Shapiro-Wilk test for the matched 
  paired samples is significant, the p-value for the complete
  matched pairs using the Wilcoxon signed rank test method
  for matched pairs will be outputted. Because the 
  Shapiro-Wilk test for the two independent samples is not 
  significant and the F test to compare the variances of the
  independent samples is not significant, the two independent
  samples using the two-sample t-test method with equal 
                     variance assumption will also be outputted.\n")
        }
      }
    }else{
      test1 <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                      alternative = alternative, mu = mu, 
                      conf.level = conf.level, paired = TRUE, 
                      var.equal = FALSE)
      if(shapiro.test(unlist(cbind(data1, data2)[1, ]))$p.value <= 0.05 |
         shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
        test2 <- wilcox.test(unlist(cbind(data1, data2)[1, ]), 
                             unlist(data1[2, ]), alternative = alternative, 
                             mu = mu, paired = FALSE, conf.int = TRUE, 
                             conf.level = conf.level)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples is
  not significant, the p-value for the complete matched pairs
  using the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent samples
  using the two-sample Wilcoxon rank sum test method will also be
  outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the 
  data contains partially matched samples. The inputted data 
  contains complete matched pairs but only the first sample
  is given. Because the Shapiro-Wilk test for the matched 
  paired samples is not significant, the p-value for the 
  complete matched pairs using the matched t-test method for
  matched pairs will be outputted. Because the Shapiro-Wilk
  test for the two independent samples is not significant
  and the F test to compare the variances of the independent
  samples is significant, the two independent samples using
  the two-sample t-test method with unequal variance
  assumption will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the 
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the first sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is not significant, the p-value for the
  complete matched pairs using the matched t-test method for
  matched pairs will be outputted. Because the Shapiro-Wilk
  test for the two independent samples is not significant
  and the F test to compare the variances of the independent
  samples is not significant, the two independent samples
  using the two-sample t-test method with equal variance
  assumption will also be outputted.\n")
        }
      }
    }
    STATISTIC <- c(test1$statistic, test2$statistic)
    PARAMETER <- c(test1$parameter, test2$parameter)
    PVAL <- c(test1$p.value, test2$p.value)
    CINT <- c(paste(test1$conf.int[1], test1$conf.int[2], "and"),
              paste(test2$conf.int[1], test2$conf.int[2]))
    ESTIMATE <- c(test1$estimate, test2$estimate)
    STDERR <- c(test1$stderr, test2$stderr)
    METHOD <- paste(test1$method, "and", test2$method)
    DNAME <- "data1 and data1 + data2"
    names(STATISTIC) <- c(names(test1$statistic), names(test2$statistic))
    names(ESTIMATE) <- c(names(test1$estimate), names(test2$estimate))
    names(mu) <- paste(names(test1$null.value), "and", 
                       names(test2$null.value))
    attr(CINT, "conf.level") <- conf.level
  }else if(n1 > 0 & n2 == 0 & n3 > 0){
    if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      test1 <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                           alternative = alternative, mu = mu, paired = TRUE,
                           conf.int = TRUE, conf.level = conf.level)
      if(shapiro.test(unlist(cbind(data1, data3)[2, ]))$p.value <= 0.05 |
         shapiro.test(unlist(data1[1, ]))$p.value <= 0.05){
        test2 <- wilcox.test(unlist(cbind(data1, data3)[2, ]), 
                             unlist(data1[1, ]), alternative = alternative, 
                             mu = mu, paired = FALSE, conf.int = TRUE, 
                             conf.level = conf.level)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent samples
  using the two-sample Wilcoxon rank sum test method will also be
  outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the second sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is significant, the p-value for the complete
  matched pairs using the Wilcoxon signed rank test method
  for matched pairs will be outputted. Because the 
  Shapiro-Wilk test for the two independent samples is not
  significant and the F test to compare the variances of the
  independent samples is significant, the two independent
  samples using the two-sample t-test method with unequal
  variance assumption will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the second sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is significant, the p-value for the complete
  matched pairs using the Wilcoxon signed rank test method
  for matched pairs will be outputted. Because the 
  Shapiro-Wilk test for the two independent samples is not
  significant and the F test to compare the variances of the
  independent samples is not significant, the two independent
  samples using the two-sample t-test method with equal
  variance assumption will also be outputted.\n")
        }
      }
    }else{
      test1 <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                      alternative = alternative, mu = mu, 
                      conf.level = conf.level, paired = TRUE, 
                      var.equal = FALSE)
      if(shapiro.test(unlist(cbind(data1, data3)[2, ]))$p.value <= 0.05 |
         shapiro.test(unlist(data1[1, ]))$p.value <= 0.05){
        test2 <- wilcox.test(unlist(cbind(data1, data3)[2, ]), 
                             unlist(data1[1, ]), alternative = alternative, 
                             mu = mu, paired = FALSE, conf.int = TRUE, 
                             conf.level = conf.level)
        warning("Ekbohm's MLE-based test under homoscedasticity assumes the data
  contains partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples is
  not significant, the p-value for the complete matched pairs
  using the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent samples
  using the two-sample Wilcoxon rank sum test method will also be
  outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the 
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the second sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is not significant, the p-value for the
  complete matched pairs using the matched t-test method for
  matched pairs will be outputted. Because the Shapiro-Wilk
  test for the two independent samples is not significant
  and the F test to compare the variances of the independent
  samples is significant, the two independent samples using
  the two-sample t-test method with unequal variance
  assumption will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Ekbohm's MLE-based test under homoscedasticity assumes the 
  data contains partially matched samples. The inputted data
  contains complete matched pairs but only the second sample
  is given. Because the Shapiro-Wilk test for the matched
  paired samples is not significant, the p-value for the
  complete matched pairs using the matched t-test method for
  matched pairs will be outputted. Because the Shapiro-Wilk
  test for the two independent samples is not significant
  and the F test to compare the variances of the independent
  samples is not significant, the two independent samples
  using the two-sample t-test method with equal variance
  assumption will also be outputted.\n")
        }
      }
    }
    STATISTIC <- c(test1$statistic, test2$statistic)
    PARAMETER <- c(test1$parameter, test2$parameter)
    PVAL <- c(test1$p.value, test2$p.value)
    CINT <- c(paste(test1$conf.int[1], test1$conf.int[2], "and"),
              paste(test2$conf.int[1], test2$conf.int[2]))
    ESTIMATE <- c(test1$estimate, test2$estimate)
    STDERR <- c(test1$stderr, test2$stderr)
    METHOD <- paste(test1$method, "and", test2$method)
    DNAME <- "data1 and data1 + data3"
    names(STATISTIC) <- c(names(test1$statistic), names(test2$statistic))
    names(ESTIMATE) <- c(names(test1$estimate), names(test2$estimate))
    names(mu) <- paste(names(test1$null.value), "and", 
                       names(test2$null.value))
    attr(CINT, "conf.level") <- conf.level
  }else{
    s_t1 <- sd(unlist(data1[1,]))
    s_n1 <- sd(unlist(data1[2,]))
    s_tn1 <- cov(unlist(data1[1,]), unlist(data1[2,]))
    r <- s_tn1/(s_t1*s_n1)
    s_t <- sd(unlist(data2[1,]))
    s_n <- sd(unlist(data3[2,]))
    sigma_hat_squared <- (s_t1^2*(n1-1)+s_n1^2*(n1-1)+(1+r^2)*(s_t^2*(n2-1)+
                                    s_n^2*(n3-1)))/(2*(n1-1)+(1+r^2)*(n2+n3-2))
    f_star <- n1*(n1+n3+n2*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    g_star <- n1*(n1+n2+n3*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    v1_star <- sigma_hat_squared*(2*n1*(1-r)+(n2+n3)*(1-r^2))/
      ((n1+n2)*(n1+n3)-n2*n3*r^2)
    t1_bar <- mean(unlist(data1[1,]))
    t_bar <- mean(unlist(data2[1,]))
    n1_bar <- mean(unlist(data1[2,]))
    n_bar <- mean(unlist(data3[2,]))
    
    z_e <- (f_star*(t1_bar-t_bar)-g_star*(n1_bar-n_bar)+t_bar-n_bar-mu)/
      sqrt(v1_star)
    
    if(alternative == "two.sided"){
      p_value <- 2*(1-pt(abs(z_e), n1))
      cint <- qt(1-(1-conf.level)/2, n1)
      cint <- z_e + c(-cint, cint)
    }else if(alternative == "less"){
      p_value <- pt(z_e, n1)
      cint <- c(-Inf, z_e+qt(conf.level, n1))
    }else{
      p_value <- 1-pt(z_e, n1)
      cint <- c(z_e-qt(conf.level, n1), Inf)
    }
    
    STATISTIC <- z_e
    names(STATISTIC) <- "Z_E"
    PARAMETER <- c("df" = n1)
    PVAL <- p_value
    CINT <- mu+cint*sqrt(v1_star)
    attr(CINT, "conf.level") <- conf.level
    ESTIMATE <- mu+z_e*sqrt(v1_star)
    names(ESTIMATE) <- names(mu) <- "difference in means"
    STDERR <- sqrt(v1_star)
    METHOD = "Ekbohm's MLE-based test under homoscedasticity"
  }
  
  result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL,
                 conf.int = CINT, estimate = ESTIMATE, null.value = mu, 
                 stderr = STDERR, alternative = alternative, method = METHOD, 
                 data.name = DNAME)
  attr(result, "class") <- "htest"
  return(result)
}
