#' @title Liptak's Weighted Z-Test
#' @description \code{weighted.z} is used to perform Liptak's weighted Z-test
#'   for partially matched samples, specified by giving a data frame or matrix,
#'   testing hypothesis, true mean, the corresponding weights in Liptak's
#'   weighted Z-test and confidence level.
#' @param data a 2-by-n or n-by-2 partially matched samples data frame or
#'   matrix.
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" or "less". You can specify
#'   just the initial letter.
#' @param mu a number indicating the true difference in means.
#' @param w1 a number indicating the weight for the first p-value. If the weight
#'   is not inputted, a default weight is assigned to to be the square root of
#'   two times the paired sample size.
#' @param w2 a number indicating the weight for the second p-value. If the
#'   weight is not inputted, a default weight is assigned to to be the square
#'   root of the sum of the two independent samples size.
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
#'   Liptak T. On the combination of independent tests. \emph{Magyar Tudom
#'   Aanyos Akad Aemia Matematikai Kutat Ao Intezetenek Kozlemenyei}. 1958;
#'   3:171-197.
#' @seealso \code{\link[stats]{t.test}} for the one and two-sample t-tests.
#'
#'   \code{\link[stats]{var.test}} for the F test to compare two variances.
#'
#'   \code{\link[stats]{shapiro.test}} for the Shapiro-Wilk test of normality.
#'
#'   \code{\link[stats]{wilcox.test}} for the one and two-sample Wilcoxon tests.
#'
#'   \code{\link{modified.t}}, \code{\link{corrected.z}},
#'   \code{\link{mle.hetero}}, and \code{\link{mle.homo}} for other statistical
#'   approaches for partially matched samples.
#' @examples
#' # pm is a sample dataset for the PMLi package
#'
#' # Liptak's weighted Z-test formula interface
#' weighted.z(pm, "less", conf.level = 0.99)
#'
#' # p-value of Liptak's weighted Z-test
#' p.value1 <- weighted.z(pm)$p.value
#' @export

weighted.z <- function(data, alternative = c("two.sided", "less", "greater"),
                       mu = 0, w1 = NULL, w2 = NULL, conf.level = 0.95, ...){
  
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
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and only the first sample is given. Because the
  Shapiro-Wilk test is significant, the p-value for the given
  sample will be outputted using the one-sample Wilcoxon signed
  rank test method.\n")
    }else{
      test <- t.test(data2[1, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and only the first sample is given. Because the
  Shapiro-Wilk test is not significant, the p-value for the given
  sample will be outputted using the one-sample t-test method.\n")
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
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and only the second sample is given. Because the
  Shapiro-Wilk test is significant, the p-value for the given
  sample will be outputted using the one-sample Wilcoxon signed
  rank test method.\n")
    }else{
      test <- t.test(data3[2, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and only the second sample is given. Because the
  Shapiro-Wilk test is not significant, the p-value for the given
  sample will be outputted using the one-sample t-test method.\n")
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
      warning("Liptak's weighted Z-test assumes the data contains partially
               matched samples. The inputted data does not have any matched
               pairs, and two independent samples are given. Because the
               Shapiro-Wilk test is significant for at least for one of the
               samples, the p-value for the given independent samples will be
               outputted using the two-sample Wilcoxon rank sum test method.\n")
    }else{  
      if(var.test(unlist(data2[1, ]), unlist(data3[2, ]))$p.value <= 0.05){
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE)
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and two independent samples are given. Because the
  Shapiro-Wilk test is not significant for both independent
  samples, and the F test to compare the variances of the
  independent samples is significant, the p-value for the given
  independent samples will be outputted using the two-sample
  t-test method with unequal variance assumption.\n")
      }else{
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = TRUE)
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data does not have any matched
  pairs, and two independent samples are given. Because the
  Shapiro-Wilk test is not significant for both independent
  samples, and the F test to compare the variances of the
  independent samples is not significant, the p-value for the
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
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data is completly matched.
  Because the Shapiro-Wilk test is significant for at least one
  of the samples, the p-value for the given data will be
  outputted using Wilcoxon signed rank test method for matched
  pairs.\n")
    }else{
      test <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                     alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = TRUE, var.equal = FALSE)
      warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data is completly matched. Because
  the Shapiro-Wilk test is not significant for both independent
  samples, the p-value for the given data will be outputted using
  the matched t-test method.\n")
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
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent
  samples using the two-sample Wilcoxon rank sum test method
  will also be outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is
  significant, the two independent samples using the two-sample
  t-test method with unequal variance assumption will also be
  outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the two-sample
  t-test method with equal variance assumption will also be
  outputted.\n")
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
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent
  samples using the two-sample Wilcoxon rank sum test method
  will also be outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is 
  significant, the two independent samples using the
  two-sample t-test method with unequal variance assumption
  will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the first sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the
  two-sample t-test method with equal variance assumption
  will also be outputted.\n")
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
    names(STATISTIC) <- c(names(test1$statistic), 
                          names(test2$statistic))
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
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the 
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent
  samples using the two-sample Wilcoxon rank sum test method
  will also be outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is
  significant, the two independent samples using the two-sample
  t-test method with unequal variance assumption will also be
  outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is
  significant, the p-value for the complete matched pairs using
  the Wilcoxon signed rank test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the two-sample
  t-test method with equal variance assumption will also be
  outputted.\n")
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
        warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be 
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is significant, the two independent
  samples using the two-sample Wilcoxon rank sum test method
  will also be outputted.\n")
      }else{
        if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
           p.value <= 0.05){
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = FALSE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is
  significant, the two independent samples using the 
  two-sample t-test method with unequal variance assumption
  will also be outputted.\n")
        }else{
          test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                          alternative = alternative, mu = mu, 
                          conf.level = conf.level, paired = FALSE, 
                          var.equal = TRUE)
          warning("Liptak's weighted Z-test assumes the data contains partially
  matched samples. The inputted data contains complete matched
  pairs but only the second sample is given. Because the
  Shapiro-Wilk test for the matched paired samples is not
  significant, the p-value for the complete matched pairs using
  the matched t-test method for matched pairs will be
  outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the
  two-sample t-test method with equal variance assumption
  will also be outputted.\n")
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
    names(STATISTIC) <- c(names(test1$statistic), 
                          names(test2$statistic))
    names(ESTIMATE) <- c(names(test1$estimate), names(test2$estimate))
    names(mu) <- paste(names(test1$null.value), "and", 
                       names(test2$null.value))
    attr(CINT, "conf.level") <- conf.level
  }else{
    if(alternative == "two.sided"){
      alternative1 <- "greater"
    }else{
      alternative1 <- alternative
    }
    if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |  
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      p_1i <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                          alternative = alternative1, mu = mu, paired = TRUE,
                          conf.int = TRUE, conf.level = conf.level)$p.value
    }else{
      p_1i <- t.test(data1[1, ], data1[2, ], alternative = alternative1,
                     mu = mu, conf.level = conf.level, paired = TRUE, 
                     var.equal = FALSE)$p.value
    }
    
    if(shapiro.test(unlist(data2[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data3[2, ]))$p.value <= 0.05){
      p_2i <- wilcox.test(unlist(data2[1, ]), unlist(data3[2, ]),
                          alternative = alternative1, mu = mu, paired = FALSE,
                          conf.int = TRUE, conf.level = conf.level)$p.value
    }else{
      if(var.test(unlist(data2[1, ]), unlist(data3[2, ]))$p.value <= 0.05){
        p_2i <- t.test(data2[1, ], data3[2, ], alternative = alternative1,
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE)$p.value
      }else{
        p_2i <- t.test(data2[1, ], data3[2, ], alternative = alternative1,
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = TRUE)$p.value
      }
    }
    
    if(is.null(w1) == TRUE){
      w1 <- sqrt(2*n1)
    }
    if(is.null(w2) == TRUE){
      w2 <- sqrt(n2+n3)
    }
    z_1i <- qnorm(1-p_1i)
    z_2i <- qnorm(1-p_2i)
    
    p_ci <- 1-pnorm((w1*z_1i+w2*z_2i)/sqrt(w1^2+w2^2))
    if(alternative == "two.sided"){
      if(p_ci < 1/2){
        p_ci_star <- 2*p_ci
        p_value <- p_ci_star
      }else{
        p_ci_star <- 2*(1-p_ci)
        p_value <- p_ci_star
      }
    }else{
      p_value <- p_ci
    }
    
    STATISTIC <- c(z_1i, z_2i)
    names(STATISTIC) <- c("Z_1i", "Z_2i")
    PARAMETER <- NULL
    PVAL <- p_value
    CINT <- NULL
    ESTIMATE <- NULL
    names(ESTIMATE) <- NULL
    names(mu) <- "difference in means"
    STDERR <- NULL
    METHOD = "Liptak's weighted Z-test"
  }

result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL,
               conf.int = CINT, estimate = ESTIMATE, null.value = mu, 
               stderr = STDERR, alternative = alternative, method = METHOD, 
               data.name = DNAME)
attr(result, "class") <- "htest"
return(result)
}
