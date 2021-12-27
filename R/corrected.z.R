#' @title Looney and Jones's Corrected Z-Test
#' @description \code{corrected.z} is used to perform Looney and Jones's
#'   corrected Z-test for partially matched samples, specified by giving a data
#'   frame or matrix, testing hypothesis, true mean and confidence level.
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
#'   Looney S, Jones P. A method for comparing two normal means using combined
#'   samples of correlated and uncorrelated data. \emph{Statistics in Medicine}.
#'   2003; 22(9):1601-1610. [PubMed: 12704618]
#' @seealso \code{\link[stats]{t.test}} for the one and two-sample t-tests.
#'
#'   \code{\link[stats]{var.test}} for the F test to compare two variances.
#'
#'   \code{\link[stats]{shapiro.test}} for the Shapiro-Wilk test of normality.
#'
#'   \code{\link[stats]{wilcox.test}} for the one and two-sample Wilcoxon tests.
#'
#'   \code{\link{weighted.z}}, \code{\link{modified.t}},
#'   \code{\link{mle.hetero}}, and \code{\link{mle.homo}} for other statistical
#'   approaches for partially matched samples.
#' @examples
#' # pm is a sample dataset for the PMLi package
#'
#' # Looney and Jones's corrected Z-test formula interface
#' corrected.z(pm, "less", conf.level = 0.99)
#'
#' # p-value of Looney and Jones's corrected Z-test
#' p.value1 <- corrected.z(pm)$p.value
#' @export

corrected.z <- function(data, alternative = c("two.sided", "less", "greater"),
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
      warning("Looney and Jones's corrected Z-test assumes the data contains 
  partially matched samples. The inputted data does not have any 
  matched pairs, and only the first sample is given. Because the 
  Shapiro-Wilk test is significant, the p-value for the given 
  sample will be outputted using the one-sample Wilcoxon signed 
  rank test method.\n")
    }else{
      test <- t.test(data2[1, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data does not have any 
  matched pairs, and only the first sample is given. Because the
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
      warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data does not have any
  matched pairs, and only the second sample is given. Because the
  Shapiro-Wilk test is significant, the p-value for the given
  sample will be outputted using the one-sample Wilcoxon signed
  rank test method.\n")
    }else{
      test <- t.test(data3[2, ], alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = FALSE, 
                     var.equal = FALSE)
      warning("Looney and Jones's corrected Z-test assumes the data contains 
  partially matched samples. The inputted data does not have any
  matched pairs, and only the second sample is given. Because the
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
      warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data does not have any
  matched pairs, and two independent samples are given. Because the
  Shapiro-Wilk test is significant for at least for one of the
  samples, the p-value for the given independent samples will be
  outputted using the two-sample Wilcoxon rank sum test method.\n")
    }else{  
      if(var.test(unlist(data2[1, ]), unlist(data3[2, ]))$p.value <= 0.05){
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE)
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data does not have any
  matched pairs, and two independent samples are given. Because
  the Shapiro-Wilk test is not significant for both independent
  samples, and the F test to compare the variances of the
  independent samples is significant, the p-value for the given
  independent samples will be outputted using the two-sample
  t-test method with unequal variance assumption.\n")
      }else{
        test <- t.test(data2[1, ], data3[2, ], alternative = alternative, 
                       mu = mu, conf.level = conf.level, paired = FALSE, 
                       var.equal = TRUE)
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data does not have any
  matched pairs, and two independent samples are given. Because
  the Shapiro-Wilk test is not significant for both independent
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
      warning("Looney and Jones's corrected Z-test assumes the data contains 
  partially matched samples. The inputted data is completly
  matched. Because the Shapiro-Wilk test is significant for at
  least one of the samples, the p-value for the given data
  will be outputted using Wilcoxon signed rank test method for
  matched pairs.\n")
    }else{
      test <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                     alternative = alternative, mu = mu, 
                     conf.level = conf.level, paired = TRUE, var.equal = FALSE)
      warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data is completly
  matched. Because the Shapiro-Wilk test is not significant for
  both independent samples, the p-value for the given data will
  be outputted using the matched t-test method.\n")
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
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
                      conf.level = conf.level, paired = TRUE, var.equal = FALSE)
      if(shapiro.test(unlist(cbind(data1, data2)[1, ]))$p.value <= 0.05 |
         shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
        test2 <- wilcox.test(unlist(cbind(data1, data2)[1, ]), 
                             unlist(data1[2, ]), alternative = alternative, 
                             mu = mu, paired = FALSE, conf.int = TRUE, 
                             conf.level = conf.level)
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the first sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the two-sample
  t-test method with equal variance assumption will also be
  outputted.\n")
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
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is significant, the p-value for the complete matched pairs
  using the Wilcoxon signed rank test method for matched pairs
  will be outputted. Because the Shapiro-Wilk test for the two
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
        warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
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
          warning("Looney and Jones's corrected Z-test assumes the data contains
  partially matched samples. The inputted data contains
  complete matched pairs but only the second sample is given.
  Because the Shapiro-Wilk test for the matched paired samples
  is not significant, the p-value for the complete matched
  pairs using the matched t-test method for matched pairs will
  be outputted. Because the Shapiro-Wilk test for the two
  independent samples is not significant and the F test to
  compare the variances of the independent samples is not
  significant, the two independent samples using the two-sample
  t-test method with equal variance assumption will also be
  outputted.\n")
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
    t_bar_star <- mean(unlist(cbind(data1[1,], data2[1,])))
    n_bar_star <- mean(unlist(cbind(data1[2,], data3[2,])))
    s_t_star <- sd(unlist(cbind(data1[1,], data2[1,])))
    s_n_star <- sd(unlist(cbind(data1[2,], data3[2,])))
    s_tn1 <- cov(unlist(data1[1,]), unlist(data1[2,]))
    
    z_corr <- (t_bar_star-n_bar_star-mu)/sqrt(s_t_star^2/(n1+n2)+s_n_star^2/
                                          (n1+n3)-2*n1*s_tn1/((n1+n2)*(n1+n3)))
    if(alternative == "two.sided"){
      p_value <- 2*(1- pnorm(abs(z_corr)))
      cint <- qnorm(1-(1-conf.level)/2)
      cint <- z_corr + c(-cint, cint)
    }else if(alternative == "less"){
      p_value <- pnorm(z_corr)
      cint <- c(-Inf, z_corr+qnorm(conf.level))
    }else{
      p_value <- 1-pnorm(z_corr)
      cint <- c(z_corr-qnorm(conf.level), Inf)
    }
    
    STATISTIC <- z_corr
    names(STATISTIC) <- "Z_corr"
    PARAMETER <- NULL
    PVAL <- p_value
    CINT <- mu+cint*sqrt(s_t_star^2/(n1+n2)+s_n_star^2/
                           (n1+n3)-2*n1*s_tn1/((n1+n2)*(n1+n3)))
    attr(CINT, "conf.level") <- conf.level
    ESTIMATE <- mu+z_corr*sqrt(s_t_star^2/(n1+n2)+s_n_star^2/
                                 (n1+n3)-2*n1*s_tn1/((n1+n2)*(n1+n3)))
    names(ESTIMATE) <- names(mu) <- "difference in means"
    STDERR <- sqrt(s_t_star^2/(n1+n2)+s_n_star^2/
                     (n1+n3)-2*n1*s_tn1/((n1+n2)*(n1+n3)))
    METHOD = "Looney and Jones's corrected Z-test"
  }
  
  result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL,
                 conf.int = CINT, estimate = ESTIMATE, null.value = mu, 
                 stderr = STDERR, alternative = alternative, method = METHOD, 
                 data.name = DNAME)
  attr(result, "class") <- "htest"
  return(result)
}
