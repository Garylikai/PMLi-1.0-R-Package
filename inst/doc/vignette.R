## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  
)

## ----message=FALSE, warning=FALSE---------------------------------------------
install.packages("PMLi_1.0.tar.gz", type = "source")
library(PMLi)

## -----------------------------------------------------------------------------
data <- pm

## -----------------------------------------------------------------------------
if(length(data) == 0){
  stop("The data is empty. Please input a nonempty dataset to begin.\n")
}

## -----------------------------------------------------------------------------
if(is.vector(data)){
    stop("The data is a vector. Please input a nonempty dataset to begin.\n")
}else if(length(data[1, ]) == 2){
  data <- as.data.frame(t(data))
}else if(length(data[, 1]) == 2){
  data <- as.data.frame(data)
}else{
  stop("The data inputted is not partially matched.\n")
}

## -----------------------------------------------------------------------------
data1 <- data[which(is.na(data[1, ]) == FALSE & is.na(data[2, ]) == FALSE)]
data2 <- data[which(is.na(data[1, ]) == FALSE & is.na(data[2, ]))]
data3 <- data[which(is.na(data[1, ]) & is.na(data[2, ]) == FALSE)]

n1 <- length(data1)
n2 <- length(data2)
n3 <- length(data3)

## -----------------------------------------------------------------------------
if(n1 == 0 & n2 == 0 & n3 == 0){
    stop("The data only contains missing values. Please input a valid dataset
          to begin.\n")
}

## -----------------------------------------------------------------------------
if(n1 == 0 & n2 > 0 & n3 == 0){
  if(shapiro.test(unlist(data2[1, ]))$p.value <= 0.05){
    test <- wilcox.test(unlist(data2[1, ]), alternative = alternative, mu = mu,
                        paired = FALSE, conf.int = TRUE, conf.level = conf.level)
  }else{
    test <- t.test(data2[1, ], alternative = alternative, mu = mu, 
                   conf.level = conf.level, paired = FALSE, var.equal = FALSE)
  }
}

## -----------------------------------------------------------------------------
if(n1 == 0 & n2 == 0 & n3 > 0){
  if(shapiro.test(unlist(data3[2, ]))$p.value <= 0.05){
    test <- wilcox.test(unlist(data3[2, ]), alternative = alternative, mu = mu,
                        paired = FALSE, conf.int = TRUE, conf.level = conf.level)
  }else{
    test <- t.test(data3[2, ], alternative = alternative, mu = mu, 
                   conf.level = conf.level, paired = FALSE, var.equal = FALSE)
  }
}

## -----------------------------------------------------------------------------
if(n1 == 0 & n2 > 0 & n3 > 0){
  if(shapiro.test(unlist(data2[1, ]))$p.value <= 0.05 | 
       shapiro.test(unlist(data3[2, ]))$p.value <= 0.05){
      test <- wilcox.test(unlist(data2[1, ]), unlist(data3[2, ]),
                          alternative = alternative, mu = mu, paired = FALSE, 
                          conf.int = TRUE, conf.level = conf.level)
  }else{
    if(var.test(unlist(data2[1, ]), unlist(data3[2, ]))$p.value <= 0.05){
      test <- t.test(data2[1, ], data3[2, ], alternative = alternative, mu = mu,
                     conf.level = conf.level, paired = FALSE, var.equal = FALSE)
    }else{
      test <- t.test(data2[1, ], data3[2, ], alternative = alternative, mu = mu,
                     conf.level = conf.level, paired = FALSE, var.equal = TRUE)
    }
  }
}

## -----------------------------------------------------------------------------
if(length(data) == n1){
  if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
     shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
    test <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                        alternative = alternative, mu = mu, paired = TRUE,
                        conf.int = TRUE, conf.level = conf.level)
  }else{
    test <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                   alternative = alternative, mu = mu, conf.level = conf.level,
                   paired = TRUE, var.equal = FALSE)
  }
}

## -----------------------------------------------------------------------------
if(n1 > 0 & n2 > 0 & n3 == 0){
  if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
     shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
    test1 <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                         alternative = alternative, mu = mu, paired = TRUE,
                         conf.int = TRUE, conf.level = conf.level)
    if(shapiro.test(unlist(cbind(data1, data2)[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      test2 <- wilcox.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]),
                           alternative = alternative, mu = mu, paired = FALSE,
                           conf.int = TRUE, conf.level = conf.level)  
    }else{
      if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
         p.value <= 0.05){
        test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE,
                        var.equal = FALSE)
      }else{
        test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE, 
                        var.equal = TRUE)
      }
    }
  }else{
    test1 <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                    alternative = alternative, mu = mu, conf.level = conf.level,
                    paired = TRUE, var.equal = FALSE)
    if(shapiro.test(unlist(cbind(data1, data2)[1, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
      test2 <- wilcox.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]),
                           alternative = alternative, mu = mu, paired = FALSE,
                           conf.int = TRUE, conf.level = conf.level)
    }else{
      if(var.test(unlist(cbind(data1, data2)[1, ]), unlist(data1[2, ]))$
         p.value <= 0.05){
        test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE, 
                        var.equal = FALSE)
      }else{
        test2 <- t.test(cbind(data1, data2)[1, ], data1[2, ],
                            alternative = alternative, mu = mu, 
                            conf.level = conf.level, paired = FALSE, 
                            var.equal = TRUE)
      }
    }
  }
}

## -----------------------------------------------------------------------------
if(n1 > 0 & n2 == 0 & n3 > 0){
  if(shapiro.test(unlist(data1[1, ]))$p.value <= 0.05 |
     shapiro.test(unlist(data1[2, ]))$p.value <= 0.05){
    test1 <- wilcox.test(unlist(data1[1, ]), unlist(data1[2, ]),
                         alternative = alternative, mu = mu, paired = TRUE,
                         conf.int = TRUE, conf.level = conf.level)
    if(shapiro.test(unlist(cbind(data1, data3)[2, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[1, ]))$p.value <= 0.05){
      test2 <- wilcox.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]),
                           alternative = alternative, mu = mu, paired = FALSE,
                           conf.int = TRUE, conf.level = conf.level)  
    }else{
      if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
         p.value <= 0.05){
        test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE,
                        var.equal = FALSE)
      }else{
        test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE, 
                        var.equal = TRUE)
      }
    }
  }else{
    test1 <- t.test(unlist(data1[1, ]), unlist(data1[2, ]), 
                    alternative = alternative, mu = mu, conf.level = conf.level,
                    paired = TRUE, var.equal = FALSE)
    if(shapiro.test(unlist(cbind(data1, data3)[2, ]))$p.value <= 0.05 |
       shapiro.test(unlist(data1[1, ]))$p.value <= 0.05){
      test2 <- wilcox.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]),
                           alternative = alternative, mu = mu, paired = FALSE,
                           conf.int = TRUE, conf.level = conf.level)
    }else{
      if(var.test(unlist(cbind(data1, data3)[2, ]), unlist(data1[1, ]))$
         p.value <= 0.05){
        test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                        alternative = alternative, mu = mu, 
                        conf.level = conf.level, paired = FALSE, 
                        var.equal = FALSE)
      }else{
        test2 <- t.test(cbind(data1, data3)[2, ], data1[1, ],
                            alternative = alternative, mu = mu, 
                            conf.level = conf.level, paired = FALSE, 
                            var.equal = TRUE)
      }
    }
  }
}

## -----------------------------------------------------------------------------
# Run the following command if you have questions about weighted.z() function.
# ?weighted.z()

## -----------------------------------------------------------------------------
# p-value is returned.
weighted.z(pm)$p.value

## -----------------------------------------------------------------------------
# Formula interface is outputted.
modified.t(pm)

## -----------------------------------------------------------------------------
# Test statistic is returned.
corrected.z(pm, "greater")$statistic

## -----------------------------------------------------------------------------
# Confidence interval is returned with specified arguments.
mle.hetero(pm, alternative = "less", conf.level = 0.90)$conf.int

## -----------------------------------------------------------------------------
# Formula interface is returned with specified arguments.
mle.homo(pm, alternative = "greater", mu = 0.01, conf.level = 0.99)

