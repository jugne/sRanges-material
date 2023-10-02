######################################################
######################################################
# 
######################################################
######################################################
library(ggplot2)
library(coda)
options(digits=7)

# clear workspace
rm(list = ls())

test_vals <- function(prob, true=0, estimated=c()){
    lower <- HPDinterval(as.mcmc(estimated), prob=prob)[1]
    upper <- HPDinterval(as.mcmc(estimated), prob=prob)[2]
    test <- c(as.numeric(true >= lower & true <= upper))
  return(test)
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Set the directory to the directory of the file
wd<-"/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_priors/"
setwd(wd)
dir.create("figures")

# read in the true rates
true_rates_file <- "/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_priors/true_rates.csv"
true.rates <- read.table(true_rates_file, header=TRUE, sep=",")

n<-200
remove_rows<-c()

diversificationRate = data.frame(true=numeric(n), estimated=numeric(n),
                                 upper=numeric(n), lower=numeric(n),
                                 rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                                 rel.err.mode=numeric(n), cv=numeric(n),
                                 hpd.lower=numeric(n), hpd.upper=numeric(n),
                                 hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                 test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                 test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                 test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                 test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                 test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                 test_90=numeric(n), test_95=numeric(n),
                                 test_100=numeric(n))

turnover = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                      lower=numeric(n), rel.err.meadian=numeric(n),
                      rel.err.mean=numeric(n), rel.err.mode=numeric(n),
                      cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n),
                      hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                      test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                      test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                      test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                      test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                      test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                      test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingProportion = data.frame(true=numeric(n), estimated=numeric(n),
                                upper=numeric(n), lower=numeric(n),
                                rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                                rel.err.mode=numeric(n), cv=numeric(n),
                                hpd.lower=numeric(n), hpd.upper=numeric(n),
                                hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))


samplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n),
                                   upper=numeric(n), lower=numeric(n),
                                   rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                                   rel.err.mode=numeric(n), cv=numeric(n),
                                   hpd.lower=numeric(n), hpd.upper=numeric(n),
                                   hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                   test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                   test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                   test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                   test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                   test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                   test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

origin = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                    lower=numeric(n), rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                    rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n),
                    hpd.upper=numeric(n), hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                    test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                    test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                    test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                    test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                    test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                    test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

mrca = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                  lower=numeric(n), rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                  rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n),
                  hpd.upper=numeric(n), hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                  test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                  test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                  test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                  test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                  test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                  test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

rsamplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n),
                                    upper=numeric(n), lower=numeric(n),
                                    rel.err.meadian=numeric(n), rel.err.mean=numeric(n),
                                    rel.err.mode=numeric(n), cv=numeric(n),
                                    hpd.lower=numeric(n), hpd.upper=numeric(n),
                                    hpd.rel.width=numeric(n), test_0=numeric(n),
                                    test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                    test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                    test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                    test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                    test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                    test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                    test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))



post_ess<- numeric(n)
for (i in 182:n){
  
  log <- list.files(path=paste0("run_",i,"/inf"), pattern=paste0("sRanges.",i,".log"), full.names = TRUE)
  
  # read in the log file
  if (length(log)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  
  used.rates = true.rates[which(true.rates$X==i),];
  
  t <- read.table(log[[1]], header=TRUE, sep="\t")

  # take a 10% burnin
  t <- t[-seq(1,ceiling(length(t$diversificationRate)/10)), ]
  
  if (nrow(t)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  ess<-  effectiveSize(as.mcmc(t))
  post_ess[i] = as.numeric(ess["posterior"])
  
  if(post_ess[i]<100){
    remove_rows = cbind(remove_rows, i)
    next
  }
  
  
  diversificationRate$true[i] <- used.rates$div_rate
  diversificationRate$estimated[i] <- median(t$diversificationRate)
  diversificationRate$upper[i] <- quantile(t$diversificationRate,0.975)
  diversificationRate$lower[i] <- quantile(t$diversificationRate,0.025)
  
  
  diversificationRate$rel.err.meadian[i] <- abs(median(t$diversificationRate)-used.rates$div_rate)/used.rates$div_rate
  diversificationRate$rel.err.mean[i] <-abs(mean(t$diversificationRate)-used.rates$div_rate)/used.rates$div_rate
  diversificationRate$rel.err.mode[i] <- abs(Mode(t$diversificationRate)-used.rates$div_rate)/used.rates$div_rate
  diversificationRate$cv[i] <- sqrt(exp(sd(log(t$diversificationRate))**2)-1)
  diversificationRate$hpd.lower[i] <- HPDinterval(as.mcmc(t$diversificationRate))[1]
  diversificationRate$hpd.upper[i] <- HPDinterval(as.mcmc(t$diversificationRate))[2]
  diversificationRate$hpd.rel.width[i] <- (diversificationRate$hpd.upper[i]-diversificationRate$hpd.lower[i])/used.rates$div_rate
  diversificationRate[i,(ncol(diversificationRate)-21+1):ncol(diversificationRate)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate$true[i], t$diversificationRate)
  
  
  turnover$true[i] <- used.rates$turnover
  turnover$estimated[i] <- median(t$turnover)
  turnover$upper[i] <- quantile(t$turnover,0.975)
  turnover$lower[i] <- quantile(t$turnover,0.025)
  
  
  turnover$rel.err.meadian[i] <- abs(median(t$turnover)-used.rates$turnover)/used.rates$turnover
  turnover$rel.err.mean[i] <-abs(mean(t$turnover)-used.rates$turnover)/used.rates$turnover
  turnover$rel.err.mode[i] <- abs(Mode(t$turnover)-used.rates$turnover)/used.rates$turnover
  turnover$cv[i] <- sqrt(exp(sd(log(t$turnover))**2)-1)
  turnover$hpd.lower[i] <- HPDinterval(as.mcmc(t$turnover))[1]
  turnover$hpd.upper[i] <- HPDinterval(as.mcmc(t$turnover))[2]
  # turnover$test[i] <- c(as.numeric(used.rates$turnover >= turnover$hpd.lower[i] & used.rates$turnover <= turnover$hpd.upper[i]))
  turnover$hpd.rel.width[i] <- (turnover$hpd.upper[i]-turnover$hpd.lower[i])/used.rates$turnover
  turnover[i,(ncol(turnover)-21+1):ncol(turnover)] <- lapply(seq(0.,1,0.05), test_vals, turnover$true[i], t$turnover)
  
  samplingProportion$true[i] <- used.rates$sampling_prop
  samplingProportion$estimated[i] <- median(t$samplingProportion)
  samplingProportion$upper[i] <- quantile(t$samplingProportion,0.975)
  samplingProportion$lower[i] <- quantile(t$samplingProportion,0.025)
  
  
  samplingProportion$rel.err.meadian[i] <- abs(median(t$samplingProportion)-used.rates$sampling_prop)/used.rates$sampling_prop
  samplingProportion$rel.err.mean[i] <-abs(mean(t$samplingProportion)-used.rates$sampling_prop)/used.rates$sampling_prop
  samplingProportion$rel.err.mode[i] <- abs(Mode(t$samplingProportion)-used.rates$sampling_prop)/used.rates$sampling_prop
  samplingProportion$cv[i] <- sqrt(exp(sd(log(t$samplingProportion))**2)-1)
  samplingProportion$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingProportion))[1]
  samplingProportion$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingProportion))[2]
  # samplingProportion$test[i] <- c(as.numeric(used.rates$sampling_prop >= samplingProportion$hpd.lower[i] & used.rates$sampling_prop <= samplingProportion$hpd.upper[i]))
  samplingProportion$hpd.rel.width[i] <- (samplingProportion$hpd.upper[i]-samplingProportion$hpd.lower[i])/used.rates$sampling_prop
  samplingProportion[i,(ncol(samplingProportion)-21+1):ncol(samplingProportion)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion$true[i], t$samplingProportion)
  
  samplingAtPresentProb$true[i] <- used.rates$rho
  samplingAtPresentProb$estimated[i] <- median(t$samplingAtPresentProb)
  samplingAtPresentProb$upper[i] <- quantile(t$samplingAtPresentProb,0.975)
  samplingAtPresentProb$lower[i] <- quantile(t$samplingAtPresentProb,0.025)
  
  
  samplingAtPresentProb$rel.err.meadian[i] <- abs(median(t$samplingAtPresentProb)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$rel.err.mean[i] <-abs(mean(t$samplingAtPresentProb)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$rel.err.mode[i] <- abs(Mode(t$samplingAtPresentProb)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$cv[i] <- sqrt(exp(sd(log(t$samplingAtPresentProb))**2)-1)
  samplingAtPresentProb$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingAtPresentProb))[1]
  samplingAtPresentProb$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingAtPresentProb))[2]
  # samplingAtPresentProb$test[i] <- c(as.numeric(used.rates$rho >= samplingAtPresentProb$hpd.lower[i] & used.rates$rho <= samplingAtPresentProb$hpd.upper[i]))
  samplingAtPresentProb$hpd.rel.width[i] <- (samplingAtPresentProb$hpd.upper[i]-samplingAtPresentProb$hpd.lower[i])/used.rates$rho
  samplingAtPresentProb[i,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)] <- lapply(seq(0.,1,0.05), test_vals, samplingAtPresentProb$true[i], t$samplingAtPresentProb)
  
  
  
  origin$true[i] <- used.rates$origin
  origin$estimated[i] <- median(t$origin)
  origin$upper[i] <- quantile(t$origin,0.975)
  origin$lower[i] <- quantile(t$origin,0.025)
  
  
  origin$rel.err.meadian[i] <- abs(median(t$origin)-used.rates$origin)/used.rates$origin
  origin$rel.err.mean[i] <-abs(mean(t$origin)-used.rates$origin)/used.rates$origin
  origin$rel.err.mode[i] <- abs(Mode(t$origin)-used.rates$origin)/used.rates$origin
  origin$cv[i] <- sqrt(exp(sd(log(t$origin))**2)-1)
  origin$hpd.lower[i] <- HPDinterval(as.mcmc(t$origin))[1]
  origin$hpd.upper[i] <- HPDinterval(as.mcmc(t$origin))[2]
  # origin$test[i] <- c(as.numeric(used.rates$origin >= origin$hpd.lower[i] & used.rates$origin <= origin$hpd.upper[i]))
  origin$hpd.rel.width[i] <- (origin$hpd.upper[i]-origin$hpd.lower[i])/used.rates$origin
  origin[i,(ncol(origin)-21+1):ncol(origin)] <- lapply(seq(0.,1,0.05), test_vals, origin$true[i], t$origin)
  
  
  mrca$true[i] <- used.rates$mrca
  mrca$estimated[i] <- median(t$TreeHeight)
  mrca$upper[i] <- quantile(t$TreeHeight,0.975)
  mrca$lower[i] <- quantile(t$TreeHeight,0.025)
  
  
  mrca$rel.err.meadian[i] <- abs(median(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$rel.err.mean[i] <-abs(mean(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$rel.err.mode[i] <- abs(Mode(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$cv[i] <- sqrt(exp(sd(log(t$TreeHeight))**2)-1)
  mrca$hpd.lower[i] <- HPDinterval(as.mcmc(t$TreeHeight))[1]
  mrca$hpd.upper[i] <- HPDinterval(as.mcmc(t$TreeHeight))[2]
  # mrca$test[i] <- c(as.numeric(used.rates$mrca >= mrca$hpd.lower[i] & used.rates$mrca <= mrca$hpd.upper[i]))
  mrca$hpd.rel.width[i] <- (mrca$hpd.upper[i]-mrca$hpd.lower[i])/used.rates$mrca
  mrca[i,(ncol(mrca)-21+1):ncol(mrca)] <- lapply(seq(0.,1,0.05), test_vals, mrca$true[i], t$TreeHeight)
  
}

save.image(file="figures/validationPlotEnvironment.RData")
# load(file="figures/validationPlotEnvironment.RData")

if (!is.null(remove_rows)){
  diversificationRate <- diversificationRate[-remove_rows, ]
  turnover <- turnover[-remove_rows, ]
  mrca <- mrca[-remove_rows, ]
  origin <- origin[-remove_rows, ]
  samplingAtPresentProb <- samplingAtPresentProb[-remove_rows, ]
  samplingProportion <- samplingProportion[-remove_rows, ]
}


#unlist(c((samplingProportion[,(ncol(samplingProportion)-21+1):ncol(samplingProportion)])))
d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate[,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/length(diversificationRate$test_0))

l <- length(diversificationRate$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(diversificationRate[,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate")

# p.qq_div_rate <- p.qq_div_rate+
#   stat_summary(data=dd, aes(x=p, y=vals), fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
#   ggtitle("Diversification rate")

ggsave(plot=p.qq_div_rate,paste("figures/qq_diversificationRate", ".pdf", sep=""),width=6, height=5)


x.min_diversificationRate = min(diversificationRate$true)
x.max_diversificationRate = max(diversificationRate$true)

y.min_diversificationRate = min(diversificationRate$lower)
y.max_diversificationRate = max(diversificationRate$upper)


lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate <- ggplot(diversificationRate)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


ggsave(plot=p.diversificationRate,paste("figures/diversificationRate", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log <- ggplot(diversificationRate)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.diversificationRate_log,paste("figures/diversificationRate_logScale", ".pdf", sep=""),width=6, height=5)

# dt_rel.err.median.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$rel.err.meadian)
# p_rel_err_diversificationRate<- ggplot(dt_rel.err.median.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_rel_err_diversificationRate,paste("figures/diversificationRate_rel_err_median", ".pdf", sep=""),width=6, height=5)
# 
# dt_rel.err.mean.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$rel.err.mean)
# p_rel_err_mean_diversificationRate<- ggplot(dt_rel.err.mean.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_rel_err_mean_diversificationRate,paste("figures/diversificationRate_rel_err_mean", ".pdf", sep=""),width=6, height=5)
# 
# dt_hpd.rel.width.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$hpd.rel.width)
# p_hpd.rel.width.diversificationRate <- ggplot(dt_hpd.rel.width.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_hpd.rel.width.diversificationRate,paste("figures/diversificationRate_rel_hpd_width", ".pdf", sep=""),width=6, height=5)

########## turnover ################


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover[,(ncol(turnover)-21+1):ncol(turnover)], sum)/length(turnover$test_0))

l <- length(turnover$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(turnover[,(ncol(turnover)-21+1):ncol(turnover)]))

p.qq_turnover <- ggplot(d)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_turnover <- p.qq_turnover+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Turnover")

ggsave(plot=p.qq_turnover,paste("figures/qq_turnover", ".pdf", sep=""),width=6, height=5)

x.min_turnover = min(turnover$true)
x.max_turnover = max(turnover$true)

y.min_turnover = min(turnover$lower)
y.max_turnover = max(turnover$upper)


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover <- ggplot(turnover)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.turnover,paste("figures/turnover", ".pdf", sep=""),width=6, height=5)

p.turnover_log <- ggplot(turnover)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.turnover_log,paste("figures/turnover_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('turnover'),
                      "HPD coverage"=c(mean(turnover$test)))
write.csv(hpd.test, file = paste("figures/turnover_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('turnover'),
                       "Average relative HPD width"=c(mean(turnover$hpd.rel.width)))
write.csv(hpd.width, file = paste("figures/turnover_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('turnover'),
                "Average Coefficient of Variation"=c(mean(turnover$cv)))
write.csv(cv, file = paste("figures/turnover_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('turnover'),
                     "Medians"=c(mean(turnover$rel.err.meadian)))
write.csv(medians, file = paste("figures/turnover_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('turnover')", "Means"=c(mean(turnover$rel.err.mean)))
write.csv(means, file = paste("figures/turnover_rel_error_mean",".csv", sep=""))


dt_rel.err.median.turnover = data.frame(param='turnover', value=turnover$rel.err.meadian)
p_rel_err_turnover<- ggplot(dt_rel.err.median.turnover, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.turnover$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_turnover,paste("figures/turnover_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.turnover = data.frame(param='turnover', value=turnover$rel.err.mean)
p_rel_err_mean_turnover<- ggplot(dt_rel.err.mean.turnover, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.turnover$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_turnover,paste("figures/turnover_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.turnover = data.frame(param='turnover', value=turnover$hpd.rel.width)
p_hpd.rel.width.turnover <- ggplot(dt_hpd.rel.width.turnover, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.turnover$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.turnover,paste("figures/turnover_rel_hpd_width", ".pdf", sep=""),width=6, height=5)


########## mrca ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(mrca[,(ncol(mrca)-21+1):ncol(mrca)], sum)/length(mrca$test_0))
p.qq_mrca <- ggplot(d)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") +
  theme_minimal()

x.min_mrca = min(mrca$true)
x.max_mrca = max(mrca$true)

y.min_mrca = min(mrca$lower)
y.max_mrca = max(mrca$upper)


lim.min = min(x.min_mrca, y.min_mrca)
lim.max = max(x.max_mrca, y.max_mrca)

p.mrca <- ggplot(mrca)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


ggsave(plot=p.mrca,paste("figures/mrca", ".pdf", sep=""),width=6, height=5)

p.mrca_log <- ggplot(mrca)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.mrca_log,paste("figures/mrca_logScale", ".pdf", sep=""),width=6, height=5)


########## origin ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(origin[,(ncol(origin)-21+1):ncol(origin)], sum)/length(origin$test_0))

l <- length(origin$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(origin[,(ncol(origin)-21+1):ncol(origin)]))


p.qq_origin <- ggplot(d)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_origin <- p.qq_origin+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Origin")

ggsave(plot=p.qq_origin,paste("figures/qq_origin", ".pdf", sep=""),width=6, height=5)

x.min_origin = min(origin$true)
x.max_origin = max(origin$true)

y.min_origin = min(origin$lower)
y.max_origin = max(origin$upper)


lim.min = min(x.min_origin, y.min_origin)
lim.max = max(x.max_origin, y.max_origin)

p.origin <- ggplot(origin)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


ggsave(plot=p.origin,paste("figures/origin", ".pdf", sep=""),width=6, height=5)

p.origin_log <- ggplot(origin)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.origin_log,paste("figures/origin_logScale", ".pdf", sep=""),width=6, height=5)


########## samplingAtPresentProb ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingAtPresentProb[,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)], sum)/length(samplingAtPresentProb$test_0))


l <- length(samplingAtPresentProb$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(samplingAtPresentProb[,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)]))

p.qq_samplingAtPresentProb <- ggplot(d)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") +
  theme_minimal()

p.qq_samplingAtPresentProb <- p.qq_samplingAtPresentProb+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("samplingAtPresentProb ")

ggsave(plot=p.qq_samplingAtPresentProb,paste("figures/qq_samplingAtPresentProb", ".pdf", sep=""),width=6, height=5)

x.min_samplingAtPresentProb = min(samplingAtPresentProb$true)
x.max_samplingAtPresentProb = max(samplingAtPresentProb$true)

y.min_samplingAtPresentProb = min(samplingAtPresentProb$lower)
y.max_samplingAtPresentProb = max(samplingAtPresentProb$upper)


lim.min = min(x.min_samplingAtPresentProb, y.min_samplingAtPresentProb)
lim.max = max(x.max_samplingAtPresentProb, y.max_samplingAtPresentProb)

p.samplingAtPresentProb <- ggplot(samplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


ggsave(plot=p.samplingAtPresentProb,paste("figures/samplingAtPresentProb", ".pdf", sep=""),width=6, height=5)

p.samplingAtPresentProb_log <- ggplot(samplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingAtPresentProb_log,paste("figures/samplingAtPresentProb_logScale", ".pdf", sep=""),width=6, height=5)


########## samplingProportion ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion[,(ncol(samplingProportion)-21+1):ncol(samplingProportion)], sum)/length(samplingProportion$test_0))

l <- length(samplingProportion$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(samplingProportion[,(ncol(samplingProportion)-21+1):ncol(samplingProportion)]))


p.qq_samplingProportion <- ggplot(d)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") +
  theme_minimal()

p.qq_samplingProportion <- p.qq_samplingProportion+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("samplingProportion")

ggsave(plot=p.qq_samplingProportion,paste("figures/qq_samplingProportion", ".pdf", sep=""),width=6, height=5)

x.min_samplingProportion = min(samplingProportion$true)
x.max_samplingProportion = max(samplingProportion$true)

y.min_samplingProportion = min(samplingProportion$lower)
y.max_samplingProportion = max(samplingProportion$upper)


lim.min = min(x.min_samplingProportion, y.min_samplingProportion)
lim.max = max(x.max_samplingProportion, y.max_samplingProportion)

p.samplingProportion <- ggplot(samplingProportion)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


ggsave(plot=p.samplingProportion,paste("figures/samplingProportion", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_log <- ggplot(samplingProportion)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingProportion_log,paste("figures/samplingProportion_logScale", ".pdf", sep=""),width=6, height=5)



######### tables #############

hpd.test = data.frame("Parameter"=c('diversificationRate',
                                    'turnover',
                                    'samplingProportion',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                      "HPD coverage"=c(mean(diversificationRate$test_95),
                                       mean(turnover$test_95),
                                       mean(samplingProportion$test_95),
                                       mean(samplingAtPresentProb$test_95),
                                       mean(mrca$test_95),
                                       mean(origin$test_95)))
write.csv(hpd.test, file = paste("figures/HPD_test",".csv", sep="" ))


hpd.width = data.frame("Parameter"=c('diversificationRate',
                                     'turnover',
                                     'samplingProportion',
                                     'samplingAtPresentProb',
                                     'mrca',
                                     'origin'),
                       "Average relative HPD width"=c(mean(diversificationRate$hpd.rel.width),
                                                      mean(turnover$hpd.rel.width),
                                                      mean(samplingProportion$hpd.rel.width),
                                                      mean(samplingAtPresentProb$hpd.rel.width),
                                                      mean(mrca$hpd.rel.width),
                                                      mean(origin$hpd.rel.width)))
write.csv(hpd.width, file = paste("figures/HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('diversificationRate',
                              'turnover',
                              'samplingProportion',
                              'samplingAtPresentProb',
                              'mrca',
                              'origin'),
                "Average Coefficient of Variation"=c(mean(diversificationRate$cv),
                                                     mean(turnover$cv),
                                                     mean(samplingProportion$cv),
                                                     mean(samplingAtPresentProb$cv),
                                                     mean(mrca$cv),
                                                     mean(origin$cv)))
write.csv(cv, file = paste("figures/CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('diversificationRate',
                                   'turnover',
                                   'samplingProportion',
                                   'samplingAtPresentProb',
                                   'mrca',
                                   'origin'),
                     "Relative Error at Medians"=c(mean(diversificationRate$rel.err.meadian),
                                 mean(turnover$rel.err.meadian),
                                 mean(samplingProportion$rel.err.meadian),
                                 mean(samplingAtPresentProb$rel.err.meadian),
                                 mean(mrca$rel.err.meadian),
                                 mean(origin$rel.err.meadian)))
write.csv(medians, file = paste("figures/rel_error_median",".csv", sep="" ))

means = data.frame("Parameter"=c('diversificationRate',
                                   'turnover',
                                   'samplingProportion',
                                   'samplingAtPresentProb',
                                   'mrca',
                                   'origin'),
                   "Relative Error at Means"=c(mean(diversificationRate$rel.err.mean),
                                                 mean(turnover$rel.err.mean),
                                                 mean(samplingProportion$rel.err.mean),
                                                 mean(samplingAtPresentProb$rel.err.mean),
                                                 mean(mrca$rel.err.mean),
                                                 mean(origin$rel.err.mean)))
write.csv(means, file = paste("figures/rel_error_mean",".csv", sep=""))

