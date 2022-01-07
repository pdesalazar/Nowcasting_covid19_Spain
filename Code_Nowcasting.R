
#Code for reproduction of analyse and Figure 1 plots in De Salazar et al, intermediate period of analysis
rm(list = ls())

require(MASS)
require(dplyr)
require(vcd)
require(ggplot2)
require(cowplot)
require(rjags) #required for NobBs, see https://cran.r-project.org/web/packages/NobBS/index.html
require(EpiEstim) #see https://cran.r-project.org/web/packages/EpiEstim/index.html

setwd("~/Nowcasting")
source('NobBS_post.R') #source NobBs modified for extraction posterior samples


# Obtain raw data dataset
load("data_Nowcasting.RData")

#Parameter selection 

r<-2 #region of analysis, 1 for Murcia, 2 for Madrid
N_impute<-100 # number imputations for the incomplete time series use 10 for test, 100 for reliable estimates
N_samples <- 10000  #  samples in NobBS, use 1000 for test, 10000 for accurate estimates
#percent_window<-0.95 #percentile of delay distribution for window in NobBs for dynamic window
today<-as.Date("2020-04-09") #select last day of reported cases for window analysis
first_day<-as.Date("2020-03-01") #select first day of reported cases for window analysis
days_analysis<-as.numeric(today-first_day+1) #compute number of days analysed
max_delay_int<-0.95 #percentile for maximum delay upper bound

#parameters for plots

last_day_plot<-as.Date("2020-04-16")
region.m<-matrix(c("MUR","MAD",200,4000), ncol=2)
xlims <- c(as.Date(first_day-1), last_day_plot) # x plot limits
ylims<-as.numeric(c(0,region.m[r,2]))   # y plot limits
size_legend<-22 #size legend
size_text<-16 #size text in legend
b_s<-24 #base size plot


#Plot 1 complete date of sympstoms onset/all(complete +incomplete) case counts

#aggregated cases with complete information by reporting date and site
df_complete<-as.data.frame(table(total[complete.cases(total$onset_day) &
                                         complete.cases(total$report_day),
                                       c("report_day","Area")]))

#aggregated cases (complete and incomplete) by date of report and site
df_all<-as.data.frame(table(total[,c("report_day","Area")]))

df_plot1 <- left_join(df_complete, df_all,by=c("report_day","Area"))%>%
  rename(complete=Freq.x,all=Freq.y,date=report_day)%>%
  mutate(date=as.Date(date))

plot1<-ggplot(subset(df_plot1,Area==region.m[r,1]))+
  geom_col(aes(x=date,y=all),fill="grey80",colour="black",width=0.7)+
  geom_col(aes(x=date,y=complete),fill="dodgerblue3",colour="black",width=0.7,alpha=1)+
  scale_x_date(limits = xlims,date_breaks="2 days",date_labels = "%d %b")+
  ylim(ylims)+
  ylab("No. cases by DOR")+
  xlab("Date of report")+
  theme_classic(base_size=b_s)+
  theme(axis.text.x=element_text(angle=60, hjust=1,size=size_text))

plot1

# Imputation step
#add variable week day 
region_analysis<-total%>%subset(Area==region.m[r,1])%>%
  mutate(wday=weekdays(report_day))%>%
  mutate(wday=ifelse(wday=="Sunday"|wday=="Saturday",1,0))%>%
  subset(report_day>=first_day &report_day<=today,c(onset_day,report_day,wday))%>%
  mutate(day=as.numeric(report_day-min(unique(report_day))))

#subset data with complete information for parametric fitting
total_complete<- subset(region_analysis,(complete.cases(onset_day)&complete.cases(report_day)),
  c(report_day,onset_day,day,wday))%>%
  mutate(delay=as.numeric(report_day-onset_day))%>%
  subset(delay>=0)

#imputation (as loops for clarity)
df_imputed<-data.frame(NULL)
for (i in region_analysis$day%>%unique()%>%sort()){
  day_incomplete<- subset(region_analysis,(!complete.cases(onset_day)&
                                            complete.cases(report_day)&
                                             day==i),
                          c(report_day,onset_day,day,wday))
  
  day_complete<- subset(region_analysis,(complete.cases(onset_day)&
                                         complete.cases(report_day)&
                                         day==i),
                        c(report_day,onset_day,day,wday))
  
  ifelse(nrow(subset(total_complete,day<=i&day>=i-2))>=30,
         subset_total_complete<-subset(total_complete,day<=i&day>=i-2),
         ifelse(nrow(subset(total_complete,day<=i&day>=i-6))>=30,
          subset_total_complete<-subset(total_complete,day<=i&day>=i-6),
          subset_total_complete<-subset(total_complete,day<=i&day>=i-9)))
         
  max_delays<-as.data.frame(table(subset_total_complete$delay))%>%
  mutate(Freq=round(cumsum(as.numeric(table(subset_total_complete$delay)))/
                      sum(as.numeric(table(subset_total_complete$delay))),2))
           
  max_delay<-as.numeric(as.vector(max_delays$Var1)[which.min(abs(max_delays$Freq-max_delay_int))])
  subset_total_complete<-subset(subset_total_complete,
                                delay<=(as.numeric(as.vector(max_delays$Var1)[which.min(abs(max_delays$Freq-max_delay_int))])))
  
  ndata<-subset(day_incomplete,day==i,c(day,wday))[1,]

  ifelse(nrow(subset_total_complete)>=30,
         {nb.out<-glm.nb(delay~wday,data=subset_total_complete);
         pred.means<-predict(nb.out,newdata=ndata,se.fit=TRUE,type="response")},
         out<-median(subset_total_complete$delay))

  estimates<-NULL
  repeat{
    ifelse(nrow(subset_total_complete)>=30,
           delay<-rnbinom(nrow(day_incomplete),mu=(rnorm(1, mean=pred.means$fit,
                                                         sd=pred.means$se.fit)),
                          size=(rnorm(1,mean=nb.out$theta,sd=(nb.out$SE.theta)))),
           delay<-rep(out,nrow(day_incomplete)))
    
    day_incomplete<-day_incomplete%>%mutate(onset_day=report_day-delay)
    imputed<-rbind(day_complete,day_incomplete)
    estimates<-rbind(estimates,imputed)
    if (nrow(estimates)/nrow(imputed) == N_impute){
      break}
  }
  df_imputed<-rbind(df_imputed,estimates%>%mutate(run=rep(c(1:N_impute),each=nrow(imputed))))
}

#Nowcasting
#window of analysis
#N_days<-as.integer(quantile(total_complete$delay, probs=c(percent_window))) #dynamic
N_days<-as.numeric(ceiling(days_analysis*3/4)) #constant (75% of period of analysis)
mult_imp_post <- NULL

for (i in c(1:N_impute)) {
  print(i)
  filled.data <- subset(as.data.frame(df_imputed),run==i&report_day<as.Date(today)) # replace with imputation code
  nowcast <- NobBS.post(data=filled.data, now=as.Date(today), units="1 day",
                        onset_date="onset_day", 
                        report_date="report_day",moving_window=N_days,
                        specs = list(dist="NB",nSamp= N_samples, 
                                     nBurnin=N_samples/10, nAdapt=N_samples/10))
  
  mult_imp_post = cbind(mult_imp_post, nowcast$full.post.n)
}

#Dataframe for ggplot
full_dates=data.frame(onset_day=seq.Date(from=first_day,to=today,by="days"))

#1.Observed time series by date or symptoms

observed<-region_analysis%>%subset((complete.cases(onset_day) &
                            complete.cases(report_day)),onset_day)%>%
                         table()%>%as.data.frame()%>%
                         mutate(onset_day=as.Date(as.character(.)))%>%
                         subset(onset_day>=first_day &onset_day<=today)%>%
                         right_join(full_dates)%>%
                         mutate(Freq=ifelse(is.na(Freq),0,Freq))

observed_cases<-c(observed[order(observed$onset_day),]$Freq,rep (NA,ifelse(today-max(observed$onset_day)>=1,today-max(observed$onset_day),0)))

#2.Imputed+observed time series and intervals by date of symptoms
imputed<-df_imputed%>%subset(report_day<=as.Date(today)&onset_day>=first_day)

imputed_matrix<-matrix(as.vector(t(table(imputed$onset_day,imputed$run))),
                ncol=max(imputed$run),byrow=TRUE)

imputed_cases<- imputed_matrix%>%
  apply( 1, quantile, probs=c(0.0275, 0.5, 0.975))%>%
  t()%>%data.frame()%>%rename(lower_imp= `X2.75.` ,higher_imp=`X97.5.`,center_imp=`X50.`)%>%
  mutate(onset_day=as.Date(sort(unique(imputed$onset_day))))%>%
  right_join(full_dates)

#3.Nowcasted time series (posterior samples from NobBS)
nowcasted_cases<-mult_imp_post%>%apply( 1,quantile, probs=c(0.0275, 0.5, 0.975))%>%
  t()%>%data.frame()%>%rename(lower_now= `X2.75.` ,higher_now=`X97.5.`,center_now=`X50.`)

if(N_days<=as.numeric(days_analysis+1))
  {extra_days<-rep (NA,as.numeric(days_analysis)-N_days)
   extra<-data.frame(lower_now= extra_days ,higher_now=extra_days,center_now=extra_days)
   nowcasted_cases<-rbind(extra,nowcasted_cases)}

#Plot2 Nowcasting/imputed + observed cases by DOS
df_plot2<-cbind(nowcasted_cases,imputed_cases,observed_cases)%>%rename(time=onset_day)

plot2<-ggplot(df_plot2) + 
  geom_line(aes(time,center_now),
            alpha=1,colour="#E69F00",linetype=1,size=0.8)+
  geom_ribbon(aes(time,ymin =lower_now, ymax = higher_now), fill = "#E69F00",alpha=0.2)+
  geom_line(aes(time,center_imp),
            alpha=1,colour="grey60",linetype="solid",size=0.8) +
  geom_col(aes(time,center_imp),fill="grey80",colour="black",width=0.7)+
  geom_ribbon(aes(time,ymin =lower_imp, ymax = higher_imp), fill = "grey60",alpha=0.3)+
  geom_col(aes(time,observed_cases),fill="dodgerblue3",colour="black",width=0.7,alpha=1)+
  scale_x_date(limits = xlims,date_breaks="2 days",date_labels = "%d %b")+
  ylim(c(ylims))+
  xlab("Date of symptoms onset")+
  ylab("No. cases by DOS")+
  theme_classic(base_size=b_s)+
  theme(axis.text.x=element_text(angle=60, hjust=1,size=size_text))
plot2

#estimation time varying R
N_sample_R<-100 #use 10 for test, 100 for more reliable estimates
mean_si<-5
std_si<-1.9
std_mean_si<-1
min_mean_si<-3
max_mean_si<-7
std_std_si<-0.5
min_std_si<-1
max_std_si<-3
#Alternative params.
# mean_si<-7
# std_si<-3.4
# std_mean_si<-2
# min_mean_si<-4
# max_mean_si<-11
# std_std_si<-1
# min_std_si<-1
# max_std_si<-7

#simpler parametric approximation
# config_cori = list(t_start = seq(2,days_analysis ),
#                    t_end = seq(2,days_analysis),
#                    mean_si = mean_si, std_si = std_si)

#uncertainty around the mean and sd estimates ;
config_cori <- make_config(list(t_start = seq(2,days_analysis ),
                                t_end = seq(2,days_analysis),
                                mean_si = mean_si, std_mean_si = std_mean_si,
                                min_mean_si = min_mean_si, max_mean_si = max_mean_si,
                                std_si = std_si, std_std_si = std_std_si,
                                min_std_si = min_std_si, max_std_si = max_std_si))
#simpler parametric approximation
# config_wt = list(t_start = seq(1, days_analysis),
#                  t_end = seq(1, days_analysis),
#                  mean_si = mean_si, std_si = std_si,
#                  n_sim = 100)

# uncertainty around the mean and sd estimates
config_wt = list(t_start = seq(1, days_analysis),
                 t_end = seq(1, days_analysis),
                 mean_si = mean_si, std_mean_si = std_mean_si,
                 min_mean_si = min_mean_si, max_mean_si = max_mean_si,
                 std_si = std_si, std_std_si = std_std_si,
                 min_std_si = min_std_si, max_std_si = max_std_si,n_sim = 100)

#extract posterior from Cori
post_estimates_C<-matrix(NA,ncol=1,nrow=days_analysis)
repeat{
  incidence<-c(imputed_matrix[1:(days_analysis-nrow(mult_imp_post)), 
                       sample(ncol(imputed), 1)],
               mult_imp_post[, sample(ncol(mult_imp_post), 1)])
  R_ts <- estimate_R(incidence,
                     method="uncertain_si",
                     config = config_cori)
  post<-matrix(sample_posterior_R(R_ts,n=length(incidence)*100,
                                  window=1:length(incidence)),ncol=100)
  
  post_estimates_C<-cbind(post_estimates_C,post)
  
  if ((ncol(post_estimates_C)) >= 100*N_sample_R){
    break}
}

estimate_C<-rename(as.data.frame(t(apply(post_estimates_C[,2:(N_sample_R*100)], 
                                         1, quantile, probs=c(0.0275, 0.5, 0.975),na.rm=TRUE)))
                   ,lower_now= `2.75%` ,higher_now=`97.5%`,center_now=`50%`)

#posterior from WT
post_estimates_WT<-matrix(NA,ncol=1,nrow=days_analysis)
repeat{
  incidence<-c(imputed_matrix[0:(days_analysis-nrow(mult_imp_post)), sample(ncol(imputed), 1)],
               mult_imp_post[, sample(ncol(mult_imp_post), 1)])
  R_ts <-  wallinga_teunis( incidence,method="parametric_si",
                            config = config_wt)
  post<-matrix(sample_posterior_R(R_ts,n=length(incidence)*100,
                                  window=1:length(incidence)),ncol=100)
  
  post_estimates_WT<-cbind(post_estimates_WT,post)
  
  if ((ncol(post_estimates_WT)) >= 100*N_sample_R){
    break}
}

estimate_WT<-rename(as.data.frame(t(apply(post_estimates_WT[,2:(N_sample_R*100)], 
                                          1, quantile, probs=c(0.0275, 0.5, 0.975),na.rm=TRUE))),
                    lower_now= `2.75%` ,higher_now=`97.5%`,center_now=`50%`)
#plot3 data frame
df_plot3<-data.frame(mean_wt=estimate_WT$center_now, 
                     higher_wt=estimate_WT$higher_now,
                     lower_wt=estimate_WT$lower_now,
                     mean_c=estimate_C$center_now, 
                     higher_c=estimate_C$higher_now,
                     lower_c=estimate_C$lower_now,
                     time=seq.Date(from=first_day,to=today,by = "days"))

#Plot 3 Rt estimates
plot3<-ggplot(df_plot3)+
  geom_line(aes(time,mean_wt),colour="#D55E00",size=1)+
  geom_ribbon(fill="#D55E00",alpha=0.3,aes(x =time,ymin=lower_wt, ymax=higher_wt ),size=2)+
  geom_line(aes(time,y=mean_c),colour="purple",size=1,alpha=1)+
  geom_ribbon(fill="purple",alpha=0.1,aes(x =time,ymin=lower_c,ymax=higher_c),size=2)+
  geom_hline(yintercept=1,linetype=2,colour="darkred",size=1,alpha=0.6)+
  scale_x_date(limits = xlims,date_breaks="2 days",date_labels = "%d %b")+
  ylab("time-varying R")+
  xlab("Date")+
  theme_classic(base_size=b_s)+
  theme(axis.text.x=element_text(angle=60, hjust=1,size=size_text))

plot3
#Figure 1
Final_plot<-plot_grid(plot1, plot2,plot3, ncol=1, align = "v")
print(Final_plot)


