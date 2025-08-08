
rm(list=ls())

setwd("~/GitHub/PSD_decomposition_ARMA")
source("~/GitHub/PSD_decomposition_ARMA/Functions/arma_kalman.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/kalman_opti_sin_me.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/kalman_opti_con_me_unknown.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/kalman_opti_con_me_known.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/Residuos.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/Criterios.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/calculate_AR_coeff.R")
source("~/GitHub/PSD_decomposition_ARMA/Functions/PSD_ARMA.R")

AGN_data<-read.table("xtej1550.dat",header=FALSE)
dim(AGN_data)
head(AGN_data)
summary(AGN_data)

AGN_data<-cbind(AGN_data,log(AGN_data[,2]))
colnames(AGN_data) <- c('TIME','RATE','LOG_RATE')

mean=mean(log(AGN_data$RATE))
LOG_rates=log(AGN_data$RATE)-mean
AGN_data=cbind(AGN_data,LOG_rates)

dim(AGN_data)
#425981      4
head(AGN_data)
dt=diff(AGN_data$TIME)[1]


#PRUEBA: SEGMENTO 72
T=5000
s=271
ind = ((s-1)*T/4+1):((s-1)*T/4+T+1)
times_dt = AGN_data$TIME[ind]
times=times_dt-times_dt[1]
times=times[-c(1)]/dt


rates=AGN_data$LOG_RATE[ind]-mean(AGN_data$LOG_RATE[ind])
rates=rates[-c(1)]

#SERIE DE TIEMPO FINAL
plot(times,rates,t='l')


#Promedio de periodogramas
PER_int=list()
Num_div=T/1000
lng=floor(T/Num_div)
lng

freqs=(1:(T/2))/T
N_f=length(freqs)
freqs=freqs[seq(1,N_f,by=1)]
N_f=length(freqs)
freqs=matrix(freqs,ncol=1,nrow=N_f)
N_f

PER=rep(0,N_f)
for(j in 1:(Num_div)){
  
  PER_int[[j]]=rep(0,N_f)
  
  for(i in 1:length(freqs)){
    
    ind=((j-1)*lng+1):(j*lng)
    
    A = sum(rates[ind]*cos(2*pi*times[ind]*freqs[i]))
    B = sum(rates[ind]*sin(2*pi*times[ind]*freqs[i]))
    PER_int[[j]][i]=(A^2+B^2)/length(rates[ind])
    
  }
  PER=cbind(PER,PER_int[[j]])
}

PER_RAW=PER[,-c(1)]
PER_RAW_mean=apply(PER_RAW,1,mean)
plot(log10(freqs),log10(PER_RAW_mean),t='l',col='gray',ylim=c(-3,1),
     ylab='PSD',xlab='Frecuency',main='Periodogram mean')

#Varianza del error estimada, encontrada con PSD que ajusta los dos peak (ARMA5,4 no pasa los supuestos)
noise_level0=0.01075112
par(mfrow=c(1,2))
#Gráfico periodograma sin escala
plot(freqs/dt,PER_RAW_mean,t='l',xlab='Frequency',ylab='PSD')
abline(h=noise_level0,col=2)
#Gráfico periodograma escala log
plot(log10(freqs/dt),log10(PER_RAW_mean),t='l',xlab='log10 of Frequency',ylab='log10 of PSD')
abline(h=log10(noise_level0),col=2)

##########################################################################################
#########################################################################################

#Ajuste ARMA5,4 sin error de medición 

p=5
q=3
r1=0.9814427+0.1445698i
r2=0.9814427-0.1445698i
r3=0.9367932+0.2912987i
r4=0.9367932-0.2912987i
r5=-0.9488249
phi=calculate_AR_coeff(c(r1,r2,r3,r4,r5))
order=c(p,q)

#Calculo del punto inicial usando la función arima de R y fijando phi
fit_arma=arima(rates,order = c(p,0,q),include.mean = FALSE,method = "ML",
               optim.control=list(maxit = 20000),fixed=c(phi,rep(NA,q)))
x0=c(fit_arma$coef,sqrt(fit_arma$sigma2))

#Ajuste del modelo ARMA5,4 con varianza de errores de medición desconocida (estimada)
aux=optim(x0, kalman_opti_sin_me,control = list(maxit = 30000),
          y=rbind(rates),order=order)
#Convergencia
aux$convergence

aux$par#Causalidad
min(Mod(polyroot(c(1,-aux$par[1:p]))))>1

#Coeficientes estimados
ar_coef=aux$par[1:p]
ma_coef=aux$par[(p+1):(p+q)]
sigma2z=aux$par[sum(order)+1]^2

#Cálculo AIC,BIC,AICC
criterios_ajuste=Criterios(x0=x0,order=order,
                           y=rbind(rates),sigmaeta=numeric())
criterios_ajuste$AIC
criterios_ajuste$BIC

#Ajuste en el dominio de la frecuencia
PSD_hat=PSD_ARMA(phi_coefs = ar_coef, theta_coefs = ma_coef, sigmaz = sqrt(sigma2z), 
                 freq = 2*pi*freqs,components=FALSE) 

plot(log10(freqs),log10(PER_RAW_mean),t='l',col='gray',ylim=c(-3,1),
     ylab='PSD',xlab='Frecuency',main=paste('Fit PSD ARMA(',p,',',q,') with Vt unknown'))
lines(log10(freqs),log10(PSD_hat$PSD*2*pi),lwd=1,col=2)
legend("bottomleft", legend=c("Periodogram","PSD ARMA(5,3)"),
       col=c("gray",2), lty=1, cex=0.8)

#Residuos
resid=Residuos(x0=x0,order=order,
               y=rbind(rates),sigmaeta=numeric())

plot(resid,t='l',xlab='Time',ylab='Residuos')
acf(resid)
pacf(resid)
hist(resid)

##########################################################################################
#########################################################################################

#Ajuste ARMA5,4 con error de medición conocido (varianza)

p=5
q=4
r1=0.9814427+0.1445698i
r2=0.9814427-0.1445698i
r3=0.9367932+0.2912987i
r4=0.9367932-0.2912987i
r5=-0.9488249
phi=calculate_AR_coeff(c(r1,r2,r3,r4,r5))
order=c(p,q)
sigma2eta=noise_level0

#Calculo del punto inicial usando la función arima de R y fijando phi
fit_arma=arima(rates,order = c(p,0,q),include.mean = FALSE,method = "ML",
               optim.control=list(maxit = 20000),fixed=c(phi,rep(NA,q)))
x0=c(fit_arma$coef,sqrt(fit_arma$sigma2))

#Ajuste del modelo ARMA5,4 con varianza de errores de medición desconocida (estimada)
aux=optim(x0, kalman_opti_con_me_known,control = list(maxit = 30000),
          y=rbind(rates),sigmaeta=sqrt(sigma2eta),order=order)
#Convergencia
aux$convergence

#Causalidad
min(Mod(polyroot(c(1,-aux$par[1:p]))))>1

#Coeficientes estimados
ar_coef=aux$par[1:p]
ma_coef=aux$par[(p+1):(p+q)]
sigma2z=aux$par[sum(order)+1]^2

#Cálculo AIC,BIC,AICC
criterios_ajuste=Criterios(x0=x0,order=order,
                           y=rbind(rates),sigmaeta=sqrt(sigma2eta))
criterios_ajuste$AIC
criterios_ajuste$BIC

#Ajuste en el dominio de la frecuencia
PSD_hat=PSD_ARMA(phi_coefs = ar_coef, theta_coefs = ma_coef, sigmaz = sqrt(sigma2z), 
                 freq = 2*pi*freqs,components=FALSE) 

plot(log10(freqs),log10(PER_RAW_mean),t='l',col='gray',ylim=c(-3,1),
     ylab='PSD',xlab='Frecuency',main=paste('Fit PSD ARMA(',p,',',q,') with Vt unknown'))
lines(log10(freqs),log10(PSD_hat$PSD*2*pi+sigma2eta),lwd=3,col=2)
lines(log10(freqs),log10(PSD_hat$PSD*2*pi),lwd=1,col=3)
abline(h=log10(sigma2eta),lwd=1,col=4)
legend("bottomleft", legend=c("Periodogram", "PSD ARMA(5,4)+ME","PSD ARMA(5,4)","ME"),
       col=c("gray",2,3,4), lty=1, cex=0.8)

#Residuos
resid=Residuos(x0=x0,order=order,
               y=rbind(rates),sigmaeta=sqrt(sigma2eta))

plot(resid,t='l',xlab='Time',ylab='Residuos')
acf(resid)
pacf(resid)
hist(resid)


##########################################################################################
#########################################################################################

#Ajuste ARMA5,4 con error de medición desconocido (estimado)
#Este modelo captura los dos peaks pero no se validan los supuestos del modelo

### vt unknown two peaks
sigmaeta=sqrt(0.03163)
p=5
q=4
r1=0.9814427+0.1445698i
r2=0.9814427-0.1445698i
r3=0.9367932+0.2912987i
r4=0.9367932-0.2912987i
r5=-0.9488249
phi=calculate_AR_coeff(c(r1,r2,r3,r4,r5))
order=c(p,q)

#Calculo del punto inicial usando la función arima de R y fijando phi
fit_arma=arima(rates,order = c(p,0,q),include.mean = FALSE,method = "ML",
               optim.control=list(maxit = 20000),fixed=c(phi,rep(NA,q)))
x0=c(fit_arma$coef,sqrt(fit_arma$sigma2),sigmaeta)

#Ajuste del modelo ARMA5,4 con varianza de errores de medición desconocida (estimada)
aux=optim(x0, kalman_opti_con_me_unknown,control = list(maxit = 30000),
          y=rbind(rates),order=order)
#Convergencia
aux$convergence

#Causalidad
min(Mod(polyroot(c(1,-aux$par[1:p]))))>1

#Coeficientes estimados
ar_coef=aux$par[1:p]
ma_coef=aux$par[(p+1):(p+q)]
sigma2z=aux$par[p+q+1]^2
sigma2eta=aux$par[p+q+2]^2

#Cálculo AIC,BIC,AICC
criterios_ajuste=Criterios(x0=aux$par[-length(aux$par)],order=c(p,q),
                y=rbind(rates),sigmaeta=aux$par[sum(order)+2])
criterios_ajuste$AIC
criterios_ajuste$BIC

#Ajuste en el dominio de la frecuencia
PSD_hat=PSD_ARMA(phi_coefs = ar_coef, theta_coefs = ma_coef, sigmaz = sqrt(sigma2z), 
                 freq = 2*pi*freqs,components=FALSE) 

plot(log10(freqs),log10(PER_RAW_mean),t='l',col='gray',ylim=c(-3,1),
     ylab='PSD',xlab='Frecuency',main=paste('Fit PSD ARMA(',p,',',q,') with Vt unknown'))
lines(log10(freqs),log10(PSD_hat$PSD*2*pi+sigma2eta),lwd=3,col=2)
lines(log10(freqs),log10(PSD_hat$PSD*2*pi),lwd=1,col=3)
abline(h=log10(sigma2eta),lwd=1,col=4)
legend("bottomleft", legend=c("Periodogram", "PSD ARMA(5,4)+ME","PSD ARMA(5,4)","ME"),
       col=c("gray",2,3,4), lty=1, cex=0.8)

#Residuos
resid=Residuos(x0=aux$par[-length(aux$par)],order=order,
               y=rbind(rates),sigmaeta=aux$par[sum(order)+2])

plot(resid,t='l',xlab='Time',ylab='Residuos')
acf(resid)
pacf(resid)
hist(resid)



