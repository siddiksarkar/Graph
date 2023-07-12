
library(Hmisc)
library(drc)
list.files()
# upload the .csv or excel file with data#
df1<- read.csv("D:/Siddik Sarkar/Documents/IICB/R results/GraphPadPrism/GraphPadPrism/Doseresponse curve_RHTBH in uM.csv", row.names=1)

par(mfrow= c(1, 2)) # Arrange the graphs with 1 row 4 col#
ID8<- data.frame("dose"= rep(df1$Var.1, 3), "response"= c(df1$ID.8.WT_1, df1$ID.8.WT_2, df1$ID.8.WT_3) )


dr.ID8<- drm(formula=  ID8$response~ID8$dose,
               
               fct = LL.4(fixed= c(NA,0,100,NA),names = c("Hill slope", "Min", "Max", "EC50"))
)


plot(dr.ID8,
     type= "all",
     col="black", #auto color selection for curves
     pch=16,
     xlim= c(0.01,100),
     lwd=1,
     xlab= expression(paste("Re(htbh) ",mu,"M")),
     ylab= "Cell viability (%)",
     main= "PA-1 treated with Re(htbh)",
)
minor.tick(nx=2, ny=2)
axis(1, at= c(0.01,0.1,1,5,10,50,100), labels=c(0.01,0.1,1,5,10,50,100))
abline(h= 50, v=dr.ID8$coefficients[2],  lty=3 )
mtext(text= paste("EC50:", round(dr.ID8$coefficients[2],digits=2)),side=3, at= 25, line= -3, cex= 0.75)

plot(dr.ID8,
     type= "confidence",
     col="black", #auto color selection for curves
     add=T
)
 # 2nd plot
#PA-1

ID8p53KO<- data.frame("dose"= rep(df1$Var.1, 3), "response"= c(df1$ID.8.P53.knockout_1, df1$ID.8.P53.knockout_2, df1$ID.8.P53.knockout_3) )

dr.ID8p53KO<-drm(formula=  ID8p53KO$response~ID8p53KO$dose,
            fct = LL.4(fixed= c(NA,0,100,NA),names = c("Hill slope", "Min", "Max", "EC50"))
            
)

plot(dr.ID8p53KO,
     type= "all",
     col="blue", #auto color selection for curves
     pch=16,
     lwd=1,
     xlab= expression(paste("Re(htbh) ", mu,"M")),
     ylab= "Cell viability (%)",
     main= "ID8p53KO treated with Re(htbh)",
     ylim= c(0,100)
)
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
axis(1, at= c(0.01,0.1,1,5,10,50), labels=c(0.01,0.1,1,5,10,50) )
abline(h= 50,v=dr.ID8p53KO$coefficients[2], lty=3 )
mtext(text= paste("EC50:", round(dr.ID8p53KO$coefficients[2], 2)),side=3, at= 25, line= -3, cex=0.75)
#Add confidence line#
plot(dr.ID8p53KO,
     type= "confidence",
     col="blue", #auto color selection for curves
     pch=16,
     lwd=1,
     add=T)
######################################
# Two sigmoidal plot in a single graph#
######################################
par(mfrow= c(1,1))
plot(dr.ID8, col= "purple", pch=15,xlab= "Conc (uM)", ylab= "Viability (%)",
     xlim= c(0.01, 100),
     ylim= c(0, 125),
     main= "Dose response curve of Re(htbh) in cancer cells",
     cex.main= 0.8)
plot(dr.ID8p53KO,col="red", pch=16, lty=1, add=T)
abline(h=50, lty=2)
minor.tick(nx=2, ny=2)
axis(1, at= c(0.01,0.1,1,5,10,50, 100), labels=c(0.01,0.1,1,5,10,50, 100) )
legend("bottomleft", legend= c("ID8", "ID8p53KO"), col= c("purple", "red"),
       lty=c(1,1), pch= c(15,16), cex= 0.75)

par(mfrow=c(2,2))
EC50<- data.frame("Cell lines"= c("ID8", "ID8P53KO"), "EC50 in uM"=c(dr.ID8$coefficients["EC50:(Intercept)"],
                                                                     dr.ID8$coefficients["EC50:(Intercept)"]))
EC50$EC50.in.uM <- round(EC50$EC50.in.uM, digits=2)

