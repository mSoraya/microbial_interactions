
#Supplementary information of "Metabolic resource overlap impacts competition among phyllosphere bacteria"
#available at https://doi.org/10.1038/s41396-023-01459-0.
source("Metabolic resource overlap impacts competition among phyllosphere bacteria\\paper_cusper_competition-main\\scripts\\cusper_libraries.R")



finaldata=read.csv("D:\\PhD_thesis\\Articles\\MRO\\Metabolic resource overlap impacts competition among phyllosphere bacteria\\paper_cusper_competition-main\\data\\finaldata.csv", fill = TRUE)

colors = c(rep("#FFF68F",7),rep("#FFD700",7),rep("#FFC125",7),rep("#DAA520",7))

names =c("ArthL145","MethL85","mono","PkP19E3","PssB728a","RhodL225","SmFR1",
         "ArthL145","MethL85","mono","PkP19E3","PssB728a","RhodL225","SmFR1",
         "ArthL145","MethL85","mono","PkP19E3","PssB728a","RhodL225","SmFR1",
         "ArthL145","MethL85","mono","PkP19E3","PssB728a","RhodL225","SmFR1")

boxplot(finaldata$cfu ~ finaldata$trt*finaldata$time,col=colors,names=names,las=2,main="Population density of Pe299R over time (0, 24, 36, 48h)",
                                    ylab = " ",xlab = " ",cex.axis=0.7)
title(ylab=expression(paste("Population density (CFU gF"," W"^-1,")")), line=2, cex.lab=1.2)
legend("bottomleft", legend = c("* 0.1< P < 0.5 ","**  0.5 ≤ P ≤ 1","time 0","time 24","time 36","time 48") ,
       col = c(" "," ","#FFF68F" , "#FFD700","#FFC125","#DAA520"), bty = "n", pch=c(NA,NA, 20,20,20,20) , pt.cex =3, cex =1, horiz = FALSE, inset = c(-0.07,0.7),y.intersp=0.5,fill=, x.intersp=0.25,border = "black")
   

