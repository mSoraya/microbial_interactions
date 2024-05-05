library(circlize)

col.pal=c(Bifidobacterium_asteroides="#E3319D",Gilliamella_apicola="green",
          Lactobacillus_apis="blue",Lactobacillus_kullabergensis="red",
          Lactobacillus_mellifer_Bin="#FFD801",Lactobacillus_mellis_Hon="#7FFFD4",
          Snodgrassella_alvi="#FFA62F")



mat1<- matrix(c(0,	0.17,0.29,	0.25,	0.32,	0.31,	0.23,	
                0,	0,	0.2,	0.22,	0.28,	0.25,	0.13,
                0,	0,	0,	0.26,	0.32,	0.33,	0.25,
                0,	0,	0,	0	,0.31	,0.32,	.031,
                0,	0,	0,	0,	0	,0.32,	0.4,
                0,	0,	0,	0,	0,	0	,0.37,
                0,	0,	0,	0,	0,	0,	0),nrow=7,ncol = 7)

sample_Matrix<-matrix(mat1,nrow = 7,ncol=7,
                      dimnames = list(c("Bifidobacterium_asteroides","Gilliamella_apicola",
                                        "Lactobacillus_apis","Lactobacillus_kullabergensis",
                                        "Lactobacillus_mellifer_Bin","Lactobacillus_mellis_Hon",
                                        "Snodgrassella_alvi"),
                                      c("Bifidobacterium_asteroides","Gilliamella_apicola",
                                        "Lactobacillus_apis","Lactobacillus_kullabergensis",
                                        "Lactobacillus_mellifer_Bin","Lactobacillus_mellis_Hon",
                                        "Snodgrassella_alvi")))
#col_fun = function(x) ifelse(x == 32, "#001639", ifelse(x==31,"#003382",ifelse(x==30,"#003d9b",ifelse(x==29,"#0043ab",ifelse(x==28,"#0050cc",ifelse(x==25,"#0668ff",ifelse(x==24,"#0060f4",ifelse(x==23,"#1672ff",ifelse(x==21,"#4f94ff",ifelse(x==19,"#70a8ff",ifelse(x==17,"#90bcff",ifelse(x==16,"#b9d5fe",ifelse(x==14,"#d2e3fe","#eaf2ff")))))))))))))

col_fun = function(x) ifelse(x == 0.4, "#001639", ifelse(x==0.37,"#003382",ifelse(x==0.33,"#003d9b",ifelse(x==0.32,"#0043ab",ifelse(x==0.31,"#0050cc",ifelse(x==0.29,"#0668ff",ifelse(x==0.28,"#0060f4",ifelse(x==0.26,"#1672ff",ifelse(x==0.25,"#4f94ff",ifelse(x==0.23,"#70a8ff",ifelse(x==0.22,"#90bcff",ifelse(x==0.2,"#b9d5fe",ifelse(x==0.17,"#d2e3fe","#eaf2ff")))))))))))))

#col_fun = function(x) ifelse(x == 32, "#0056dc", ifelse(x==31,"#0060f4",ifelse(x==30,"#0e6dff",ifelse(x==29,"#2e81ff",ifelse(x==28,"#478fff",ifelse(x==25,"#0668ff",ifelse(x==24,"#90bcff",ifelse(x==23,"#1672ff",ifelse(x==21,"#4f94ff",ifelse(x==19,"#70a8ff",ifelse(x==17,"#90bcff",ifelse(x==16,"#b9d5fe",ifelse(x==14,"#d2e3fe","#eaf2ff")))))))))))))
par(cex = 1.77,mar = c(0, 0, 0, 0))
circos.par("track.height" = 0.01)
chordDiagram(sample_Matrix,grid.col = col.pal,annotationTrack = c("name","grid"), scale = 0.1,link.lwd = 0.1,col = col_fun, transparency = 0.4)
circos.clear()

mat2<- matrix(c(0,	102,	84,	90,	94,	94,	63,
               0,	0,	110,	117,	103,	100,	98,
               0,	0,	0,	111,	91,	90,	64,
               0,	0,	0,	0	,97	,92,	63,
               0,	0,	0,	0,	0	,102,	61,
               0,	0,	0,	0,	0,	0	,23,
               0,	0,	0,	0,	0,	0,	0),nrow=7,ncol = 7)


sample_Matrix<-matrix(mat2,nrow = 7,ncol=7,
                      dimnames = list(c("Bifidobacterium_asteroides","Gilliamella_apicola",
                                        "Lactobacillus_apis","Lactobacillus_kullabergensis",
                                        "Lactobacillus_mellifer_Bin","Lactobacillus_mellis_Hon",
                                        "Snodgrassella_alvi"),
                                      c("Bifidobacterium_asteroides","Gilliamella_apicola",
                                        "Lactobacillus_apis","Lactobacillus_kullabergensis",
                                        "Lactobacillus_mellifer_Bin","Lactobacillus_mellis_Hon",
                                        "Snodgrassella_alvi")))
col_fun = function(x) ifelse(x == 117, "#001639", ifelse(x==111,"#003382",ifelse(x==110,"#003d9b",ifelse(x==103,"#0043ab",ifelse(x==102,"#0050cc",ifelse(x==100,"#0668ff",ifelse(x==98,"#0060f4",ifelse(x==97,"#1672ff",ifelse(x==94,"#4f94ff",ifelse(x==90,"#70a8ff",ifelse(x==84,"#90bcff",ifelse(x==64,"#b9d5fe",ifelse(x==63,"#d2e3fe",ifelse(x==92,"#5f9eff",ifelse(x==91,"#68a3ff","#eaf2ff")))))))))))))))
par(cex = 1.77,mar = c(0, 0, 0, 0))
circos.par("track.height" = 0.01)
chordDiagram(sample_Matrix,grid.col = col.pal,annotationTrack = c("name","grid"), scale = 0.1,link.lwd = 0.1,col = col_fun, transparency = 0.4)

dev.off()


