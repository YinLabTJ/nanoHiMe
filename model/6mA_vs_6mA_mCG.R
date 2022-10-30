library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
dat <- read.table("6mA_vs_6mA_mCG.input",header = TRUE)
p<-ggplot(data=dat)+geom_hline(aes(yintercept=0),size=0.2)+geom_vline(aes(xintercept=0),size=0.2)+geom_point(mapping=aes(x=mA,y=mA_mCG,colour=Class,shape=Class,size=Class))+scale_shape_manual(values=c(16, 16, 16))+labs(title=args[1])+scale_color_manual(values=c('#F8766D','#00BA38','#619CFF'))+scale_size_manual(values=c(1, 1, 1))+geom_abline(slope=1,intercept=0, colour="#696969", linetype="dashed")+scale_x_continuous(limit=c(55,125),breaks=c(60,75,90,105,120))+scale_y_continuous(limit=c(55,125),breaks=c(60,75,90,105,120))+ coord_fixed(ratio=1) + theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line.x=element_line(linetype=1,color="black",size=0.2),
  axis.line.y=element_line(linetype=1,color="black",size=0.2),
#  axis.ticks.x=element_line(size=0),
#  axis.ticks.y=element_line(size=0),
  plot.title = element_text(hjust = 0.5)
  #axis.text.y = element_text(hjust = 40),
  #axis.text.x = element_text(vjust = 5)
  )
#ver <- data.frame(x=c(55,55,55,55,55,55,55,125),y=c(55,55,55,55,55,55,55,125))
#p + geom_segment(data=ver,aes(x=c(55,55,55,55,120,105,90,75),y=c(120,105,90,75,55,55,55,55),xend=c(50,50,50,50,120,105,90,75),yend=c(125,105,90,75,50,50,50,50)),size=0.2)
#p + labs(title = args[1])
p + ggsave(file="color.png", width=6, height=6)

