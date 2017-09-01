require(dplyr)
require(tidyr)
require(IgorR)
require(ggplot2)
require(data.table)

folder <-"//nanobud/QC/2017/Experiments/0 - Active Experiments/AFM Short Loop Devices/Short Loop Wafer 2/Stack 2/Raw Data/"
Savedir<-"//nanobud/QC/2017/Experiments/0 - Active Experiments/AFM Short Loop Devices/Short Loop Wafer 2/Stack 2/"
FN<-"115 nm Bit"

Read_AFM_File <- function(file){
  filepath <- file
  #filepath<-"//nanobud/QC/2017/Experiments/0 - Active Experiments/AFM Short Loop Devices/Short Loop Wafer 2/Stack 4/Raw Data/slw2_s4_3v0029.ibw"
  sourceData <- read.ibw(filepath,Verbose = FALSE, ReturnTimeSeries = FALSE,
                         MakeWave = FALSE, HeaderOnly = FALSE) 
  
  data <- as.data.frame(sourceData)[1:6]
  colnames(data)<-c("Zsensor_m", "Defl_m", "Curr_A", "Curr2_A", "Bias_V", "Raw")
  #data$File <- filepath
  
  filename<-tools::file_path_sans_ext(basename(filepath))
  filename2<-as.character(filepath)
  filename2<-strsplit(filename2, split="_")
  filename2<-unlist(filename2)
  filename2<-filename2[3]
  filename3<-as.character(filename2)
  filename3<-strsplit(filename2, split=".", fixed=TRUE)
  filename3<-unlist(filename3)
  filename3<-filename3[1]
  filename4<-as.character(filename3)
  filename4<-strsplit(filename4, split="v", fixed=TRUE)
  filename4<-unlist(filename4)
  filename4<-filename4[2]
  data$File<-filename4
  
  
  data3<-as.data.frame(unlist(attributes(sourceData)))
  note<-as.data.frame(strsplit(as.character(data3[4,1]), "\r")) 
  colnames(note) <- c("Values")
  Infot <- note %>% separate(Values, c("Parameter", "Values"), sep = ":") %>% 
    t() %>% as.data.frame(stringsAsFactors=FALSE)
  colnames(Infot)<-Infot[1,]
  Infot<-Infot[-1,]
  
  ##Find Indexes to parse data file by approach, dwell, and retract portion of curves
  Indexes<-unlist(strsplit(as.character(Infot$Indexes), ","))
  StartIndex<-as.numeric(Indexes[1])
  EndApproachIndex<-as.numeric(Indexes[2])
  EndDwellIndex<-as.numeric(Indexes[3])
  EndRetractIndex<-as.numeric(Indexes[4])
  
  
  ##Split data into approach, dwell, and retract portion of curves
  ApproachData<-data[1:EndApproachIndex,]
  ApproachData$Portion<-"Approach"
  DwellData<-data[(EndApproachIndex+1):EndDwellIndex,]
  DwellData$Portion<-"Dwell"
  RetractData<-data[(EndDwellIndex+1):EndRetractIndex,]
  RetractData$Portion<-"Retract"
  
  Dataframe1<-rbind(ApproachData, DwellData, RetractData)
  
  
  
  SpringConstant<-as.numeric(Infot$SpringConstant)
  DeflOffset<-mean(Dataframe1$Defl_m[1:50])
  NoBiasCurr<-DwellData
  NoBiasCurr$Bias_V<-round(NoBiasCurr$Bias_V, 3)
  NoBiasCurr<-subset(NoBiasCurr, NoBiasCurr$Bias_V==0)
  CurrOffset<-mean(ApproachData$Curr_A[1:50])
  Curr2Offset<-mean(ApproachData$Curr2_A[1:50])
  
  Dataframe1$Defl_m<-as.numeric(Dataframe1$Defl_m)
  Dataframe1$Defl_Corr_m<-Dataframe1$Defl_m-DeflOffset
  Dataframe1$Curr_Corr_A<-Dataframe1$Curr_A-CurrOffset
  Dataframe1$Curr_Corr_nA<-Dataframe1$Curr_Corr_A*1e9
  Dataframe1$Curr2_Corr_A<-Dataframe1$Curr2_A-Curr2Offset
  Dataframe1$Curr2_Corr_uA<-Dataframe1$Curr2_Corr_A*1e6
  Dataframe1$Defl_Corr_nm<- Dataframe1$Defl_Corr_m*1e9
  Dataframe1$Resistance_Ohm<-Dataframe1$Bias_V/Dataframe1$Curr_Corr_A
  Dataframe1$Resistance_kOhm<-Dataframe1$Resistance_Ohm/1000
  Dataframe1$Force_N<-Dataframe1$Defl_Corr_m*SpringConstant
  Dataframe1$Force_nN<-Dataframe1$Force_N*1e9
  Interval<-1/as.numeric(Infot$NumPtsPerSec)
  Dataframe1$Time_s<- seq(from=0, to=((nrow(Dataframe1)-1)*Interval), by=Interval)
  
  
  Dataframe2<-Dataframe1[Dataframe1$Portion =="Dwell",]
  minDwelltime<-min(Dataframe2$Time_s)
  Dataframe2$Corr_Bias_V<-Dataframe2$Bias_V*-1
  Dataframe2$Curr2_Corr_uA<-Dataframe2$Curr2_Corr_uA*-1
  Dataframe2$Corr_Time_s<-Dataframe2$Time_s-minDwelltime
  Dataframe2$dydx<-c(diff(Dataframe2$Corr_Bias_V)/diff(Dataframe2$Corr_Time_s),0)
  
  SweepV<-subset(Dataframe2, Dataframe2$dydx!=0)
  SweepV$SetType <- ifelse(SweepV$Corr_Bias_V>0, "Reset", "Set")
  
  ReadV<-subset(Dataframe2, Dataframe2$dydx==0)
  ReadV<-ReadV[2:nrow(ReadV),]
  ReadV<-subset(ReadV, ReadV$Corr_Bias_V<1)
  ReadV<-subset(ReadV, ReadV$Corr_Bias_V>-1)
  ReadV$bin<-cut(ReadV$Corr_Time_s, seq(from=0, to=52.08, by=2.48),labels=c(0,1,1.01,2,2.01,3,3.01,4,4.01,5,5.01,6,6.01,7,7.01,8,8.01,9,9.01,10,10.01))
  ReadV<-na.omit(ReadV)
  
  
  BinnedData2<-plyr::ddply(ReadV, c("bin"), tail,500) 
  
  BinnedData2 %>% group_by(bin, File) %>%
    summarise(meanCurrent_uA=mean(Curr2_Corr_uA), meanBias_V=mean(Corr_Bias_V)) %>%
    ungroup() %>% mutate(n=seq(1:21),
                         Type=c("Virgin",rep(c("Reset", "Set"),10)),
                         Resistance_MOhm = meanBias_V/meanCurrent_uA
    ) %>%
    select(-bin) -> BinnedReadData1
  
  BinnedReadData1$n<-c(0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
  
  ResetSweep<-subset(SweepV, SweepV$SetType=="Reset")
  ResetSweep$Dir<- ifelse(ResetSweep$dydx<0, "Reverse", "Forward")
  
  SetSweep<-subset(SweepV, SweepV$SetType=="Set")
  SetSweep$Dir<- ifelse(SetSweep$dydx>0,"Reverse", "Forward")
  
  ResetSweep$bin<-cut(ResetSweep$Corr_Time_s, seq(from=0, to=48, by=4.8), labels=c(1,2,3,4,5,6,7,8,9,10))
  ResetSweep<-na.omit(ResetSweep)
  
  ResetVoltage1 <- plyr::ddply(ResetSweep, .variables = c("bin", "File", "Dir"), subset,
                               subset= Curr2_Corr_uA==max(Curr2_Corr_uA),
                               select=c(Corr_Bias_V,Curr2_Corr_uA ))
  
  SetSweep$bin<-cut(SetSweep$Corr_Time_s, seq(from=0, to=48, by=4.8), labels=c(1,2,3,4,5,6,7,8,9,10))
  SetSweep<-na.omit(SetSweep)
  SweepData1<-rbind(SetSweep, ResetSweep)  
  
  minV<-round(min(SweepData1$Corr_Bias_V), 1)
  maxV<-round(max(SweepData1$Corr_Bias_V), 1)
  SweepData1$bin2<-cut(SweepData1$Corr_Bias_V, seq(from=minV, to=maxV, by=0.1))
  SweepData1<-na.omit(SweepData1)

  
  SweepData1 %>%
    group_by(bin, bin2, File, Dir, SetType)%>% 
    summarise(Curr2_Corr_uA=mean(Curr2_Corr_uA), Corr_Bias_V=mean(Corr_Bias_V), Corr_Time_s=mean(Corr_Time_s), dydx=mean(dydx)) ->
    SweepData2

  ResetSweep2<-subset(SweepData2, SweepData2$SetType=="Reset")
  ResetSweep2$bin<-cut(ResetSweep2$Corr_Time_s, seq(from=0, to=48, by=4.8), labels=c(1,2,3,4,5,6,7,8,9,10))
  ResetSweep2<-na.omit(ResetSweep)
  
  ResetVoltage2 <- plyr::ddply(ResetSweep2, .variables = c("bin", "File"), subset,
                               subset= Curr2_Corr_uA==max(Curr2_Corr_uA),
                               select=c(Corr_Bias_V,Curr2_Corr_uA ))

  ResetVoltage3 <- plyr::ddply(ResetVoltage2, .variables = c("bin", "File"), subset,
                               subset= Corr_Bias_V==max(Corr_Bias_V),
                               select=c(Corr_Bias_V,Curr2_Corr_uA ))
  
  
  SetSweep2<-subset(SweepData2, SweepData2$SetType=="Set")
  SetSweep2$bin<-cut(SetSweep2$Corr_Time_s, seq(from=0, to=48, by=4.8), labels=c(1,2,3,4,5,6,7,8,9,10))
  SetSweep2<-na.omit(SetSweep2)
  
  SetVoltage1 <-plyr::ddply(SetSweep2, .variables = c("bin", "File"), subset,
                            subset= Corr_Bias_V==min(Corr_Bias_V),
                            select=c(Corr_Bias_V,Curr2_Corr_uA ))
  

  
  
  return(list(Dataframe=Dataframe1, SweepData=SweepData2, 
              BinnedReadData=BinnedReadData1, ResetVoltage=ResetVoltage3, SetVoltage=SetVoltage1))
}


Read_AFM_Folder <- function(folder){
  file_names <- list.files(path = folder , full.names = T )
  
  ReadAll <- Map(c, lapply(file_names, Read_AFM_File))
  
  SweepData <- data.table::rbindlist(sapply(1:length(ReadAll), FUN = function(X){
    return(ReadAll[[X]][2])
    }))
  print("2")
  
  BinnedReadData <- data.table::rbindlist(sapply(1:length(ReadAll), FUN = function(X){
    return(ReadAll[[X]][3])
      }))
  print("3")
  
  ResetVoltage <- data.table::rbindlist(sapply(1:length(ReadAll), FUN = function(X){
    return(ReadAll[[X]][4])
     }))
  print("4")
  gc()
  
  SetVoltage <- data.table::rbindlist(sapply(1:length(ReadAll), FUN = function(X){
    return(ReadAll[[X]][5])
  }))
  print("5")
  rm(ReadAll)
  gc()
  
  
  Dataframe <- data.table::rbindlist(lapply(file_names, FUN=function(x) return(Read_AFM_File(x)$Dataframe)))
  print("done 1")
  
  return(list(Dataframe=Dataframe, SweepData=SweepData, 
              BinnedReadData=BinnedReadData, ResetVoltage=ResetVoltage, SetVoltage=SetVoltage))
  
}

system.time(OUTPUT <- Read_AFM_Folder(folder))

Dataframe<-as.data.frame(OUTPUT[1])
SweepData<-as.data.frame(OUTPUT[2])
BinnedReadData<-as.data.frame(OUTPUT[3])
ResetVoltage<-as.data.frame(OUTPUT[4])
SetVoltage<-as.data.frame(OUTPUT[5])


##DataAdjust
{


colnames(Dataframe)<-c("Zsensor_m", "Defl_m", "Curr_A", "Curr2_A", "Bias_V", "Raw", "File", "Portion", "Defl_Corr_m", "Curr_Corr_A","Curr_Corr_nA", "Curr2_Corr_A", "Curr2_Corr_uA", "Defl_Corr_nm", "Resistance_Ohm", "Resistance_kOhm", "Force_N", "Force_nN", "Time_s")
colnames(BinnedReadData)<-c("File","meanCurrent_uA", "meanBias_V", "n", "Type", "Resistance_MOhm")
colnames(ResetVoltage)<-c("bin", "File","Reset_Voltage_V", "Reset_Current_uA")
colnames(SetVoltage)<-c("bin", "File", "Corr_Bias_V", "Curr2_Corr_uA")
colnames(SweepData)<-c("bin", "bin2", "File", "Dir", "SetType", "Curr2_Corr_uA", "Corr_Bias_V", "Corr_Time_s", "dydx")

ResetSweepData<-subset(SweepData, SweepData$SetType=="Reset")
ResetSweepData$Resistance_MOhm<-ResetSweepData$Corr_Bias_V/ResetSweepData$Curr2_Corr_uA
SetSweepData<-subset(SweepData, SweepData$SetType=="Set")
SetSweepData$Curr2_Corr_uA<-SetSweepData$Curr2_Corr_uA*-1
SetSweepData$Corr_Bias_V<-SetSweepData$Corr_Bias_V*-1

AvgI<-mean(ResetVoltage$Reset_Current_uA)
cutoffI<-AvgI/20
ResetVoltage$Dead<-ifelse(ResetVoltage$Reset_Current_uA>cutoffI & ResetVoltage$Reset_Voltage_V>3.94, "Stuck Set", ifelse(ResetVoltage$Reset_Current_uA<cutoffI, "Dead", "Alive"))

Deadstats<-as.data.frame.matrix(xtabs(~ResetVoltage$Dead+ResetVoltage$bin))

AliveSweepData<-merge(ResetVoltage, SweepData, by=c("File", "bin"), all=TRUE)
AliveSweepData<-subset(AliveSweepData, AliveSweepData$Dead=="Alive")
AliveSweepData<-subset(AliveSweepData, AliveSweepData$Dir=="Forward")
FullSweepData<-AliveSweepData
FullSweepData<-FullSweepData[,c(1,2,8,9,10)]
colnames(FullSweepData)<-c("Bit", "Cycle", "Sweep_Type", "Current_uA", "Voltage_V")
FullSweepData$Resistance_MOhm<-abs(FullSweepData$Voltage_V/FullSweepData$Current_uA)
AliveSweepDataReset<-subset(AliveSweepData, AliveSweepData$SetType=="Reset")
AliveSweepDataSet<-subset(AliveSweepData, AliveSweepData$SetType=="Set")


AliveSweepDataReset %>%
  group_by(bin, bin2)%>% 
  summarise(mean_Curr2_Corr_uA=mean(Curr2_Corr_uA), mean_Corr_Bias_V=mean(Corr_Bias_V), Curr_SD=sd(Curr2_Corr_uA, na.rm=TRUE), Bias_SD=sd(Corr_Bias_V, na.rm=TRUE)) ->
  BinnedSweepData

BinnedSweepData$Curr_UL<-BinnedSweepData$mean_Curr2_Corr_uA+BinnedSweepData$Curr_SD
BinnedSweepData$Curr_LL<-BinnedSweepData$mean_Curr2_Corr_uA-BinnedSweepData$Curr_SD
BinnedSweepData$Curr_LL<-ifelse(BinnedSweepData$Curr_LL<0, 0.1, BinnedSweepData$Curr_LL)
}

##Save Data
{
  sweep2savefile<-paste0(FN, " Set and Reset Sweeps.csv" )
  sweep2csvfilesave<-paste0(Savedir,sweep2savefile )
  
  fwrite(FullSweepData, file = sweep2csvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
  statssavefile<-paste0(FN, " Dead Stats.csv" )
  statscsvfilesave<-paste0(Savedir,statssavefile )
  
  fwrite(Deadstats, file = statscsvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
  
  sweepsavefile<-paste0(FN, " Sweep Data.csv" )
  sweepcsvfilesave<-paste0(Savedir,sweepsavefile )
  
  fwrite(SweepData, file = sweepcsvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
  
  readsavefile<-paste0(FN, " Read Data.csv" )
  readcsvfilesave<-paste0(Savedir,readsavefile )
  
  fwrite(BinnedReadData, file = readcsvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
  
  resetsavefile<-paste0(FN, " Reset Data.csv" )
  resetcsvfilesave<-paste0(Savedir,resetsavefile )
  
  fwrite(ResetVoltage, file = resetcsvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
  
  setsavefile<-paste0(FN, " Set Data.csv" )
  setcsvfilesave<-paste0(Savedir,setsavefile )
  
  fwrite(SetVoltage, file = setcsvfilesave, append = FALSE, quote = "auto",
         sep = ",", eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double", logicalAsInt = FALSE, dateTimeAs = "write.csv",
         showProgress = getOption("datatable.showProgress"),
         verbose = getOption("datatable.verbose"))
  
}

##Plot Data
{
  binnedinplot<-qplot(data = BinnedReadData, x=n, y=meanCurrent_uA, col=Type, ylab = "Read Current / uA")+facet_wrap(~File, nrow = 10)+scale_x_continuous(breaks=seq(from=0, to=10, by=1))
  binnedinplot   
  binnedinplotsave<-paste0(Savedir,FN, "-read current.png")
  ggsave(binnedinplotsave, binnedinplot, width=16, height=8, dpi = 300)
  
  
  binnedinplot<-qplot(data = FullSweepData, x=Voltage_V, y=log10(Current_uA), col=Cycle, xlim = c(0,4), ylim=c(-3,3))+geom_line()+facet_wrap(~Bit, nrow = 2)
  binnedinplot   
  
  
  DCIVplot<-qplot(data = SweepData, x=Corr_Bias_V, y=Curr2_Corr_uA, col=bin, xlab = "Voltage / V", ylab = "Current / uA")+facet_wrap(~File, nrow = 4)
  DCIVplot  
  DCIVplotsave<-paste0(Savedir,FN, "-DCIV.png")
  ggsave(DCIVplotsave, DCIVplot, width=16, height=8, dpi = 300)
  
  logDCIVplot<-qplot(data = ResetSweepData, x=Corr_Bias_V, y=log10(Curr2_Corr_uA), col=bin, ylim=c(-2, 3), xlab = "Voltage / V", ylab = "log(Current) / uA")+facet_wrap(~File, nrow = 10)
  logDCIVplot  
  logDCIVplotsave<-paste0(Savedir,FN, "-log DCIV Reset.png")
  ggsave(logDCIVplotsave, logDCIVplot, width=16, height=8, dpi = 300)
  
  logDCIVplot<-ggplot(AliveSweepData, aes(x=Corr_Bias_V, y=log10(Curr2_Corr_uA), color=bin))+xlab("Voltage / V") + ylab ("log(Current) / uA")+ ylim(-2, 3)+ xlim(0,4)+ geom_line(size=0.5)+facet_wrap(~bin, nrow = 1)
  logDCIVplot  

  logDCIVplot<-ggplot(BinnedSweepData, aes(x=mean_Corr_Bias_V, y=log10(mean_Curr2_Corr_uA), color=bin))+xlab("Voltage / V") + ylab ("log(Current) / uA")+ ylim(-2, 3)+ xlim(0,4)+ geom_line(size=0.5)+facet_wrap(~bin, nrow = 1)
  logDCIVplot  
  
  pd<-position_dodge(.5)
  logerrorplot<-ggplot(BinnedSweepData, aes(x=mean_Corr_Bias_V, y=log10(mean_Curr2_Corr_uA), color=bin))+geom_errorbar(aes(ymin=log10(Curr_LL), ymax=log10(Curr_UL)), color="gray", width=0.1, position=pd)+geom_line(position=pd)+geom_point(position=pd)
  logerrorplot<-logerrorplot+xlab("Voltage / V") + ylab ("log(Current) / uA")+xlim(0,4)+ylim(-1,2.5)+geom_line(size=0.5)+facet_wrap(~bin, nrow = 1)
  logerrorplot
  logerrorplotsave<-paste0(Savedir,FN, "-log DCIV Reset with error.png")
  ggsave(logerrorplotsave, logerrorplot, width=8, height=4, dpi = 300)
  
  errorplot<-ggplot(BinnedSweepData, aes(x=mean_Corr_Bias_V, y=mean_Curr2_Corr_uA, color=bin))+geom_errorbar(aes(ymin=Curr_LL, ymax=Curr_UL), color="gray", width=0.1, position=pd)+geom_line(position=pd)+geom_point(position=pd)
  errorplot<-errorplot+xlab("Voltage / V") + ylab("Current / uA")+xlim(0,4)+ geom_line(size=0.5)+facet_wrap(~bin, nrow = 1)
  errorplot
  errorplotsave<-paste0(Savedir,FN, "-DCIV Reset with error.png")
  ggsave(errorplotsave, errorplot, width=8, height=4, dpi = 300)
  
  
  # logDCRVplot<-qplot(data = ResetSweepData[ResetSweepData$Dir=="Forward",], x=Corr_Bias_V, y=(Resistance_MOhm), col=bin, ylim=c(0, 10))+facet_wrap(~File, nrow = 10)
  # logDCRVplot  
  # logDCRVplotsave<-paste0(Savedir,FN, "-log DCRV Reset.png")
  # ggsave(logDCRVplotsave, logDCRVplot, width=16, height=8, dpi = 300)
  
  
  logDCIVplot2<-qplot(data = SetSweepData, x=Corr_Bias_V, y=log10(Curr2_Corr_uA), col=bin, ylim=c(-2, 3), xlab = "Voltage / V", ylab = "log(Current) / uA")+facet_wrap(~File, nrow = 4)
  #logDCIVplot2  
  logDCIVplot2save<-paste0(Savedir,FN, "-log DCIV Set.png")
  ggsave(logDCIVplot2save, logDCIVplot2, width=16, height=8, dpi = 300)
 
  AliveData<-subset(ResetVoltage,ResetVoltage$Dead=="Alive" )
   
  ResetVplot<-qplot(data = AliveData, x=bin, y=Reset_Voltage_V, xlab = "n", ylab = "Reset Voltage / V")+facet_wrap(~File, nrow = 4)
  ResetVplot  
  ResetVplotsave<-paste0(Savedir,FN, "-Reset Voltage.png")
  ggsave(ResetVplotsave, ResetVplot, width=16, height=8, dpi = 300)
  
  ResetIplot<-qplot(data = AliveData, x=bin, y=Reset_Current_uA, xlab = "n", ylab = "Reset Current / uA")+facet_wrap(~File, nrow = 4)
  #ResetIplot  
  ResetIplotsave<-paste0(Savedir,FN, "-Reset Current.png")
  ggsave(ResetIplotsave, ResetIplot, width=16, height=8, dpi = 300)
  
  SetIplot<-qplot(data = SetVoltage, x=bin, y=-1*Curr2_Corr_uA, xlab = "n", ylab = "Set Current / uA")+facet_wrap(~File, nrow = 4)
  #SetIplot  
  SetIplotsave<-paste0(Savedir,FN, "-Set Current.png")
  ggsave(SetIplotsave, SetIplot, width=16, height=8, dpi = 300)



Resetsummary <- plyr::ddply(AliveData, .variables=c("bin"), summarise, medV = round(median(Reset_Voltage_V),2), medI=round(median(Reset_Current_uA),2))

ResetVboxplot<-ggplot(AliveData, aes(x=bin, y=Reset_Voltage_V))+ xlab("n")+ylab("Reset Voltage / uA")+ geom_boxplot(notch = TRUE)
ResetVboxplot<-ResetVboxplot+geom_text(data = Resetsummary, aes(x = bin, y = medV, label = medV),size = 3, vjust = -1.5)
ResetVboxplot
ResetVhistplotsave<-paste0(Savedir,FN, "-average reset V.png")
ggsave(ResetVhistplotsave, ResetVboxplot, width=8, height=3, dpi = 300)


ResetIboxplot<-ggplot(AliveData, aes(x=bin, y=Reset_Current_uA))+xlab("n")+ylab("Reset Current / uA")+ geom_boxplot(notch = TRUE)
ResetIboxplot<-ResetIboxplot+geom_text(data = Resetsummary, aes(x = bin, y = medI, label = medI),size = 3, vjust = -1.5)
ResetIboxplot
ResetIhistplotsave<-paste0(Savedir,FN, "-average reset I.png")
ggsave(ResetIhistplotsave, ResetIboxplot, width=8, height=3, dpi = 300)


logDCRVplot<-qplot(data = AliveSweepDataReset, x=Corr_Bias_V, y=log10(Corr_Bias_V/Curr2_Corr_uA), col=bin, ylim=c(-2, 2), xlab = "Voltage / V", ylab = "log(Resistance / MOhm)")+geom_line()+facet_wrap(~File, nrow = 4)
logDCRVplot
logDCRVplotsave<-paste0(Savedir,FN, "-log DCRV.png")
ggsave(logDCRVplotsave, logDCRVplot, width=16, height=8, dpi = 300)

logDCRVSetplot<-qplot(data = AliveSweepDataSet, x=-Corr_Bias_V, y=log10(Corr_Bias_V/Curr2_Corr_uA), col=bin, ylim=c(-2, 4), xlab = "Voltage / V", ylab = "log(Resistance / MOhm)")+geom_line()+facet_wrap(~File, nrow = 4)
logDCRVSetplot
logDCRVSetplotsave<-paste0(Savedir,FN, "-log DCRV Set.png")
ggsave(logDCRVSetplotsave, logDCRVSetplot, width=16, height=8, dpi = 300)


}



Reset1plot<-qplot(data = FullSweepData[FullSweepData$Cycle==1,], x=Voltage_V, y=(Current_uA), col=Bit, ylim=c(-0.5,300), xlim=c(0,4), xlab = "Voltage / V", ylab = "Resistance / MOhm")+geom_line()+facet_wrap(~Bit, nrow=5)
Reset1plot  


logDCIVplot<-qplot(data = ResetSweepData, x=Corr_Bias_V, y=(Curr2_Corr_uA), col=bin, ylim=c(0, 1000), xlab = "Voltage / V", ylab = "log(Current) / uA")+facet_wrap(~File, nrow = 10)
logDCIVplot  

logDCIVplot<-qplot(data = ResetSweepData[ResetSweepData$Dir=="Forward",], x=Corr_Bias_V, y=log10(Resistance_MOhm), col=bin, ylim=c(-2, 3), xlab = "Voltage / V", ylab = "log(Current) / uA")+facet_wrap(~File, nrow = 10)
logDCIVplot

logDCIVplot<-qplot(data = AliveSweepData, x=Corr_Bias_V, y=log10(Corr_Bias_V/Curr2_Corr_uA), col=bin, ylim=c(-2, 3), xlab = "Voltage / V", ylab = "log(Current) / uA")+geom_line()+facet_wrap(~File, nrow = 4)
logDCIVplot




logDCIVplot<-qplot(data = AliveSweepData[AliveSweepData$bin==1,], x=Corr_Bias_V, y=log10(Corr_Bias_V/Curr2_Corr_uA), col=File, ylim=c(-2, 3), xlab = "Voltage / V", ylab = "log(Resistance / MOhm)")+geom_line(size=1)
logDCIVplot

 
# ##Single File Plot
{
# BinnedReadData$n2<-(BinnedReadData$n)
# for (i in (1:nrow(BinnedReadData)))
# {
#   if (BinnedReadData$Type[i]=="Reset") {BinnedReadData$n2[i]<-round(((BinnedReadData$n[i]+1)/2-0.5),0)}
#   if (BinnedReadData$Type[i]=="Set") {BinnedReadData$n2[i]<-round(((BinnedReadData$n[i])/2-0.5),0)} 
#   if (BinnedReadData$Type[i]=="Virgin") {BinnedReadData$n2[i]<-0}
# }
# 
# ResetSweepData$bin<-as.character(ResetSweepData$bin)
# SetSweepData$bin<-as.character(SetSweepData$bin)
# 
# onefile<-"Ru_80nm_cycling_0015.ibw"
# 
# readcurrentplot<-qplot(data = BinnedReadData[BinnedReadData$File==onefile,], x=n2, y=meanCurrent_uA, col=Type, xlab = "n", ylab="log(Current / uA)", main="Read Current")+scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10"))+geom_point(size=3)
# resetplot<-qplot(data = ResetSweepData[ResetSweepData$File==onefile,], x=Corr_Bias_V, y=log10(Curr2_Corr_uA), col=bin, ylim=c(-2, 2.5), ylab="Current / uA",  main="Reset DCIV")+labs(color="n")+ aes(group=rev(bin))
# setplot<-qplot(data = SetSweepData[SetSweepData$File==onefile,], x=Corr_Bias_V, y=log10(Curr2_Corr_uA), col=bin, ylim=c(-2, 2.5), ylab="Current / uA",  main="Set DCIV")+labs(color="n")+ aes(group=rev(bin))
# voltplot<-qplot(data = ResetVoltage[ResetVoltage$File==onefile,], x=bin, y=Corr_Bias_V, xlab = "n", ylab = "Reset Voltage / V", ylim=c(0,4))+geom_point(size=3)
# currentplot<-qplot(data = ResetVoltage[ResetVoltage$File==onefile,], x=bin, y=Curr2_Corr_uA, xlab = "n", ylab = "Reset Current / uA", ylim=c(0,100))+geom_point(size=3)
# 
# combineplot<-grid.arrange(resetplot,setplot,voltplot,currentplot,readcurrentplot, ncol=2)
# combineplot
# combineplotsave<-paste0(Savedir,FN, "-file 15.png")
# ggsave(combineplotsave, combineplot, width=16, height=8, dpi = 300)
# 
}