##%%%%%%%%%%%%%%%%%%##

### Frost Analysis ###

##%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%%%%%%%#

### PRELIMINARY ANALYSIS ###

#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#------------------------#
##### Load Libraries #####
#------------------------#

library(data.table)
library(dplyr)
library(terra)
library(ggplot2)
library(scales)
library(colorspace)
library(emmeans)
library(corrplot)
library(lubridate)

#-------------------#
##### Load Data #####
#-------------------#

data<-fread("Data/Results_filt.csv")
points<-vect("Data/beech_40x40.shp")
source(file="Scripts/wvioplot.r")
df_filtered<-fread("Data/df_filtered.csv")
quad_model_centered <- readRDS("Data/quad_model_centered.rds")
data_with_ndvi <- fread('Data/H2_data_allyear.csv')
weekly_NDVI<-fread("Data/all_data_week_class.csv")
data_filt<-fread("Data/df_filtered.csv")

#/////////////////////#
##### MAIN FIGURE #####
#/////////////////////#

# Define elevation band thresholds
TH<-c(0, 800, 1000, 1200, 1400, 2500)

# Initialize a vector to store the dates
first_dates <- vector("list", length(TH) - 1)

# Loop through elevation bands
for (th in 1:(length(TH)-1)){
  
  # Define elevation min and max
  hMIN<-TH[th]
  hMAX<-TH[th+1]
  
  # Subset points within current elevation band
  pts<-points[points$DTM>hMIN & points$DTM<=hMAX,]
  
  # List NDVI quantile CSV files for this elevation band
  FLgam<-list.files("Data/Quantiles/",".csv$",full.names=T)
  FLgam<-FLgam[grep("NDVI_",FLgam)]
  FLgam<-FLgam[grep(hMIN,FLgam)]
  FLgam<-FLgam[grep(hMAX,FLgam)]
  
  # Read Wilcoxon test results for NDVI and temperature
  wt<-as.data.frame(fread("Data/wilcoxon_test_NDVI.csv"))
  wt_T<-as.data.frame(fread("Data/Temp/wilcoxon_test_temp.csv"))
  
  ### --- Process NDVI Wilcoxon results --- ###
  # Sort by ID and keep only points in current band
  wt<-wt[order(wt$id),]
  wt<-wt[wt$id %in% pts$id,]
  
  pts<-pts[pts$id %in% wt$id,]
  
  # Extract date from column names
  date_wt<-as.Date(gsub("R_","",names(wt)[-c(1:2)]),format="%Y%m%d")
  
  # Keep only April–August 2019
  wt<-wt[,c(1,2,which(date_wt>=as.Date("20190401",format="%Y%m%d") & date_wt<=as.Date("20190828",format="%Y%m%d"))+2)]
  
  # Round p-values
  wt[,-c(1:2)]<-round(wt[,-c(1:2)],3)
  
  # Split into "L" (lower tail) and "G" (greater tail) subsets
  wt_l<-wt[wt$Dir=="L",-2]
  wt_g<-wt[wt$Dir=="G",-2]
  
  # Convert p-values into binary (sig=0, non-sig=1), treating NA as non-sig
  
  wt_l_bin<-wt_l
  wt_l_bin[,2:ncol(wt_l_bin)][wt_l_bin[,2:ncol(wt_l_bin)]<=0.05]<-0
  wt_l_bin[,2:ncol(wt_l_bin)][wt_l_bin[,2:ncol(wt_l_bin)]>0.05]<-1
  wt_l_bin[,2:ncol(wt_l_bin)][is.na(wt_l_bin[,2:ncol(wt_l_bin)])]<-1
  
  wt_g_bin<-wt_g
  wt_g_bin[,2:ncol(wt_g_bin)][wt_g_bin[,2:ncol(wt_g_bin)]<=0.05]<-0
  wt_g_bin[,2:ncol(wt_g_bin)][wt_g_bin[,2:ncol(wt_g_bin)]>0.05]<-1
  wt_g_bin[,2:ncol(wt_g_bin)][is.na(wt_g_bin[,2:ncol(wt_g_bin)])]<-1
  
  # Calculate percentage of significant results per date
  wt_l_perc<-round(100*((dim(wt_l)[1]-colSums(wt_l_bin[,-1],na.rm=T))/dim(wt_l)[1]),0)            
  wt_g_perc<-round(100*((dim(wt_g)[1]-colSums(wt_g_bin[,-1],na.rm=T))/dim(wt_g)[1]),0)            
  
  # Get corresponding dates
  date_wt_l<-as.Date(gsub("R_","",names(wt_l_perc)),format="%Y%m%d")
  date_wt_g<-as.Date(gsub("R_","",names(wt_g_perc)),format="%Y%m%d")
  
  ### --- Process Temperature Wilcoxon results --- ###
  
  # Sort and filter to current points
  wt_T<-wt_T[order(wt_T$id),]
  wt_T<-wt_T[wt_T$id %in% pts$id,]
  
  # Extract dates and filter April–August 2019
  date_wt_T<-as.Date(gsub("R_","",names(wt_T)[-c(1:2)]),format="%Y%m%d")
  wt_T<-wt_T[,c(1,2,which(date_wt_T>=as.Date("20190401",format="%Y%m%d") & date_wt_T<=as.Date("20190828",format="%Y%m%d"))+2)]
  wt_T[,-c(1:2)]<-round(wt_T[,-c(1:2)],3)
  
  # Split by direction (L vs G)
  wt_T_l<-wt_T[wt_T$Dir=="L",-2]
  wt_T_g<-wt_T[wt_T$Dir=="G",-2]
  
  # Binarise p-values
  wt_T_l_bin<-wt_T_l
  wt_T_l_bin[,2:ncol(wt_T_l_bin)][wt_T_l_bin[,2:ncol(wt_T_l_bin)]<=0.05]<-0
  wt_T_l_bin[,2:ncol(wt_T_l_bin)][wt_T_l_bin[,2:ncol(wt_T_l_bin)]>0.05]<-1
  wt_T_l_bin[,2:ncol(wt_T_l_bin)][is.na(wt_T_l_bin[,2:ncol(wt_T_l_bin)])]<-1
  
  wt_T_g_bin<-wt_T_g
  wt_T_g_bin[,2:ncol(wt_T_g_bin)][wt_T_g_bin[,2:ncol(wt_T_g_bin)]<=0.05]<-0
  wt_T_g_bin[,2:ncol(wt_T_g_bin)][wt_T_g_bin[,2:ncol(wt_T_g_bin)]>0.05]<-1
  wt_T_g_bin[,2:ncol(wt_T_g_bin)][is.na(wt_T_g_bin[,2:ncol(wt_T_g_bin)])]<-1
  
  # Percentages of significant results
  wt_T_l_perc<-round(100*((dim(wt_T_l)[1]-colSums(wt_T_l_bin[,-1],na.rm=T))/dim(wt_T_l)[1]),0)            
  wt_T_g_perc<-round(100*((dim(wt_T_g)[1]-colSums(wt_T_g_bin[,-1],na.rm=T))/dim(wt_T_g)[1]),0)            
  
  # Dates for plotting
  date_wt_T_l<-as.Date(gsub("R_","",names(wt_T_l_perc)),format="%Y%m%d")
  date_wt_T_g<-as.Date(gsub("R_","",names(wt_T_g_perc)),format="%Y%m%d")
  
  ### --- Load NDVI reference and 2019 quantile data --- ###
  
  meREF<-read.table(FLgam[grep(paste("NDVI_median_ref_",hMIN,"_",hMAX,sep=""),FLgam)],sep=" ",header=T)
  q25REF<-read.table(FLgam[grep(paste("NDVI_q10_reference_",hMIN,"_",hMAX,sep=""),FLgam)],header=T)
  q75REF<-read.table(FLgam[grep(paste("NDVI_q90_reference_",hMIN,"_",hMAX,sep=""),FLgam)],sep=" ",header=T)
  
  me<-read.table(FLgam[grep(paste("NDVI_median_2019_",hMIN,"_",hMAX,sep=""),FLgam)],sep=" ",header=T)
  q25<-read.table(FLgam[grep(paste("NDVI_q10_2019_",hMIN,"_",hMAX,sep=""),FLgam)],sep=" ",header=T)
  q75<-read.table(FLgam[grep(paste("NDVI_q90_2019_",hMIN,"_",hMAX,sep=""),FLgam)],sep=" ",header=T)
  
  # Convert Date columns
  me$Date<-as.Date(me$Date,format="%Y-%m-%d")
  q25$Date<-as.Date(q25$Date,format="%Y-%m-%d")
  q75$Date<-as.Date(q75$Date,format="%Y-%m-%d")
  
  meREF$Date<-as.Date(meREF$Date,format="%Y-%m-%d")
  q25REF$Date<-as.Date(q25REF$Date,format="%Y-%m-%d")
  q75REF$Date<-as.Date(q75REF$Date,format="%Y-%m-%d")
  
  ### --- TRIM LINES TO SHADED POLYGONS RANGE --- ###
  shade_start <- min(c(date_wt_T_l, date_wt_T_g, date_wt_l, date_wt_g))
  shade_end   <- max(c(date_wt_T_l, date_wt_T_g, date_wt_l, date_wt_g)) + 7  # add 7 days for last week
  
  # --- Trim NDVI data to shading range ---
  me      <- me[me$Date >= shade_start & me$Date <= shade_end, ]
  meREF   <- meREF[meREF$Date >= shade_start & meREF$Date <= shade_end, ]
  q25     <- q25[q25$Date >= shade_start & q25$Date <= shade_end, ]
  q75     <- q75[q75$Date >= shade_start & q75$Date <= shade_end, ]
  q25REF  <- q25REF[q25REF$Date >= shade_start & q25REF$Date <= shade_end, ]
  q75REF  <- q75REF[q75REF$Date >= shade_start & q75REF$Date <= shade_end, ]
  
  
  ### --- Categorise percentages into bins (1–10) for plotting --- ###
  # (applied separately for NDVI G, NDVI L, Temp G, Temp L)
  
  wt_g_perc_co<-wt_g_perc
  wt_g_perc_co[wt_g_perc<=10]<-1
  wt_g_perc_co[wt_g_perc>10 & wt_g_perc<=20]<-2
  wt_g_perc_co[wt_g_perc>20 & wt_g_perc<=30]<-3
  wt_g_perc_co[wt_g_perc>30 & wt_g_perc<=40]<-4
  wt_g_perc_co[wt_g_perc>40 & wt_g_perc<=50]<-5
  wt_g_perc_co[wt_g_perc>50 & wt_g_perc<=60]<-6
  wt_g_perc_co[wt_g_perc>60 & wt_g_perc<=70]<-7
  wt_g_perc_co[wt_g_perc>70 & wt_g_perc<=80]<-8
  wt_g_perc_co[wt_g_perc>80 & wt_g_perc<=90]<-9
  wt_g_perc_co[wt_g_perc>90]<-10
  
  wt_l_perc_co<-wt_l_perc
  wt_l_perc_co[wt_l_perc<=10]<-1
  wt_l_perc_co[wt_l_perc>10 & wt_l_perc<=20]<-2
  wt_l_perc_co[wt_l_perc>20 & wt_l_perc<=30]<-3
  wt_l_perc_co[wt_l_perc>30 & wt_l_perc<=40]<-4
  wt_l_perc_co[wt_l_perc>40 & wt_l_perc<=50]<-5
  wt_l_perc_co[wt_l_perc>50 & wt_l_perc<=60]<-6
  wt_l_perc_co[wt_l_perc>60 & wt_l_perc<=70]<-7
  wt_l_perc_co[wt_l_perc>70 & wt_l_perc<=80]<-8
  wt_l_perc_co[wt_l_perc>80 & wt_l_perc<=90]<-9
  wt_l_perc_co[wt_l_perc>90]<-10
  
  wt_T_g_perc_co<-wt_T_g_perc
  wt_T_g_perc_co[wt_T_g_perc<=10]<-1
  wt_T_g_perc_co[wt_T_g_perc>10 & wt_T_g_perc<=20]<-2
  wt_T_g_perc_co[wt_T_g_perc>20 & wt_T_g_perc<=30]<-3
  wt_T_g_perc_co[wt_T_g_perc>30 & wt_T_g_perc<=40]<-4
  wt_T_g_perc_co[wt_T_g_perc>40 & wt_T_g_perc<=50]<-5
  wt_T_g_perc_co[wt_T_g_perc>50 & wt_T_g_perc<=60]<-6
  wt_T_g_perc_co[wt_T_g_perc>60 & wt_T_g_perc<=70]<-7
  wt_T_g_perc_co[wt_T_g_perc>70 & wt_T_g_perc<=80]<-8
  wt_T_g_perc_co[wt_T_g_perc>80 & wt_T_g_perc<=90]<-9
  wt_T_g_perc_co[wt_T_g_perc>90]<-10
  
  wt_T_l_perc_co<-wt_T_l_perc
  wt_T_l_perc_co[wt_T_l_perc<=10]<-1
  wt_T_l_perc_co[wt_T_l_perc>10 & wt_T_l_perc<=20]<-2
  wt_T_l_perc_co[wt_T_l_perc>20 & wt_T_l_perc<=30]<-3
  wt_T_l_perc_co[wt_T_l_perc>30 & wt_T_l_perc<=40]<-4
  wt_T_l_perc_co[wt_T_l_perc>40 & wt_T_l_perc<=50]<-5
  wt_T_l_perc_co[wt_T_l_perc>50 & wt_T_l_perc<=60]<-6
  wt_T_l_perc_co[wt_T_l_perc>60 & wt_T_l_perc<=70]<-7
  wt_T_l_perc_co[wt_T_l_perc>70 & wt_T_l_perc<=80]<-8
  wt_T_l_perc_co[wt_T_l_perc>80 & wt_T_l_perc<=90]<-9
  wt_T_l_perc_co[wt_T_l_perc>90]<-10
  
  # Define colour palette for shading
  
  pal <- rev(                    # ensure light → dark
    adjustcolor(                # apply transparency
      sequential_hcl(12, palette = "Blues 3"),  # keep last 10 (drop lightest 3)
      alpha.f = 0.8
    )
  )
  
  pal<- pal[3:12]
  
  # Save output as PDF
  
  cairo_pdf( paste0("Figures/",
                       formatC(th, width = 2, flag = "0"), "_NDVI_", hMIN, "_", hMAX, "m.pdf"),
                width = 10.88,   
                height = 9.22,
                family = "Cambria")
  
  par(
    family = "Cambria",
    mfrow = c(1, 1),
    mar = c(5, 6, 4, 1),
    las = 1,
    xpd = TRUE
  )
  
  # --- Plot NDVI median and reference curves ---
  plot(me$Date, me$Index, type = "l", lwd = 2,
       col = 1,
       xlab = "", ylab = "", ylim = c(0.4, 1.1),
       xlim = c(shade_start, shade_end), axes = FALSE)
  title(ylab = "NDVI", line = 4.5, cex.lab = 1.5, adj = 0.35)
  
  # --- Axes ---
  axis(side = 2, at = seq(0.45, 0.95, 0.1), labels = seq(0.45, 0.95, 0.1), cex.axis = 1.5)
  xticks <- seq(shade_start, shade_end, by = "month")
  axis(side = 1, at = xticks, labels = format(xticks, "%b %d"), cex.axis = 1.5, line = 0.5) # move labels slightly down
  
  # --- Frost reference ---
  rect(xleft  = as.Date("2019-05-01"),
       xright = as.Date("2019-05-29"),
       ybottom = par("usr")[3],
       ytop    = 1.11,
       col = adjustcolor("grey", alpha.f = 0.2),
       border = adjustcolor("grey", alpha.f = 1))
  
  # --- Shaded polygons for NDVI quantiles ---
  polygon(x = c(q75$Date, rev(q25$Date), q75$Date[1]),
          y = c(q75$Index, rev(q25$Index), q75$Index[1]),
          col = adjustcolor("#F2CCE2", alpha.f = 0.6),   
          border = "#F2CCE2",                          
          lty = 1,
          lwd = 2 )
  polygon(x = c(q75REF$Date, rev(q25REF$Date), q75REF$Date[1]),
          y = c(q75REF$Index, rev(q25REF$Index), q75REF$Index[1]),
          col = alpha('white', 0), border = 'black', lty = 3, lwd = 2)
  
  # --- Weekly polygons for TEMP/NDVI ---
  for (po in 1:length(date_wt_T_l)) {
    week_start <- max(date_wt_T_l[po], shade_start)
    week_end <- min(date_wt_T_l[po] + 7, shade_end)
    if (week_start < week_end) {
      polygon(x = c(week_start, week_end, week_end, week_start, week_start),
              y = c(1.1, 1.1, 1.08, 1.08, 1.1),
              col = pal[wt_T_g_perc_co[po]], border = grey(0.7))
      polygon(x = c(week_start, week_end, week_end, week_start, week_start),
              y = c(1.08, 1.08, 1.06, 1.06, 1.08),
              col = pal[wt_T_l_perc_co[po]], border = grey(0.7))
    }
  }
  for (po in 1:length(date_wt_g)) {
    week_start <- max(date_wt_g[po], shade_start)
    week_end <- min(date_wt_g[po] + 7, shade_end)
    if (week_start < week_end) {
      polygon(x = c(week_start, week_end, week_end, week_start, week_start),
              y = c(1.04, 1.04, 1.02, 1.02, 1.04),
              col = pal[wt_g_perc_co[po]], border = grey(0.7))
    }
  }
  for (po in 1:length(date_wt_l)) {
    week_start <- max(date_wt_l[po], shade_start)
    week_end <- min(date_wt_l[po] + 7, shade_end)
    if (week_start < week_end) {
      polygon(x = c(week_start, week_end, week_end, week_start, week_start),
              y = c(1.02, 1.02, 1.00, 1.00, 1.02),
              col = pal[wt_l_perc_co[po]], border = grey(0.7))
    }
  }
  
  # --- NDVI median lines ---
  lines(me$Date, me$Index, lwd = 2, col = 1)
  lines(meREF$Date, meREF$Index, lwd = 2, col = 1, lty = 2)
  
  # --- TEMP/NDVI labels ---
  par(xpd = TRUE)
  text(x = shade_start - 5, y = 1.08, labels = "TEMP", cex = 1.2, adj = 1)
  text(x = shade_start - 5, y = 1.02, labels = "NDVI", cex = 1.2, adj = 1)
  text(x = shade_start - 2, y = 1.09, labels = "G", cex = 1.2)
  text(x = shade_start - 2, y = 1.07, labels = "L", cex = 1.2)
  text(x = shade_start - 2, y = 1.03, labels = "G", cex = 1.2)
  text(x = shade_start - 2, y = 1.01, labels = "L", cex = 1.2)
  par(xpd = FALSE)
  

  
  # --- Legends ---
  legend(x = as.Date("2019-07-15"), y = 0.53,
         legend = c("Mean 2019",
                    "Mean 2018–23",
                    "2019 Q10–Q90",
                    "2018–23 Q10–Q90",
                    "Late spring frost"),
         
         col = c(1,
                 1,
                 adjustcolor("#F2CCE2", alpha.f = 0.6),
                 1,
                 adjustcolor("grey", alpha.f = 0.5)),  # frost colour
         
         lty = c(1, 2, 1, 3, 1),
         
         lwd = c(2, 2, 10, 2, 10),  # thick band for both polygons
         
         cex = 1.2,
         ncol = 1,
         bty = "n")
  
  legend("bottom",
         legend = c("0%   <= x < 10%", "10% <= x < 20%", "20% <= x < 30%", "30% <= x < 40%",
                    "40% <= x < 50%", "50% <= x < 60%", "60% <= x < 70%", "70% <= x < 80%",
                    "80% <= x < 90%", "90% <= x <= 100%"),
         col = 'black', fill = pal, pt.cex = 8, ncol = 1, cex = 1.2, bty = "n")
  
  dev.off()
  
  
}

#---------------------------------#
##### 01 Preliminary Analysis #####
#---------------------------------#

#----------------------#
###### Prelim fig ######
#----------------------#

#%%%%%%#
# Year #
#%%%%%%#

cairo_pdf("Figures/Prelim.pdf",
          width = 10,
          height = 9,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(2, 2),
  mar = c(4,4.5,2,2),
  las = 1,
  xpd = TRUE
)

plot(1, bty="l", xlab="", ylab="Green-up (DOY)", type="n", yaxt="n", ylim=c(50,170),
     xlim=c(0.5,6.5), cex.lab=1.6, axes=FALSE)
clipplot(abline(h=0,col="grey75",lty=2), xlim=c(0.65, par("usr")[2]))
wvioplot(data$Greenup[data$Year==2018], at=1, add=TRUE, col=alpha("#F2CCE2",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Year==2019], at=2, add=TRUE, col=alpha("#6CA5C6",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Year==2020], at=3, add=TRUE, col=alpha("#D495BB",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Year==2021], at=4, add=TRUE, col=alpha("#C06FA2",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Year==2022], at=5, add=TRUE, col=alpha("#AA4188",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Year==2023], at=6, add=TRUE, col=alpha("#94006E",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)

abline(h=par("usr")[4], col="white", lwd=6); abline(h=par("usr")[3], col="white", lwd=6) 
abline(v=par("usr")[1], col="white", lwd=6); abline(v=par("usr")[2], col="white", lwd=6)
axis(2, at=seq(50,170,20), cex.axis=1.2, las=1, col="grey20", col.axis="grey20", lwd=1.5, line=-1)
axis(1, at=1:6, labels=FALSE, col="grey20", lwd=1.5, line=0.5)
mtext(side=1, text=c("2018","2019","2020","2021","2022", "2023"),
      at=1:6, padj=1.5, cex=1, col="grey20", srt=45)
mtext("Year", side=1, line=3, cex=1.4, col="grey20")

#%%%%%%%%%%%#
# Elevation #
#%%%%%%%%%%%#

plot(1, bty="l", xlab="", ylab="Green-up (DOY)", type="n", yaxt="n", ylim=c(50,170),
     xlim=c(0.5,5.5), cex.lab=1.6, axes=FALSE)
wvioplot(data$Greenup[data$Elevation_Band=="<800"], at=1, add=TRUE, col=alpha("#ACCCE4",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Elevation_Band=="800-1000"], at=2, add=TRUE, col=alpha("#7FABD3",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Elevation_Band=="1000-1200"], at=3, add=TRUE, col=alpha("#5087C1",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Elevation_Band=="1200-1400"], at=4, add=TRUE, col=alpha("#325FA2",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(data$Greenup[data$Elevation_Band==">1400"], at=5, add=TRUE, col=alpha("#273871",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)

abline(h=par("usr")[4], col="white", lwd=6); abline(h=par("usr")[3], col="white", lwd=6) 
abline(v=par("usr")[1], col="white", lwd=6); abline(v=par("usr")[2], col="white", lwd=6)

axis(2, at=seq(50,170,20), cex.axis=1.2, las=1, col="grey20", col.axis="grey20", lwd=1.5, line=-1)

# Multiline, aligned labels using expression(atop) and adjusted line
axis(1, at=1:5, labels=FALSE)  # draw tick marks only
mtext(side=1, at=1:5, text=c(
  "<800",
  expression(atop("800–","1000")),
  expression(atop("1000–","1200")),
  expression(atop("1200–","1400")),
  ">1400"
), line=2.5, cex=1, col="grey20")  # adjust 'line' to move down
mtext("Elevation band", side=1, line=4.2, cex=1.4, col="grey20")

#%%%%%%%%#
# Aspect #
#%%%%%%%%#

# Bin cosine-transformed Aspect and calculate mean ± SD
binned <- data %>%
  mutate(Aspect_cos_bin = cut(Aspect, breaks = seq(-1, 1, 0.05), include.lowest = TRUE)) %>%
  group_by(Aspect_cos_bin) %>%
  summarise(
    mean_greenup = mean(Greenup, na.rm = TRUE),
    sd_greenup = sd(Greenup, na.rm = TRUE)
  ) %>%
  mutate(
    Aspect_mid = sapply(strsplit(as.character(Aspect_cos_bin), ","), function(x) {
      as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.025  # midpoint of bin
    })
  )

# Base R plot
plot(1, type="n", xlim=c(-1,1), ylim=c(100,140),
     xlab = expression(Aspect[cos]), ylab="Mean Green-up (DOY)",
     bty="l", cex.lab=1.6, axes=FALSE)

# Add ribbon (mean ± SD)
polygon(
  x = c(binned$Aspect_mid, rev(binned$Aspect_mid)),
  y = c(binned$mean_greenup - binned$sd_greenup,
        rev(binned$mean_greenup + binned$sd_greenup)),
  col = "#D5E8F3", border = NA
)

# Add mean line
lines(binned$Aspect_mid, binned$mean_greenup, lwd=3, col="#6CA5C6")

# Axes
axis(1, at=seq(-1,1,0.5), labels=seq(-1,1,0.5), col="grey20", col.axis="grey20", lwd=1.5)
axis(2, at=seq(100,140,10), col="grey20", col.axis="grey20", lwd=1.5, las=1)

dev.off()

#---------------------#
###### Bar Plots ######
#---------------------#

#------------------------------------#
##### Check distribution of data #####
#------------------------------------#

#/////////////////////////////////#
# Observations per elevation band #
#/////////////////////////////////#

# Count number of observations per elevational band
elevational_band_counts <- table(data$Elevation_Band)

# Order elevational bands from smallest to largest
elev_order <- c("<800", "800-1000", "1000-1200", "1200-1400", ">1400")
elevational_band_counts <- elevational_band_counts[elev_order]

#///////////////////////#
# Observations per year #
#///////////////////////#

# Count the number of observations 
year_counts <- table(data$Year)

# Order years from smallest to largest
year_order <- c("2018", "2019", "2020", "2021", "2022", "2023")
year_counts <- year_counts[year_order]


# Define colors
bar_colors_elev <- alpha(c("#ACCCE4", "#7FABD3", "#5087C1", "#325FA2", "#273871"), 0.6)
bar_colors_year <- alpha(c("#F2CCE2", "#6CA5C6", "#D495BB", "#C06FA2", "#AA4188", "#94006E"), 0.6)

cairo_pdf("Figures/Sample_Distribution.pdf",
          width = 9,
          height = 5,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 2),
  mar = c(6,5,2,2),
  las = 1,
  xpd = TRUE,
  cex.axis = 0.8
)

#%%%%%%%%%%%%%%%%%%%%#
# Elevation barplot #
#%%%%%%%%%%%%%%%%%%%%#

bp1 <- barplot(elevational_band_counts,
               col = bar_colors_elev,
               border = "grey20",
               main = "",
               xlab = "",
               ylab = "",
               ylim = c(0, max(elevational_band_counts)*1.1),
               names.arg = FALSE)  # suppress labels

# Custom multiline labels
mtext(side = 1, at = bp1, text = c(
  "<800",
  expression(atop("800–","1000")),
  expression(atop("1000–","1200")),
  expression(atop("1200–","1400")),
  ">1400"
), line = 2.5, cex = 0.8, col = "grey20")

mtext("Elevational band", side = 1, line = 4, cex = 1)

# Y label
mtext("Number of Observations", side = 2, line = 4, las = 0, col = "black", cex = 1)

#%%%%%%%%%%%%%%%%#
# Year barplot #
#%%%%%%%%%%%%%%%%#

bp2 <- barplot(year_counts,
               col = bar_colors_year,
               border = "grey20",
               main = "",
               xlab = "",
               ylab = "",
               ylim = c(0, 50000),
               names.arg = FALSE)

# Compute y-offset
yrange <- diff(par("usr")[3:4])

# Add rotated labels with offset
text(x = bp2,
     y = par("usr")[3] - 0.06* yrange, 
     labels = names(year_counts),
     srt = 45,
     adj = c(1, 1),
     xpd = TRUE,
     cex = 0.8,
     col = "grey20")
mtext("Year", side = 1, line = 4, cex = 1)

# Y label
mtext("Number of Observations", side = 2, line = 4, las = 0, col = "black", cex = 1)

dev.off()

#-------------------------#
##### 02 Hypothesis 1 #####
#-------------------------#

#//////////////////////////#
###### Pre-frost NDVI ######
#//////////////////////////#

cairo_pdf("Figures/H1prefrost_aspect.pdf",
          width = 10,
          height = 5,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 2),
  mar = c(4,4.5,1,1),
  oma = c(2,0,2,0),
  las = 1,
  xpd = TRUE
)

# Ensure factor
df_filtered$Elevation_Band <- factor(df_filtered$Elevation_Band)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Pre-frost NDVI (marginal) #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

x_seq <- seq(min(df_filtered$mean_pre_frost, na.rm = TRUE),
             max(df_filtered$mean_pre_frost, na.rm = TRUE),
             length.out = 200)

x_seq_c <- x_seq - mean(df_filtered$mean_pre_frost, na.rm = TRUE)

elev_levels <- levels(df_filtered$Elevation_Band)

# Store predictions
fit_mat <- matrix(NA, nrow = length(x_seq), ncol = length(elev_levels))
se_mat  <- matrix(NA, nrow = length(x_seq), ncol = length(elev_levels))

for(i in seq_along(elev_levels)){
  
  newdat <- data.frame(
    mean_pre_frost_c  = x_seq_c,
    mean_pre_frost_c2 = x_seq_c^2,
    Aspect            = mean(df_filtered$Aspect, na.rm = TRUE),
    Elevation_Band    = factor(elev_levels[i], levels = elev_levels)
  )
  
  pred <- predict(quad_model_centered, newdata = newdat, se.fit = TRUE)
  
  fit_mat[, i] <- pred$fit
  se_mat[, i]  <- pred$se.fit
}

# Average across elevation bands
ndvi_pred <- rowMeans(fit_mat)
ndvi_se   <- sqrt(rowMeans(se_mat^2))

ndvi_lwr <- ndvi_pred - 1.96 * ndvi_se
ndvi_upr <- ndvi_pred + 1.96 * ndvi_se

# Plot NDVI
plot(df_filtered$mean_pre_frost, df_filtered$ndvi_change,
     xlab = "Pre-frost NDVI",
     ylab = "NDVI Change",
     pch = 16,
     col = adjustcolor("#273871", alpha.f = 0.05),
     cex.lab = 1.3,
     bty = "l")

polygon(c(x_seq, rev(x_seq)),
        c(ndvi_lwr, rev(ndvi_upr)),
        col = adjustcolor("#ACCCE4", alpha.f = 0.4),
        border = NA)

lines(x_seq, ndvi_pred, col = '#6CA5C6', lwd = 2)

#%%%%%%%%#
# Aspect #
#%%%%%%%%#

aspect_seq <- seq(min(df_filtered$Aspect, na.rm = TRUE),
                  max(df_filtered$Aspect, na.rm = TRUE),
                  length.out = 200)

emm_aspect <- emmeans(
  quad_model_centered,
  ~ Aspect,
  at = list(
    Aspect = aspect_seq,
    mean_pre_frost_c = 0,
    mean_pre_frost_c2 = 0
  )
)

aspect_vals <- summary(emm_aspect)

aspect_pred <- aspect_vals$emmean
aspect_lwr  <- aspect_vals$lower.CL
aspect_upr  <- aspect_vals$upper.CL

# Plot Aspect
plot(df_filtered$Aspect, df_filtered$ndvi_change,
     xlab = expression(Aspect[plain(cos)]),
     ylab = "NDVI Change",
     pch = 16,
     col = adjustcolor("#273871", alpha.f = 0.05),
     cex.lab = 1.3,
     bty = "l")

polygon(c(aspect_seq, rev(aspect_seq)),
        c(aspect_lwr, rev(aspect_upr)),
        col = adjustcolor("#ACCCE4", alpha.f = 0.4),
        border = NA)

lines(aspect_seq, aspect_pred, col = '#6CA5C6', lwd = 2)

dev.off()

#-------------------------#
###### By Elevation  ######
#-------------------------#

cairo_pdf("Figures/H1elevation.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 1),
  mar=c(6,5,2,2),
  las = 1,
  xpd = TRUE
)

plot(1, bty="l", xlab="", ylab="NDVI Change", type="n", yaxt="n", ylim=c(-0.3,0.7),
     xlim=c(0.5,5.5), cex.lab=1.2, axes=FALSE)
wvioplot(df_filtered$ndvi_change[df_filtered$Elevation_Band=="<800"], at=1, add=TRUE, col=alpha("#ACCCE4",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(df_filtered$ndvi_change[df_filtered$Elevation_Band=="800-1000"], at=2, add=TRUE, col=alpha("#7FABD3",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(df_filtered$ndvi_change[df_filtered$Elevation_Band=="1000-1200"], at=3, add=TRUE, col=alpha("#5087C1",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(df_filtered$ndvi_change[df_filtered$Elevation_Band=="1200-1400"], at=4, add=TRUE, col=alpha("#325FA2",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)
wvioplot(df_filtered$ndvi_change[df_filtered$Elevation_Band==">1400"], at=5, add=TRUE, col=alpha("#273871",0.6),
         border="grey20", wex=0.6, drawRect=TRUE, lwd=1.5, adjust=1)

abline(h=par("usr")[4], col="white", lwd=6); abline(h=par("usr")[3], col="white", lwd=6) 
abline(v=par("usr")[1], col="white", lwd=6); abline(v=par("usr")[2], col="white", lwd=6)

axis(2, at=seq(-0.3,0.7,0.2), cex.axis=1, las=1, col="grey20", col.axis="grey20", lwd=1.5, line=-1)

# Multiline, aligned labels using expression(atop) and adjusted line
axis(1, at=1:5, labels=FALSE)  # draw tick marks only
mtext(side=1, at=1:5, text=c(
  "<800",
  expression(atop("800–","1000")),
  expression(atop("1000–","1200")),
  expression(atop("1200–","1400")),
  ">1400"
), line=2.5, cex=1, col="grey20")  # adjust 'line' to move down
mtext("Elevation band", side=1, line=4.2, cex=1.2, col="grey20")

dev.off()

#---------------------------#
##### 05 Hypothesis 2/3 #####
#---------------------------#

#----------------------#
###### Line plot #######
#----------------------#

data_with_ndvi_2019 <- data_with_ndvi[data_with_ndvi$Year == 2019, ]

cairo_pdf("Figures/H2NDVILine.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 1),
  mar=c(6,5,2,2),
  las = 1,
  xpd = TRUE
)

# Summarize mean + SE NDVI per Elevation_Band and Frost_Class
ndvi_summary <- data_with_ndvi_2019 %>%
  group_by(Elevation_Band, Frost_Class) %>%
  summarise(
    mean_NDVI = mean(Mean_NDVI, na.rm = TRUE),
    sd_NDVI = sd(Mean_NDVI, na.rm = TRUE),
    n = n(),
    se_NDVI = sd_NDVI / sqrt(n),
    .groups = "drop"
  )

# Plot lines with error ribbons
ggplot(ndvi_summary, aes(x = Elevation_Band, y = mean_NDVI, color = Frost_Class, group = Frost_Class)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = mean_NDVI - se_NDVI, ymax = mean_NDVI + se_NDVI, fill = Frost_Class),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Frost" = "#6CA5C6", "Non-frost" = "#D495BB")) +
  scale_fill_manual(values = c("Frost" = "#6CA5C6", "Non-frost" = "#D495BB")) +
  scale_x_discrete(
    limits = c("<800", "800-1000", "1000-1200", "1200-1400", ">1400"),
    labels = c(
      "<800",
      "800–\n1000",
      "1000–\n1200",
      "1200–\n1400",
      ">1400"
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 12)),
    axis.title.y = element_text(size = 13, margin = margin(r = 12)),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )+
  labs(
    x = "Elevation Band (m)",
    y = "Mean NDVI",
    color = "Frost class",
    fill = "Frost class"
  )

dev.off()

#------------------------#
###### Violin plot #######
#------------------------#

# Base R approach
ndvi2018_long <- weekly_NDVI[weekly_NDVI$Year == 2018, ]
ndvi2019_long <- weekly_NDVI[weekly_NDVI$Year == 2019, ]
ndvi2020_long <- weekly_NDVI[weekly_NDVI$Year == 2020, ]
ndvi2021_long <- weekly_NDVI[weekly_NDVI$Year == 2021, ]
ndvi2022_long <- weekly_NDVI[weekly_NDVI$Year == 2022, ]
ndvi2023_long <- weekly_NDVI[weekly_NDVI$Year == 2023, ]

# Original datasets in a list
ndvi_list <- list(ndvi2018_long, ndvi2019_long, ndvi2020_long, 
                  ndvi2021_long, ndvi2022_long, ndvi2023_long)

# Names for each dataset
names(ndvi_list) <- c("ndvi2018_long", "ndvi2019_long", "ndvi2020_long", 
                      "ndvi2021_long", "ndvi2022_long", "ndvi2023_long")

# Filter each dataset to IQR and overwrite in the list
ndvi_list <- lapply(ndvi_list, function(df) {
  Q1 <- quantile(df$Mean_NDVI, 0.25, na.rm = TRUE)
  Q3 <- quantile(df$Mean_NDVI, 0.75, na.rm = TRUE)
  df[df$Mean_NDVI >= Q1 & df$Mean_NDVI <= Q3, ]
})

# Assign back to original variables
list2env(ndvi_list, envir = .GlobalEnv)

# Filter for frost
filter_frost <- function(df) {
  df %>%
    filter(format(Date, "%m-%d") >= "05-01" &
             format(Date, "%m-%d") <= "05-29")
}

ndvi2018_long <- filter_frost(ndvi2018_long)
ndvi2019_long <- filter_frost(ndvi2019_long)
ndvi2020_long <- filter_frost(ndvi2020_long)
ndvi2021_long <- filter_frost(ndvi2021_long)
ndvi2022_long <- filter_frost(ndvi2022_long)
ndvi2023_long <- filter_frost(ndvi2023_long)

range(ndvi2018_long$Mean_NDVI)
range(ndvi2019_long$Mean_NDVI)
range(ndvi2020_long$Mean_NDVI)
range(ndvi2021_long$Mean_NDVI)
range(ndvi2022_long$Mean_NDVI)
range(ndvi2023_long$Mean_NDVI)

cairo_pdf("Figures/NDVI_violin.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 1),
  mar=c(6,6,6,2),
  las = 1,
  xpd = TRUE
)

plot(1,bty="l",xlab="Year",ylab="Mean NDVI", type="n",yaxt="n",ylim=c(0.74,0.9),xlim=c(0.5,12),cex.lab=1.6,axes=F)
wvioplot(ndvi2018_long$Mean_NDVI[ndvi2018_long$Frost_Class=="Non-frost"],at=1,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2018_long$Mean_NDVI[ndvi2018_long$Frost_Class=="Frost"],at=1.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2019_long$Mean_NDVI[ndvi2019_long$Frost_Class=="Non-frost"],at=3,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2019_long$Mean_NDVI[ndvi2019_long$Frost_Class=="Frost"],at=3.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2020_long$Mean_NDVI[ndvi2020_long$Frost_Class=="Non-frost"],at=5,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2020_long$Mean_NDVI[ndvi2020_long$Frost_Class=="Frost"],at=5.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2021_long$Mean_NDVI[ndvi2021_long$Frost_Class=="Non-frost"],at=7,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2021_long$Mean_NDVI[ndvi2021_long$Frost_Class=="Frost"],at=7.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2022_long$Mean_NDVI[ndvi2022_long$Frost_Class=="Non-frost"],at=9,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2022_long$Mean_NDVI[ndvi2022_long$Frost_Class=="Frost"],at=9.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2023_long$Mean_NDVI[ndvi2023_long$Frost_Class=="Non-frost"],at=11,add=T,
         col=alpha("#D495BB",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(ndvi2023_long$Mean_NDVI[ndvi2023_long$Frost_Class=="Frost"],at=11.8,add=T,
         col=alpha("#6CA5C6",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)

# these lines get rid of extra lines, make them the same colour as background:
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6) 
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
# Define the positions of each pair (Non-frost + Frost)
x_pos <- c(1, 1.8, 3, 3.8, 5, 5.8, 7, 7.8, 9, 9.8, 11, 11.8)
# Define the positions where you want year labels (center of each pair)
year_pos <- c(1.4, 3.4, 5.4, 7.4, 9.4, 11.4)  # midpoints between Non-frost and Frost
# Define the labels
years <- 2018:2023
# Add x-axis with your year labels at the center of each pair
axis(1, at = year_pos, labels = years, cex.axis = 1.2, col = "grey20", col.axis = "grey20", lwd = 1.5)
axis(2,at=seq(0.74,0.9,0.04),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)

# Add legend
legend("topright",
       legend = c("Non-frost", "Frost"),
       fill = c(alpha("#D495BB",0.6), alpha("#6CA5C6",0.6)),
       border = "grey20",
       bty = "n",
       cex = 1.2,
       inset = c(0.02, -0.2))  # move *up* by decreasing y (negative moves outward)


dev.off()

#----------------------#
###### Line plot #######
#----------------------#

weekly_NDVI[, DOY := yday(Date)]
head(weekly_NDVI)
NDVI_frost <- weekly_NDVI[Frost_Class == "Frost"]

cairo_pdf("Figures/H3Legacy_frost.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(1, 1),
  mar=c(8,6,2,2),
  las = 1,
  xpd = TRUE
)

ggplot(NDVI_frost, aes(x = DOY, y = Mean_NDVI, color = factor(Year), group = Year, linetype = factor(Year))) +
  annotate(
    "rect",
    xmin = 121,
    xmax = 149,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey",
    alpha = 0.2,
    colour = "grey"
  ) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  scale_color_manual(values = c(
    "2018" = "#F2CCE2",
    "2019" = "#6CA5C6", 
    "2020" = "#D495BB",
    "2021" = "#C06FA2", 
    "2022" = "#AA4188", 
    "2023" = "#94006E"
  )) +
  scale_linetype_manual(values = c(
    "2018" = "solid",
    "2019" = "solid",
    "2020" = "dotted",
    "2021" = "twodash",
    "2022" = "longdash",
    "2023" = "solid"
  )) +
  labs(
    x = "Day of Year",
    y = "Mean NDVI",
    color = "Year",
    linetype = "Year"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 15)),
    axis.title.y = element_text(size = 11, margin = margin(r = 15)),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

dev.off()