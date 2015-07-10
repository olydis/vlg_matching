require(tikzDevice)
source("../../basic_functions.R")

tex_file = "gm.tex"

# Load experiment information
tc_config <- readConfig("../test_case.config",c("TC_ID","LATEX-NAME","PCH","LTY","COL"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the query time figure
plot_gm_query_times <- function( data, title=""){
  cat(title,"\n")
  #data <- aggregate(data[c('total_time_ms')], by=c(data['PATT_SAMPLE']), FUN=min)
  #data[c('total_time_ms')] <- data[c('total_time_ms')]/1000.0
  #data <- data[order(data[['PATT_SAMPLE']]), ]

  max_runtime <- max(data[['mean_time_ms']])
  max_sample <- max(data[['PATT_SAMPLE']])

    #plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
    plot(NA, NA, type="l", xlim=c(0, max_sample ), ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", xaxt="n"
       )

  tc_ids <- unique(data$TC_ID)
  for(tc in tc_ids){
    d <- data[data$TC_ID==tc,]
    lines(d[['PATT_SAMPLE']], d[['mean_time_ms']],
          lwd=1, type="b", 
          pch=tc_config[tc,"PCH"], 
          lty=tc_config[tc,"LTY"],
          col=tc_config[tc,"COL"])
    box("plot", col="grey")    
    axis( 1, at = axTicks(1), labels=T, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )

  }


    mtext("Gap size", side=1, line=2, las=0)
    axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
    mtext("Average query time in ($ms$)", side=2, line=2, las=0)

  grid(lty=1)  
  draw_figure_heading(sprintf("collection = %s",title))

  plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
  legend( "top", legend=tc_config[tc_ids,"LATEX-NAME"], pch=tc_config[tc_ids,"PCH"], col=tc_config[tc_ids,"COL"],
		    lty=tc_config[tc_ids,"LTY"], bty="n", y.intersp=1.5, ncol=2, title="Algorithm", cex=1.2)
}

data <- data_frame_from_key_value_pairs( "../results/all.txt" )
tex_doc <- paste(readLines("gm-header.tex"),collapse="\n")

colls <- unique(data$COLL_ID)
n <- length(colls)
for(coll in colls){
# Handle query time
  fig_name <- paste("fig-gm-time-",coll,".tex",sep="")
  open_tikz( fig_name )

  par(mfrow=c(1,1))
  multi_figure_style( 1, 1 )  

  d <- data[data$COLL_ID==coll,]
  plot_gm_query_times(d, title=coll)
  dev.off()
  tex_doc <- paste(tex_doc,"\\begin{figure}
               \\input{",fig_name,"}
               \\caption{Query time on \\texttt{",coll,"} depending on gap size.
               }
              \\end{figure}")
}

tex_doc <- paste(tex_doc, readLines("gm-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)


