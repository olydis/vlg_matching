require(tikzDevice)
source("../../basic_functions.R")

tex_file = "gm.tex"

# Load experiment information
tc_config <- readConfig("../test_case.config",c("TC_ID","LATEX-NAME"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the query time figure
plot_gm_query_times <- function( data, title=""){
  cat(title,"\n")
  #data <- aggregate(data[c('total_time_ms')], by=c(data['PATT_SAMPLE']), FUN=min)
  #data[c('total_time_ms')] <- data[c('total_time_ms')]/1000.0
  #data <- data[order(data[['PATT_SAMPLE']]), ]

  max_runtime <- max(data[['total_time_ms']])

  for(tc in unique(data$TC_ID)){
    d <- data[data$TC_ID==tc,]
    plot(data[['PATT_SAMPLE']], data[['total_time_ms']], type="l", ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", cex.axis = 0.8, xaxt="n",
       col = terrain.colors(6)[1], lwd=1.5
       )
    box("plot", col="grey")    
    axis( 1, at = axTicks(1), labels=xaxis, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )

  }

    xlable <- "Gap size"
    mtext(xlable, side=1, line=2, las=0)

    axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
    mtext("Query time in ($\\mu s$)", side=2, line=2, las=0)

  grid(lty=1)  

  #draw_figure_heading(sprintf("collection = %s",title))
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
               \\caption{Query time on \\texttt{",coll,"} dependent on gap size..
               }
              \\end{figure}")
}

tex_doc <- paste(tex_doc, readLines("gm-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)


