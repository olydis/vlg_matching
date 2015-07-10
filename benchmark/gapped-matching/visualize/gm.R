require(tikzDevice)
source("../../basic_functions.R")

tex_file = "gm.tex"

# Load experiment information
tc_config <- readConfig("../test_case.config",c("TC_ID","LATEX-NAME"))
compile_config <- readConfig("../compile_options.config",c("COMPILE_ID","OPTIONS"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the query time figure
plot_gm_query_times <- function( data, max_y=NA, title="", yaxis=T, xaxis=T){
  cat(title,"\n")
print(data[c('total_time_ms')])
  data <- aggregate(data[c('total_time_ms')], by=c(data['PATT_SAMPLE']), FUN=min)
  data[c('total_time_ms')] <- data[c('total_time_ms')]/1000.0
  data <- data[order(data[['PATT_SAMPLE']]), ]

  max_runtime <- data[['total_time_ms']]
  if ( !is.na(max_y) ){
    max_runtime = max_y
  }

  plot(data[['COLL_PATT_ID']], data[['total_time_ms']], type="l", ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", cex.axis = 0.8, xaxt="n",
       col = terrain.colors(6)[1], lwd=1.5
       )
  box("plot", col="grey")    
  axis( 1, at = axTicks(1), labels=xaxis, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
  if ( xaxis ){
    xlable <- "Block size K"
    mtext(xlable, side=1, line=2, las=0)
  }
  if ( yaxis ){
    axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
    mtext("Time per operation in ($\\mu s$)", side=2, line=2, las=0)
  }
  grid(lty=1)
  lines( data[['K']], data[['rank_time']], lty=1, col=terrain.colors(6)[3], lwd=1.5)    
  lines( data[['K']], data[['select_time']], lty=1, col=terrain.colors(6)[5], lwd=1.5)    

  draw_figure_heading(sprintf("bitvector = %s",title))
}

data <- data_frame_from_key_value_pairs( "../results/all.txt" )
tex_doc <- paste(readLines("gm-header.tex"),collapse="\n")

for ( compile_id in compile_config[["COMPILE_ID"]] ){
# Handle query time
  fig_name <- paste("fig-gm-time-",compile_id,".tex",sep="")
  open_tikz( fig_name )
  n <- nrow(tc_config)
  d <- subset(data, data[["COMPILE_ID"]]==compile_id)
  par(mfrow=c(n/2+1,2))
  multi_figure_style( n/2+1, 2 )  
  nr <- 0
  xlabnr <- 2*(n/2)-2+1
  for ( tc_id in tc_config[["TC_ID"]] ){
    dd <- subset(d, d[["TC_ID"]]==tc_id)
    plot_gm_query_times(dd, 2.1, 
                         title=tc_config[tc_id,"LATEX-NAME"],
                            yaxis=(nr%%2==0), xaxis=(nr>=xlabnr) )
    if ( nr == 0 ){
        plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
        legend("topleft", legend=rev(c("access","rank","select")), box.lwd=0, lty=rev(c(1,1,1)), 
              title="Operation", col=rev(terrain.colors(6)[seq(1,5,2)]), bg="white", cex=1.5)
        nr <- nr+1
    }
    nr <-nr+1 
  }
  dev.off()
  tex_doc <- paste(tex_doc,"\\begin{figure}
               \\input{",fig_name,"}
               \\caption{Runtime of \\texttt{gm\\_vector} dependent on block size K.
               Compile options:
               \\texttt{",gsub("_","\\\\_",compile_config[compile_id, "OPTIONS"]),"}.
               }
              \\end{figure}")
}

tex_doc <- paste(tex_doc, readLines("gm-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)


