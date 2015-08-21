require(tikzDevice)
source("../../basic_functions.R")

tex_file = "gm.tex"

# Load experiment information
algo_config <- readConfig("../algorithms.config",c("ALGO_ID","LATEX-NAME","PCH","LTY","COL"))
pattern_config <- readConfig("../patterns.config",c("TC_ID","SP-LEN","GAP"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the query time figure
plot_gm_query_times <- function( data, title=""){
  cat("Graph: ",title,"\n")
  #data <- aggregate(data[c('total_time_ms')], by=c(data['PATT_SAMPLE']), FUN=min)
  #data[c('total_time_ms')] <- data[c('total_time_ms')]/1000.0
  #data <- data[order(data[['PATT_SAMPLE']]), ]

  runtime <- data[['total_time_ms']] # / data[['num_results']]
  max_runtime <- max(runtime[!is.infinite(runtime) & !is.nan(runtime)])
  max_sample <- length(unique(data[['PATT_SAMPLE']]))

  plot(NA, NA, type="l", xlim=c(0, max_sample ), ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", xaxt="n")

  ALGO_IDs <- unique(data$ALGO)
  for(algo in ALGO_IDs){
    d <- data[data$ALGO==algo,]
    lines(d[['PATT_SAMPLE']], d[['total_time_ms']] # / d[['num_results']]
          ,lwd=1, type="b", 
          pch=algo_config[algo,"PCH"], 
          lty=algo_config[algo,"LTY"],
          col=algo_config[algo,"COL"])
    box("plot", col="grey")    
    axis( 1, at = axTicks(1), labels=T, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
  }

  mtext("Test case", side=1, line=2, las=0)
  axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
  mtext("Total query time in ($ms$)", side=2, line=2, las=0)
  #mtext("Average query time in ($\\mu s$)", side=2, line=2, las=0)

  grid(lty=1)  
  draw_figure_heading(sprintf("collection = %s",title))

  plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
  legend( "top", legend=algo_config[ALGO_IDs,"LATEX-NAME"], pch=algo_config[ALGO_IDs,"PCH"], col=algo_config[ALGO_IDs,"COL"],
            lty=algo_config[ALGO_IDs,"LTY"], bty="n", y.intersp=1.5, ncol=2, title="Algorithm", cex=1.2)
}

get_sa_range_for <- function(data, coll, pattern_id) {
  d <- data[data$PATT_SAMPLE==pattern_id,]
  numbers <- as.numeric(unlist(strsplit(paste(d[["info"]]),",")))
  numbers <- numbers[!is.na(numbers)]
  return(sum(numbers));
}

create_table_for <- function(data, coll, algo) {
  d <- data[data$ALGO==algo,]
  cat("Table: ",coll, " ", algo,"\n")
  # layout: rows=patterns, cols=values of interest

  PATT_IDs <- unique(d$PATT_SAMPLE)
  fig_name <- paste("fig-gm-time-",coll,"-",algo,"-table.tex",sep="")
  sink(fig_name)
  cat("\\begin{tabular}{|l|", paste(rep("r", length(PATT_IDs)), sep=" "), "|}\n", sep="")
  cat("\\hline\n")

  table = c("pattern",
            "phrase length",
            "gap size",
            "mean time ($ms$)",
            "SA range (\\#potential matches)",
            "mean time per potential match ($\\mu s$)")
            
  for(patt_id in PATT_IDs){
    dd <- d[d$PATT_SAMPLE==patt_id,]
    sa_range <- get_sa_range_for(data,coll,patt_id)
    table = paste(table, 
             c(patt_id,
               pattern_config[patt_id,"SP-LEN"],
               pattern_config[patt_id,"GAP"],
               paste("$", dd[["mean_time_ms"]], "$", sep=""),
               paste("$", sa_range, "$", sep=""),
               paste("$", sprintf("%.3f",1000 * dd[["mean_time_ms"]] / sa_range), "$", sep="")
             ), sep=" & ")
  }

  cat(table, "\\hline", sep="\\\\\n")
  cat("\\end{tabular}\n")
  sink(NULL)
  
  cat(table)
  
  return(paste("\\begin{figure}
                \\input{",fig_name,"}
                \\caption{Query times of algorithm \\texttt{",algo_config[algo,"LATEX-NAME"],"} on \\texttt{",coll,"}.}
                \\end{figure}", sep=""))
}

data <- data_frame_from_key_value_pairs( "../results/all.txt" )
tex_doc <- paste(readLines("gm-header.tex"),collapse="\n")

colls <- unique(data$COLL_ID)
n <- length(colls)
for(coll in colls){
  coll_name <- 
  
  # Graph for total query time
  fig_name <- paste("fig-gm-time-",coll,".tex",sep="")
  open_tikz( fig_name )

  multi_figure_style( 1, 1 )  

  d <- data[data$COLL_ID==coll,]
  plot_gm_query_times(d, title=coll)
  dev.off()
  tex_doc <- paste(tex_doc,"\\begin{figure}
               \\input{",fig_name,"}
               \\caption{Query time on \\texttt{",coll,"} depending on gap size.
               }
              \\end{figure}")
              
  # tables for mean query times and result numbers
  ALGO_IDs <- unique(d$ALGO)
  for(algo in ALGO_IDs){
    tex_doc <- paste(tex_doc,create_table_for(d, coll, algo))
  }
}

tex_doc <- paste(tex_doc, readLines("gm-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)


