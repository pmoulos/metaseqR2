#' MDS plot JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diagplotMds}}.
#'
#' @param obj A list holding MDS plot data. See \code{\link{diaplot.mds}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
mdsToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    x <- obj$x
    y <- obj$y
    xlim <- obj$xlim
    ylim <- obj$ylim
    samples <- obj$samples
    cols <- .getColorScheme()
    
    # Construct series
    counter <- 0
    series <- vector("list",length(samples))
    names(series) <- names(samples)
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        series[[n]]$type="scatter"
        series[[n]]$color=cols$fill[counter]
        series[[n]]$marker <- list(
            lineWidth=1,
            states=list(
                hover=list(
                    enabled=TRUE,
                    lineColor=cols$border[counter]
                ),
                select=list(
                    fillColor=cols$selected[counter],
                    lineColor=cols$border[counter],
                    lineWidth=2
                )
            )
        )
        m <- match(samples[[n]],names(x))
        if (length(m)>0) {            
            series[[n]]$data <- makeHighchartsPoints(x[m],y[m])
        }
    }
    
    switch(jl,
        highcharts = {
            point.format=paste("<strong>Sample name: </strong>{point.name}<br>",
                "<strong>Principal coordinate 1: </strong>{point.x}<br>",
                "<strong>Principal coordinate 2: </strong>{point.y}",sep="")
                
                json <- toJSON(list(
                    chart=list(
                        type="scatter",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("Multidimensional Scaling")
                    ),
                    xAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="1<sup>st</sup> Principal Coordinate",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(xlim[1],3),
                        max=round(xlim[2],3)
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="2<sup>nd</sup> Principal Coordinate",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(ylim[1],3),
                        max=round(ylim[2],3)
                    ),
                    plotOptions=list(
                        scatter=list(
                            allowPointSelect=TRUE,
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=point.format
                            )
                        )
                    ),
                    series=unname(series)
                ))
        }
    )
    return(json)
}

countsBioToJSON <- function(obj,by=c("sample","biotype"),jl=c("highcharts")) {
    by <- tolower(by[1])
    jl <- tolower(jl[1])
    samples <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    counts <- round(2^obj$user$counts - 1)
    counts[counts==0] <- 0.001
    #counts <- obj$user$counts
    
    covars <- obj$user$covars
    biotypes <- unique(as.character(covars$biotype))
    if (!is.null(altnames))
        names(altnames) <- rownames(counts)
    
    grouped <- FALSE
    if (is.null(samples)) {
        if (is.null(colnames(counts)))
            samplenames <- paste("Sample",1:ncol(counts),sep=" ")
        else
            samplenames <- colnames(counts)
        samples <- list(Samples=nams)
    }
    else if (is.list(samples)) {
        samplenames <- unlist(samples,use.names=FALSE)
        grouped <- TRUE
    }
    
    # y label formatter for logarithmic axis
    y.label.formatter <- paste('function() {if(this.value === 0.001)',
        '{return 0;} else {return Highcharts.Axis.prototype.',
        'defaultLabelFormatter.call(this);}}',sep="")
        
    tooltip.point.formatter <- paste("function() {",
        "   var min = this.low === 0.001 ? 0 : this.low;" ,
        "   var q1 = this.q1 === 0.001 ? 0 : this.q1;" ,
        "   var med = this.median === 0.001 ? 0 : this.median;",
        "   var q3 = this.q3 === 0.001 ? 0 : this.q3;",
        "   var max = this.high === 0.001 ? 0 : this.high;",
        "   var str = 'Maximum: ' + max + '<br/>' +",
        "       'Upper quartile: ' + q3 + '<br/>' +",
        "       'Median: ' + med + '<br/>' +",
        "       'Lower quartile: ' + q1 + '<br/>' +",
        "       'Minimum: ' + min + '<br/>';",
        "   return  str;",
        "}",sep="")
        
    # Legend clicker
    boxplot.onclick <- paste("function() { ",
        "   var chart =  this.chart;",
        "   var outlier_id =  chart.get(this.name);",
        "   if (!outlier_id.visible) {",
        "       outlier_id.show();",
        "   } else {",
        "       outlier_id.hide();",
        "   }",
        "}",sep="")
    
    # Outliers tooltip
    if (is.null(obj$altnames)) {
        outlier.pointformat <- paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    else {
        outlier.pointformat <- paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Gene name: {point.alt_name}<br/>',
            'Value: {point.y}<br/>',
            sep=""
        )
    }
    
    if (by=="sample") {
        cols <- .getColorScheme(length(biotypes))
        box.list <- json <- vector("list",length(samplenames))
        names(box.list) <- names(json) <- samplenames
        for (n in samplenames) {
            box.list[[n]] <- vector("list",length(biotypes))
            names(box.list[[n]]) <- biotypes
            for (b in biotypes)
                box.list[[n]][[b]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(box.list[[n]],plot=FALSE)$stats
            colnames(B) <- biotypes
            oList <- lapply(names(box.list[[n]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },box.list[[n]],B)
            names(oList) <- biotypes
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in 1:ncol(B))
                BB[,jj] <- round(B[,jj],3)
            d <- as.data.frame(BB)
            ids <- 0:(ncol(d)-1)
            d <- rbind(ids,d)
            names(ids) <- colnames(d)
            counter <- 0
            series <- vector("list",length(biotypes))
            names(series) <- biotypes
            for (s in names(series)) {
                counter <- counter + 1
                series[[s]] <- list()
                series[[s]]$name <- s
                series[[s]]$color <- cols$fill[counter]
                #series[[s]]$turboThreshold <- 10000
                series[[s]]$data <- list(unname(as.list(d[,s])))
                r <- round(d[,s])
                series[[s]]$tooltip=list(
                    pointFormat=paste('<strong>Population: ',
                        length(box.list[[n]][[s]]),'</strong><br/>',
                        'Maximum: ',r[6],'<br/>',
                        'Upper quartile: ',r[5],'<br/>',
                        'Median: ',r[4],'<br/>',
                        'Lower quartile: ',r[3],'<br/>',
                        'Minimum: ',r[2],'<br/>',sep="")
                )
            }
            
            # Outlier series (if any)
            counter <- 0
            outliers <- vector("list",length(biotypes))
            names(outliers) <- biotypes
            for (o in names(outliers)) {
                counter <- counter + 1
                outliers[[o]] <- list()
                outliers[[o]]$id <- o
                outliers[[o]]$name <- o
                outliers[[o]]$type <- "scatter"
                outliers[[o]]$showInLegend <- FALSE
                #outliers[[o]]$turboThreshold <- 10000
                outliers[[o]]$color <- cols$fill[counter]
                outliers[[o]]$marker <- list(
                    fillColor=cols$fill[counter],
                    symbol="circle",
                    lineWidth=1,
                    lineColor=cols$border[counter]
                )
                outliers[[o]]$data <- list()
                x <- rep(d[1,o],length(oList[[o]]))
                names(x) <- names(oList[[o]])
                if (is.null(obj$altnames)) {
                    outliers[[o]]$data <- 
                        makeHighchartsPoints(x,oList[[o]])
                }
                else {
                    outliers[[o]]$data <- 
                        makeHighchartsPoints(x,oList[[o]],
                            unname(altnames[names(x)]))
                }
            }
            
            json[[n]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="boxplot"
                            ),
                            title=list(
                                text=paste("Biotype detection for sample ",n,
                                    sep="")
                            ),
                            legend=list(
                                enabled=TRUE,
                                itemHoverStyle=list(
                                    color="#B40000"
                                )
                            ),
                            xAxis=list(
                                categories=biotypes,
                                title=list(
                                    text="Biotype",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                type="logarithmic",
                                showFirstLabel=FALSE,
                                min=1e-4,
                                tickInterval=1,
                                title=list(
                                    useHTML=TRUE,
                                    #text="Read count (log<sub>2</sub>)",
                                    text="Expression (read count)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    ),
                                    formatter=y.label.formatter
                                )
                            ),
                            plotOptions=list(
                                boxplot=list(
                                    fillColor="#F0F0E0",
                                    lineWidth=2,
                                    medianColor="#000000",
                                    medianWidth=3,
                                    stemColor="#000000",
                                    stemDashStyle="dash",
                                    stemWidth=1,
                                    whiskerColor="#000000",
                                    whiskerLength="75%",
                                    whiskerWidth=1,
                                    grouping=FALSE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-size:1.1em;',
                                            'color:{series.color};',
                                            'font-weight:bold">',
                                            '\u25CF </span>',
                                            '<span style="font-size:1.1em;',
                                            'font-weight:bold">',
                                            'Biotype {series.name}</span><br/>',
                                            sep=""
                                        )
                                    ),
                                    events=list(
                                        legendItemClick=boxplot.onclick
                                    )
                                ),
                                scatter=list(
                                    allowPointSelect=TRUE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Biotype {series.name}</span><br/>',
                                            sep=""
                                        ),
                                        pointFormat=outlier.pointformat
                                    ),
                                    states=list(
                                        hover=list(
                                            marker=list(
                                                enabled=FALSE
                                            )
                                        )
                                    )
                                )
                            ),
                            series=c(unname(series),unname(outliers))
                        )
                    )
                }
            )
        }
        return(.unquote_js_fun(json))
    }
    else if (by=="biotype") {
        cols <- .getColorScheme(length(samples))
        box.list <- json <- vector("list",length(biotypes))
        names(box.list) <- names(json) <- biotypes
        for (b in biotypes) {
            box.list[[b]] <- vector("list",length(samplenames))
            names(box.list[[b]]) <- samplenames
            for (n in samplenames)
                box.list[[b]][[n]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(box.list[[b]],plot=FALSE)$stats
            colnames(B) <- samplenames
            oList <- lapply(names(box.list[[b]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },box.list[[b]],B)
            names(oList) <- samplenames
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in 1:ncol(B))
                BB[,jj] <- round(B[,jj],3)
            d <- as.data.frame(BB)
            ids <- 0:(ncol(d)-1)
            d <- rbind(ids,d)
            names(ids) <- colnames(d)
            counter <- 0
            series <- vector("list",length(samples))
            names(series) <- names(samples)
            for (s in names(series)) {
                counter <- counter + 1
                series[[s]] <- list()
                series[[s]]$name=s
                if (grouped)
                    series[[s]]$color=cols$fill[counter]
                else
                    series[[s]]$color=cols$fill[1]
                m <- match(samples[[s]],colnames(d))
                series[[s]]$data <- unname(as.list(d[,m]))
            }
            
            # Outlier series (if any)
            counter <- 0
            outliers <- vector("list",length(samples))
            names(outliers) <- names(samples)
            for (o in names(outliers)) {
                counter <- counter + 1
                outliers[[o]] <- list()
                outliers[[o]]$id <- o
                outliers[[o]]$name <- o
                outliers[[o]]$type <- "scatter"
                outliers[[o]]$showInLegend <- FALSE
                if (grouped) {
                    outliers[[o]]$color <- cols$fill[counter]
                    outliers[[o]]$marker <- list(
                        fillColor=cols$fill[counter],
                        symbol="circle",
                        lineWidth=1,
                        lineColor=cols$border[counter]
                    )
                }
                else {
                    outliers[[o]]$color <- cols$fill[1]
                    outliers[[o]]$marker <- list(
                        fillColor=cols$fill[1],
                        symbol="circle",
                        lineWidth=1,
                        lineColor=cols$border[1]
                    )
                }
                outliers[[o]]$data <- list()
                m <- match(samples[[o]],colnames(d))
                if (length(m)>0) {
                    for (i in m) {
                        x <- rep(d[1,i],length(oList[[i]]))
                        names(x) <- names(oList[[i]])
                        if (is.null(obj$altnames)) {
                            outliers[[o]]$data <- 
                                makeHighchartsPoints(x,oList[[i]])
                        }
                        else {
                            outliers[[o]]$data <- c(outliers[[o]]$data,
                                makeHighchartsPoints(x,oList[[i]],
                                unname(altnames)))
                        }
                    }
                }
            }
            
            json[[b]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="boxplot"
                            ),
                            title=list(
                                text=paste("Detection for biotype ",b,
                                    " (population: ",lengths(box.list[[b]])[1],
                                    ")",sep="")
                            ),
                            legend=list(
                                enabled=TRUE
                            ),
                            xAxis=list(
                                categories=samplenames,
                                title=list(
                                    text="Sample name",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                type="logarithmic",
                                showFirstLabel=FALSE,
                                min=1e-4,
                                tickInterval=1,
                                title=list(
                                    text="Expression (read count)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    ),
                                    formatter=y.label.formatter
                                )
                            ),
                            plotOptions=list(
                                boxplot=list(
                                    fillColor="#F0F0E0",
                                    lineWidth=2,
                                    medianColor="#000000",
                                    medianWidth=3,
                                    stemColor="#000000",
                                    stemDashStyle="dash",
                                    stemWidth=1,
                                    whiskerColor="#000000",
                                    whiskerLength="75%",
                                    whiskerWidth=1,
                                    grouping=FALSE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-size:1.1em;',
                                            'color:{series.color};',
                                            'font-weight:bold">',
                                            '\u25CF </span>',
                                            '<span style="font-size:1.1em;',
                                            'font-weight:bold">',
                                            'Condition {series.name}</span>',
                                            '<br/>',
                                            '<span style="font-weight:bold;">',
                                            'Sample {point.key}',
                                            '</span><br/>',sep=""
                                        ),
                                        pointFormatter=tooltip.point.formatter
                                    ),
                                    events=list(
                                        legendItemClick=boxplot.onclick
                                    )
                                ),
                                scatter=list(
                                    allowPointSelect=TRUE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Condition {series.name}</span>',
                                            '<br/>',sep=""
                                        ),
                                        pointFormat=outlier.pointformat
                                    ),
                                    states=list(
                                        hover=list(
                                            marker=list(
                                                enabled=FALSE
                                            )
                                        )
                                    )
                                )
                            ),
                            series=c(unname(series),unname(outliers))
                        )
                    )
                }
            )
        }
        return(.unquote_js_fun(json))
    }
}

bioDetectionToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    samples <- obj$samples
    status <- obj$status
    plotdata <- obj$user$plotdata
    covars <- obj$user$covars
    
    if (!is.null(samples)&& is.list(samples)) {
        samplenames <- unlist(samples,use.names=FALSE)
        names(plotdata$biotables) <- samplenames
    }
    # Otherwise we are using the names present in the input object
    
    abu <- which(plotdata$genome>7)
    nabu <- which(plotdata$genome<=7)
    
    cols <- .getColorScheme()
    json <- vector("list",length(samplenames))
    names(json) <- samplenames
    for (n in samplenames) {
        # Data series
        series.abu <- vector("list",3)
        names(series.abu) <- c("genome","detectionVSgenome","detectionVSsample")
        series.abu$genome <- list()
        series.abu$genome$id <- "abu_genome"
        series.abu$genome$name <- "% in genome"
        series.abu$genome$color <- cols$trans[1]
        series.abu$genome$pointPlacement <- -0.2
        series.abu$genome$data <- round(as.numeric(plotdata$genome[abu]),3)
        series.abu$detectionVSgenome <- list()
        series.abu$detectionVSgenome$id <- "abu_detected"
        series.abu$detectionVSgenome$name <- "% detected"
        series.abu$detectionVSgenome$color <- cols$trans[2]
        series.abu$detectionVSgenome$pointPlacement <- 0
        series.abu$detectionVSgenome$data <- round(as.numeric(
            plotdata$biotables[[n]][1,abu]),3)
        series.abu$detectionVSsample <- list()
        series.abu$detectionVSsample$id <- "abu_sample"
        series.abu$detectionVSsample$name <- "% in sample"
        series.abu$detectionVSsample$color <- cols$trans[3]
        series.abu$detectionVSsample$pointPlacement <- 0.2
        series.abu$detectionVSsample$data <- round(as.numeric(
            plotdata$biotables[[n]][2,abu]),3)
        series.nabu <- vector("list",3)
        names(series.nabu) <- c("genome","detectionVSgenome",
            "detectionVSsample")
        series.nabu$genome <- list()
        series.nabu$genome$name <- "% in genome"
        series.nabu$genome$yAxis <- 1
        series.nabu$genome$pointStart <- length(abu)
        series.nabu$genome$linkedTo <- "abu_genome"
        series.nabu$genome$color <- cols$trans[1]
        series.nabu$genome$pointPlacement <- -0.2
        series.nabu$genome$data <- round(as.numeric(plotdata$genome[nabu]),3)
        series.nabu$detectionVSgenome <- list()
        series.nabu$detectionVSgenome$name <- "% detected"
        series.nabu$detectionVSgenome$yAxis <- 1
        series.nabu$detectionVSgenome$pointStart <- length(abu)
        series.nabu$detectionVSgenome$linkedTo <- "abu_detected"
        series.nabu$detectionVSgenome$color <- cols$trans[2]
        series.nabu$detectionVSgenome$pointPlacement <- 0
        series.nabu$detectionVSgenome$data <- round(as.numeric(
            plotdata$biotables[[n]][1,nabu]),3)
        series.nabu$detectionVSsample <- list()
        series.nabu$detectionVSsample$name <- "% in sample"
        series.nabu$detectionVSsample$yAxis <- 1
        series.nabu$detectionVSsample$pointStart <- length(abu)
        series.nabu$detectionVSsample$linkedTo <- "abu_sample"
        series.nabu$detectionVSsample$color <- cols$trans[3]
        series.nabu$detectionVSsample$pointPlacement <- 0.2
        series.nabu$detectionVSsample$data <- round(as.numeric(
            plotdata$biotables[[n]][2,nabu]),3)
        
        json[[n]] <- switch(jl,
            highcharts = {
                toJSON(
                    list(
                        chart=list(
                            type="column",
                            alignTicks=FALSE
                        ),
                        title=list(
                            text=paste("Comparative biotype detection for ",
                                "sample ",n,sep="")
                        ),
                        legend=list(
                            enabled=TRUE,
                            itemHoverStyle=list(
                                color="#B40000"
                            )
                        ),
                        tooltip=list(
                            shared=TRUE
                        ),  
                        xAxis=list(
                            categories=names(plotdata$genome)[c(abu,nabu)],
                            title=list(
                                text="Biotype",
                                margin=25,
                                style=list(
                                    color="#000000",
                                    fontSize="1.2em"
                                )
                            ),
                            labels=list(
                                style=list(
                                    color="#000000",
                                    fontWeight="bold"
                                )
                            ),
                            plotLines=list(
                                list(
                                    color="#8A8A8A",
                                    width=1.5,
                                    dashStyle="Dash",
                                    value=length(abu)-0.5
                                )
                            ),
                            plotBands=list(
                                list(
                                    color="#FFFFE0",
                                    from=-0.5,
                                    to=length(abu)-0.5
                                ),
                                list(
                                    color="#FFECEB",
                                    from=length(abu)-0.5,
                                    to=length(plotdata$genome)
                                )
                            )
                        ),
                        yAxis=list(
                            list(
                                min=0,
                                max=70,
                                title=list(
                                    text="% of abundant features",
                                    margin=20,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),                          
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            list(
                                min=0,
                                max=7,
                                title=list(
                                    text="% of non-abundant features",
                                    margin=20,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),                          
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                ),
                                opposite=TRUE
                            )
                        ),
                        plotOptions=list(
                            column=list(
                                grouping=FALSE,
                                shadow=FALSE,
                                groupPadding=0.3,
                                pointPadding=0.25,
                                tooltip=list(
                                    headerFormat=paste(
                                        '<span style="font-size:1.1em;',
                                        'font-weight:bold">',
                                        '{point.key}</span><br/>',sep=""
                                    )
                                )
                            )
                        ),
                        series=c(unname(series.abu),unname(series.nabu))
                    )
                )
            }
        )
    }
    return(json)
}

bioSaturationToJSON <- function(obj,by=c("sample","biotype"),
    jl=c("highcharts")) {
    
    by <- tolower(by[1])
    jl <- tolower(jl[1])
    samples <- obj$samples
    plotdata <- obj$user$plotdata
    
    if (!is.null(samples)&& is.list(samples)) {
        samplenames <- unlist(samples,use.names=FALSE)
        names(plotdata) <- samplenames
    }
    
    if (by=="sample") {
        json <- vector("list",length(samplenames))
        names(json) <- samplenames
        for (n in samplenames) {
            depth <- round(plotdata[[n]][,1]/1e+6)
            global <- round(plotdata[[n]][,2])
            M <- plotdata[[n]][,3:ncol(plotdata[[n]])]
            
            # To determine the separation
            ord <- sort(M[nrow(M),],decreasing=TRUE,index.return=TRUE)
            abu <- ord$ix[1:2]
            names(abu) <- names(ord$x[1:2])
            nabu <- ord$ix[3:length(ord$ix)]
            names(nabu) <- names(ord$x[3:length(ord$x)])
    
            cols <- .getColorScheme(ncol(plotdata[[n]])-1)
            counter <- 1
            
            global.series <- list(
                global=list(
                    id="global",
                    name="global",
                    color=cols$fill[counter],
                    data=makeHighchartsPoints(depth,global)
                )
            )
            
            abu.series <- vector("list",2)
            names(abu.series) <- names(abu)
            for (s in names(abu.series)) {
                counter <- counter + 1
                abu.series[[s]] <- list()
                abu.series[[s]]$id <- s
                abu.series[[s]]$name <- s
                abu.series[[s]]$color <- cols$fill[counter]
                abu.series[[s]]$data <- makeHighchartsPoints(depth,
                    round(M[,s]))
            }
            
            nabu.series <- vector("list",length(3:ncol(M)))
            names(nabu.series) <- names(nabu)
            for (s in names(nabu.series)) {
                counter <- counter + 1
                nabu.series[[s]] <- list()
                nabu.series[[s]]$id <- s
                nabu.series[[s]]$name <- s
                nabu.series[[s]]$color <- cols$fill[counter]
                nabu.series[[s]]$data <- makeHighchartsPoints(depth,
                    round(M[,s]))
            }
            
            json[[n]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="scatter",
                                zoomType="xy"
                            ),
                            title=list(
                                text=paste("Biotype saturations for sample ",n,
                                    sep="")
                            ),
                            legend=list(
                                enabled=TRUE,
                                itemHoverStyle=list(
                                    color="#B40000"
                                )
                            ),
                            xAxis=list(
                                title=list(
                                    text="Read depth (millions of reads)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                min=0,
                                max=max(global),
                                title=list(
                                    text="Detected features",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            plotOptions=list(
                                series=list(
                                    lineWidth=2
                                ),
                                scatter=list(
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Biotype {series.name}</span><br/>',
                                            sep=""
                                        ),
                                        pointFormat=paste(
                                            "Depth: {point.x}M<br/>",
                                            "Detected features: {point.y}",
                                            sep="")
                                    )
                                )
                            ),
                            series=c(unname(global.series),
                                unname(abu.series),unname(nabu.series))
                        )
                    )
                }
            )
        }
        return(json)
    }
    else if (by=="biotype") {
        biotypes <- colnames(plotdata[[1]])[2:ncol(plotdata[[1]])]
        depths <- vector("list",length(plotdata))
        names(depths) <- samplenames
        for (n in samplenames)
            depths[[n]] <- round(plotdata[[n]][,1]/1e+6)
        json <- vector("list",length(biotypes))
        names(json) <- biotypes
        
        for (b in biotypes) {
            series <- vector("list",length(plotdata))
            names(series) <- samplenames
            cols <- .getColorScheme(length(samplenames))
            counter <- 0
            for (s in names(series)) {
                counter <- counter + 1
                series[[s]] <- list()
                series[[s]]$id <- s
                series[[s]]$name <- s
                series[[s]]$color <- cols$fill[counter]
                series[[s]]$data <- makeHighchartsPoints(depths[[s]],
                    round(plotdata[[s]][,b]))
            }
            
            json[[b]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="scatter",
                                zoomType="xy"
                            ),
                            title=list(
                                text=paste("Sample saturations for biotype ",b,
                                    sep="")
                            ),
                            legend=list(
                                enabled=TRUE,
                                itemHoverStyle=list(
                                    color="#B40000"
                                )
                            ),
                            xAxis=list(
                                title=list(
                                    text="Read depth (millions of reads)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                title=list(
                                    text="Detected features",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            plotOptions=list(
                                series=list(
                                    lineWidth=2
                                ),
                                scatter=list(
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Sample {series.name}</span><br/>',
                                            sep=""
                                        ),
                                        pointFormat=paste(
                                            "Depth: {point.x}M<br/>",
                                            "Detected features: {point.y}",
                                            sep="")
                                    )
                                )
                            ),
                            series=c(unname(series))
                        )
                    )
                }
            )
        }
        return(json)
    }
}

readNoiseToJSON <- function(obj,jl=c("highcharts"),seed=42) {
    jl <- tolower(jl[1])
    d <- obj$user
    samples <- obj$samples
    
    # Too many points for a lot of curves of interactive data
    if (nrow(d)>1000) {
        set.seed(seed)
        ii <- sort(sample(1:nrow(d),998))
        ii <- c(1,ii,nrow(d))
        d <- cbind(d[ii,1],d[ii,2:ncol(d)])
    }

    if (is.null(samples)) 
        samples <- 1:(ncol(d)-1)
    if (is.numeric(samples)) 
        samplenames = colnames(dat)[samples+1]
    if (is.list(samples))
        samplenames <- unlist(samples)
    
    cols <- .getColorScheme(length(samplenames))
    
    # Construct series
    counter <- 0
    series <- vector("list",length(samplenames))
    names(series) <- samplenames
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        series[[n]]$color=cols$fill[counter]
        series[[n]]$data <- makeHighchartsPoints(d[,1],d[,n])
        series[[n]]$tooltip=list(
            headerFormat=paste("<span style=",
                "\"font-size:1.1em;color:{series.color};",
                "font-weight:bold\">{series.name}<br>",
                sep=""),
            pointFormat=NULL
        )
    }
    
    switch(jl,
        highcharts = {
                json <- toJSON(list(
                    chart=list(
                        type="line",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("RNA-Seq mapped reads noise")
                    ),
                    xAxis=list(
                        title=list(
                            text="% detected features",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=0,
                        max=100
                    ),
                    yAxis=list(
                        title=list(
                            text="% of total reads",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        tickPositions=seq(0,110,10)
                    ),
                    plotOptions=list(
                        line=list(
                            allowPointSelect=TRUE,
                            lineWidth=1,
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=NULL
                            ),
                            turboThreshold=50000
                        )
                    ),
                    series=unname(series)
                )
            )
        }
    )
    return(json)
}

boxplotToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    b <- obj$plot
    name <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    oList <- obj$user
    
    grouped <- FALSE
    if (is.null(name)) {
        if (is.null(colnames(b$stats)))
            nams <- paste("Sample",1:ncol(b$stats),sep=" ")
        else
            nams <- colnames(b$stats)
        name <- list(Samples=nams)
    }
    else if (length(name)==1 && name=="none") {
        nams <- rep("",ncol(b$stats))
        name <- list(Samples=nams)
    }
    else if (is.list(name)) { # Is sampleList
        nams <- unlist(name,use.names=FALSE)
        grouped <- TRUE
    }
    cols <- .getColorScheme()
    
    # Data series
    d <- as.data.frame(round(b$stat,3))
    ids <- 0:(ncol(d)-1)
    d <- rbind(ids,d)
    colnames(d) <- nams
    names(ids) <- colnames(d)
    counter <- 0
    series <- vector("list",length(name))
    names(series) <- names(name)
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        if (grouped)
            series[[n]]$color=cols$fill[counter]
        else
            series[[n]]$color=cols$fill[1]
        m <- match(name[[n]],colnames(d))
        series[[n]]$data <- unname(as.list(d[,m]))
    }
    # Outlier series (if any)
    counter <- 0
    outliers <- vector("list",length(name))
    names(outliers) <- names(name)
    for (n in names(outliers)) {
        counter <- counter + 1
        outliers[[n]] <- list()
        outliers[[n]]$id <- n
        outliers[[n]]$name <- n
        outliers[[n]]$type <- "scatter"
        outliers[[n]]$showInLegend <- FALSE
        if (grouped) {
            outliers[[n]]$color <- cols$fill[counter]
            outliers[[n]]$marker <- list(
                fillColor=cols$fill[counter],
                symbol="circle",
                lineWidth=1,
                lineColor=cols$border[counter]
            )
        }
        else {
            outliers[[n]]$color <- cols$fill[1]
            outliers[[n]]$marker <- list(
                fillColor=cols$fill[1],
                symbol="circle",
                lineWidth=1,
                lineColor=cols$border[1]
            )
        }
        outliers[[n]]$data <- list()
        m <- match(name[[n]],colnames(d))
        if (length(m)>0) {
            for (i in m)
                x <- rep(d[1,i],length(oList[[i]]))
                names(x) <- names(oList[[i]])
                outliers[[n]]$data <- c(outliers[[n]]$data,
                    makeHighchartsPoints(x,oList[[i]],unname(altnames)))
        }
    }
        
    # Boxplot tooltip point formatter for the case of zeros
    tooltip.point.formatter <- paste("function() {",
        "   var min = this.low === 0.001 ? 0 : this.low;" ,
        "   var q1 = this.q1 === 0.001 ? 0 : this.q1;" ,
        "   var med = this.median === 0.001 ? 0 : this.median;",
        "   var q3 = this.q3 === 0.001 ? 0 : this.q3;",
        "   var max = this.high === 0.001 ? 0 : this.high;",
        "   var str = 'Maximum: ' + max + '<br/>' +",
        "       'Upper quartile: ' + q3 + '<br/>' +",
        "       'Median: ' + med + '<br/>' +",
        "       'Lower quartile: ' + q1 + '<br/>' +",
        "       'Minimum: ' + min + '<br/>';",
        "   return  str;",
        "}",sep="")
    
    # Legend clicker
    boxplot.onclick <- paste("function() {",
        "   var chart =  this.chart;",
        "   var outlier_id =  chart.get(this.name);",
        "   if (!outlier_id.visible) {",
        "       outlier_id.show();",
        "   } else {",
        "       outlier_id.hide();",
        "   }",
        "}",sep="")
    
    if (is.null(obj$altnames)) {
        outlier.pointformat=paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    else {
        outlier.pointformat=paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Gene name: {point.alt_name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    
    json <- switch(jl,
        highcharts = {
            toJSON(
                list(
                    chart=list(
                        type="boxplot"
                    ),
                    title=list(
                        text=paste("Boxplot ",status,sep="")
                    ),
                    legend=list(
                        enabled=TRUE
                    ),
                    xAxis=list(
                        categories=nams,
                        title=list(
                            text="Sample name",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontWeight="bold"
                            )
                        )
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="Read count (log<sub>2</sub>)",
                            margin=30,
                            style=list(
                                color="#000000",
                                fontSize="1.1em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        )
                    ),
                    plotOptions=list(
                        boxplot=list(
                            fillColor="#F0F0E0",
                            lineWidth=2,
                            medianColor="#000000",
                            medianWidth=3,
                            stemColor="#000000",
                            stemDashStyle="dash",
                            stemWidth=1,
                            whiskerColor="#000000",
                            whiskerLength="75%",
                            whiskerWidth=1,
                            grouping=FALSE,
                            tooltip=list(
                                headerFormat=paste(
                                    '<span style="font-size:1.1em;',
                                    'color:{series.color};',
                                    'font-weight:bold">',
                                    '\u25CF </span>',
                                    '<span style="font-size:1.1em;',
                                    'font-weight:bold">',
                                    'Condition {series.name}</span><br/>',
                                    '<span style="font-weight:bold">',
                                    'Sample {point.key}</span><br/>',sep=""
                                ),
                                pointFormatter=tooltip.point.formatter
                            ),
                            events=list(
                                legendItemClick=boxplot.onclick
                            )
                        ),
                        scatter=list(
                            allowPointSelect=TRUE,
                            tooltip=list(
                                headerFormat=paste(
                                    '<span style="font-weight:bold;',
                                    'color:{series.color};">',
                                    '\u25CF </span>',
                                    '<span style="font-weight:bold">',
                                    'Condition {series.name}</span><br/>',
                                    sep=""
                                ),
                                pointFormat=outlier.pointformat
                            ),
                            states=list(
                                hover=list(
                                    marker=list(
                                        enabled=FALSE
                                    )
                                )
                            )
                            
                        )
                    ),
                    series=c(unname(series),unname(outliers))
                )
            )
        }
    )
    return(.unquote_js_fun(json))
}

biasPlotToJSON <- function(obj,jl=c("highcharts"),seed=1) {
    jl <- tolower(jl[1])
    counts <- round(nat2log(obj$user$counts),3)
    status <- obj$status
    covar <- obj$user$covar
    covarname <- obj$user$covarname
    samples <- obj$samples
    
    # Too many points for a lot of curves of interactive data
    if (nrow(counts)>2000) {
        set.seed(seed)
        ii <- sample(1:nrow(counts),2000)
        counts <- counts[ii,]
        covar <- covar[ii]
    }
    
    # If length bias, not nice to have x-axis at -200k
    minX <- ifelse(max(covar>100),0,"undefined")

    if (is.null(samples)) {
        if (is.null(colnames(x)))
            samplenames <- paste("Sample",1:ncol(counts),sep=" ")
        else
            samplenames <- colnames(counts)
        samples <- list(Samples=nams)
        cols <- .getColorScheme(length(samples))
    }
    else if (is.list(samples)) { # Is sampleList
        samplenames <- unlist(samples,use.names=FALSE)
        grouped <- TRUE
    }
    colnames(counts) <- samplenames
    
    # Construct series
    counter <- 0
    series <- vector("list",length(samplenames))
    names(series) <- samplenames
    for (n in names(series)) {
        counter <- counter + 1
        x <- counts[,n]
        fit <- lowess(covar,x)
        series[[n]] <- list()
        series[[n]]$name <- n
        series[[n]]$color <- cols$fill[counter]
        series[[n]]$data <- lapply(1:length(x),function(i,x,y) {
            return(c(x[i],y[i])) },round(fit$x,3),round(fit$y,3))
    }
    
    switch(jl,
        highcharts = {
                json <- toJSON(list(
                    chart=list(
                        type="line",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste(covarname," bias detection - ",status)
                    ),
                    xAxis=list(
                        min=minX,
                        title=list(
                            text=covarname,
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="Read count (log<sub>2</sub>)",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE
                    ),
                    plotOptions=list(
                        line=list(
                            marker=list(
                                enabled=FALSE,
                                states=list(
                                    hover=list(
                                        enabled=FALSE
                                    )
                                )
                            ),
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=NULL
                            ),
                            turboThreshold=50000
                        )
                    ),
                    series=unname(series)
                )
            )
        }
    )
    return(json)
}

filteredToJSON <- function(obj,by=c("chromosome","biotype"),
    jl=c("highcharts")) {
    
    jl <- tolower(jl[1])
    by <- tolower(by[1])
    filtered <- obj$user$filtered
    total <- obj$user$total
    cols <- .getColorScheme(2)
    
    if (by=="chromosome") {
        chr <- table(as.character(filtered$chromosome))
        chr.all <- table(as.character(total$chromosome))
        barlab.chr <- as.character(chr)        
        per.chr <- round(chr/chr.all[names(chr)],3)
        per.chr[per.chr>1] <- 1
        
        series <- vector("list",2)
        names(series) <- c("number","fraction")
                
        series$number <- list()
        series$number$id <- "chr_number"
        series$number$name <- "Number of genes"
        series$number$color <- cols$fill[1]
        series$number$pointPlacement <- -0.2
        series$number$data <- unname(chr)
        
        series$fraction <- list()
        series$fraction$id <- "chr_fraction"
        series$fraction$name <- "Fraction of total genes"
        series$fraction$color <- cols$fill[2]
        series$fraction$pointPlacement <- 0.2
        series$fraction$yAxis <- 1
        series$fraction$data <- unname(per.chr)
        
        what <- chr
    }
    else if (by=="biotype") {
        bt <- table(as.character(filtered$biotype))
        bt.all <- table(as.character(total$biotype))
        barlabBt <- as.character(bt)
        per.bt <- round(bt/bt.all[names(bt)],3)
        per.bt[per.bt>1] <- 1
        
        series <- vector("list",2)
        names(series) <- c("number","fraction")
                
        series$number <- list()
        series$number$id <- "bt_number"
        series$number$name <- "Number of genes"
        series$number$color <- cols$fill[1]
        series$number$pointPlacement <- -0.2
        series$number$data <- unname(bt)
        
        series$fraction <- list()
        series$fraction$id <- "bt_fraction"
        series$fraction$name <- "Fraction of total genes"
        series$fraction$color <- cols$fill[2]
        series$fraction$pointPlacement <- 0.2
        series$fraction$yAxis <- 1
        series$fraction$data <- unname(per.bt)
        
        what <- bt
    }
    
    json <- switch(jl,
        highcharts = {
            toJSON(
                list(
                    chart=list(
                        type="column",
                        alignTicks=FALSE
                    ),
                    title=list(
                        text=paste("Filtered genes per ",by,sep="")
                    ),
                    legend=list(
                        enabled=TRUE,
                        itemHoverStyle=list(
                            color="#B40000"
                        )
                    ),
                    tooltip=list(
                        shared=TRUE
                    ),
                    xAxis=list(
                        categories=names(what),
                        title=list(
                            text=paste(toupper(substr(by,1,1)),substr(by,2,
                                nchar(by)),sep=""),
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontWeight="bold"
                            )
                        )
                    ),
                    yAxis=list(
                        list(
                            lineColor=cols$fill[1],
                            lineWidth=2,
                            min=0,
                            tickAmount=11,
                            title=list(
                                text="Number of genes",
                                margin=20,
                                style=list(
                                    color="#000000",
                                    fontSize="1.2em"
                                )
                            ),                 
                            labels=list(
                                style=list(
                                    color="#000000",
                                    fontSize="1.1em",
                                    fontWeight="bold"
                                )
                            ),
                            offset=10
                        ),
                        list(
                            lineColor=cols$fill[2],
                            lineWidth=2,
                            min=0,
                            max=1,
                            tickAmount=11,
                            #tickInterval=0.1,
                            title=list(
                                text="Fraction of total genes",
                                margin=20,
                                style=list(
                                    color="#000000",
                                    fontSize="1.2em"
                                )
                            ),                          
                            labels=list(
                                style=list(
                                    color="#000000",
                                    fontSize="1.1em",
                                    fontWeight="bold"
                                )
                            ),
                            opposite=TRUE,
                            offset=10
                        )
                    ),
                    plotOptions=list(
                        column=list(
                            grouping=FALSE,
                            shadow=FALSE,
                            groupPadding=0.2,
                            pointPadding=0.2,
                            tooltip=list(
                                headerFormat=paste(
                                    '<span style="font-size:1.1em;',
                                    'font-weight:bold">',
                                    '{point.key}</span><br/>',sep=""
                                )
                            )
                        )
                    ),
                    series=c(unname(series))
                )
            )
        }
    )
    
    return(json)
}

volcanoToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    f <- obj$x
    p <- obj$y
    xlim <- obj$xlim
    ylim <- obj$ylim
    altNames <- obj$altnames
    pcut <- obj$pcut
    fcut <- obj$fcut
    up <- obj$user$up
    down <- obj$user$down
    ff <- obj$user$unf
    pp <- obj$user$unp
    altNamesNeutrall <- obj$user$ualt
    con <- obj$user$con
    
    switch(jl,
        highcharts = {
            if (is.null(altNames))
                point.format=paste("<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change: </strong>{point.x}<br>",
                    "<strong>Significance: </strong>{point.y}",sep="")
            else
                point.format=paste("<strong>Gene name: </strong>",
                    "{point.alt_name}<br>",
                    "<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change: </strong>{point.x}<br>",
                    "<strong>Significance: </strong>{point.y}",sep="")
                json <- toJSON(
                    list(
                        chart=list(
                        type="scatter",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("Volcano plot for",con)
                    ),
                    xAxis=list(
                        title=list(
                            text="Fold change",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(xlim[1],3),
                        max=round(xlim[2],3)
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="Significance (-log<sub>10</sub>(p-value))",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(ylim[1]-2,3),
                        max=round(ylim[2],3)
                    ),
                    #legend=list(
                    #    layout="vertical",
                    #    align="left",
                    #    verticalAlign="top",
                    #    floating=TRUE,
                    #    backgroundColor="#FFFFFF",
                    #    borderWidth=1
                    #),
                    plotOptions=list(
                        scatter=list(
                            allowPointSelect=TRUE,
                            marker=list(
                                radius=2,
                                states=list(
                                    hover=list(
                                        enabled=TRUE,
                                        lineColor="#333333"
                                    )
                                )
                            ),
                            states=list(
                                hover=list(
                                    marker=list(
                                        enabled=FALSE
                                    )
                                )
                            ),
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=point.format
                            ),
                            turboThreshold=50000
                        )
                    ),
                    series=list(
                        list(
                            name="Up-regulated",
                            color="#EE0000",
                            marker=list(
                                symbol="circle"
                            ),
                            data=makeHighchartsPoints(f[up],-log10(p[up]),
                                unname(altNames[up]))
                        ),
                        list(
                            name="Down-regulated",
                            marker=list(
                                symbol="circle"
                            ),
                            color="#00CD00",
                            data=makeHighchartsPoints(f[down],-log10(p[down]),
                                unname(altNames[down]))
                        ),
                        list(
                            name="Unregulated",
                            marker=list(
                                symbol="circle"
                            ),
                            color="#0000EE",
                            data=makeHighchartsPoints(ff,-log10(pp),
                                unname(altNamesNeutrall))
                        ),
                        list(
                            name="Downfold threshold",
                            color="#000000",
                            type="line",
                            dashStyle="dash",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.x}<br/>',sep="")
                            ),
                            data=list(round(c(-fcut,ylim[1]-5),3),
                                round(c(-fcut,ylim[2]),3))
                        ),
                        list(
                            name="Upfold threshold",
                            color="#000000",
                            type="line",
                            dashStyle="Dash",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.x}<br/>',sep="")
                            ),
                            data=list(round(c(fcut,ylim[1]-5),3),
                                round(c(fcut,ylim[2]),3))
                        ),
                        list(
                            name="Significance threshold",
                            color="#000000",
                            type="line",
                            dashStyle="DashDot",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.y}<br/>',sep="")
                            ),
                            data=list(round(c(xlim[1],-log10(pcut)),3),
                                round(c(xlim[2],-log10(pcut)),3))
                        )
                    )
                )
            )
        }
    )
    return(json)
}

.unquote_js_fun <- function(js) {
    if (is.list(js))
        js <- lapply(js,.unquote_js_fun)
    else {
        op <- gregexpr(pattern="function",js)
        cl <- gregexpr(pattern="}\\\"",js)
        if (length(op)>0) {
            starts <- as.numeric(op[[1]])
            for (i in 1:length(starts))
                substr(js,starts[i]-1,starts[i]-1) <- " "
            ends <- as.numeric(cl[[1]])
            for (i in 1:length(starts))
                substr(js,ends[i]+1,ends[i]+1) <- " "
        }
    }
    return(js)
}

.getGroupColorScheme <- function(group) {
    cols <- .getColorScheme(length(group))
    classes <- as.factor(asClassVector(group))
    design <- as.numeric(classes)
    return(lapply(cols,function(x,classes,design) {
        return(x[1:length(levels(classes))][design])
    },classes,design))
}

.getColorScheme <- function(n=NULL) {
    if (missing(n) || is.null(n))
        return(.getColors())
    else {
        cols <- .getColors()
        if (n > length(cols$fill)) {
            cols$fill <- rep(cols$fill,length.out=n)
            cols$border <- rep(cols$border,length.out=n)
            cols$select <- rep(cols$select,length.out=n)
            cols$trans <- rep(cols$trans,length.out=n)
        }
        return(cols)
    }
}

.getColors <- function() {
    return(list(
        fill=c("#CD0000","#00CD00","#0000EE","#FFD700","#87CEEB","#CD8500",
            "#DEB887","#FF0000","#0000FF","#00FF00","#FFA500","#A9A9A9",
            "#008B00","#313131","#FFC0CB","#A52A2A","#FF00FF","#9ACD32",
            "#8B636C","#2E8B57","#008B8B"),
        border=c("#850000","#006B00","#000085","#927C00","#156280","#5A3A00",
            "#8B7457","#935E18","#000080","#008500","#603E00","#454545",
            "#073E07","#000000","#896067","#691111","#7C007C","#3A4D14",
            "#5B1726","#0C2517","#062A2A"),
        selected=c("#FF0000","#00FF00","#0066FF","#FFD700","#FFEB77","#FFB428",
            "#FFD9A5","#FF326D","#0089FF","#B3FF00","#FFC352","#D9D9D9",
            "#00EC00","#8E8E8E","#FFDAE0","#F94444","#FF87FF","#C2FF45",
            "#EA889D","#4EE590","#00DADA"),
        trans=c("rgba(205,0,0,0.6)","rgba(0,205,0,0.6)","rgba(0,0,238,0.6)",
            "rgba(255,215,0,0.6)","rgba(135,206,235,0.6)","rgba(205,133,0,0.6)",
            "rgba(222,184,135,0.6)","rbga(255,0,0,0.5)","rgba(0,0,255,0.5)",
            "rgba(0,255,0,0.5)","rgba(255,165,0,0.6)","rgba(169,169,169,0.5)",
            "rgba(0,139,0,0.6)","rgba(49,49,49,0.6)","rgba(255,192,203,0.5)",
            "rgba(165,42,42,0.6)","rgba(255,0,255,0.6)","rgba(154,205,50,0.6)",
            "rgba(139,99,108,0.6)","rgba(46,139,87,0.6)","rgba(0,139,139,0.6)")
    ))
}
