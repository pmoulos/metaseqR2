mdsToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- out[1]
    
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
                
                json <- list(
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
                )
        }
    )
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else
        return(json)
}

countsBioToJSON <- function(obj,by=c("sample","biotype"),jl=c("highcharts"),
    out=c("json","list")) {
    by <- tolower(by[1])
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    samples <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    counts <- round(2^obj$user$counts - 1)
    #counts[counts==0] <- 0.001
    counts <- counts + 1 # Offset for log
    #counts <- obj$user$counts
    
    covars <- obj$user$covars
    biotypes <- unique(as.character(covars$biotype))
    if (!is.null(altnames))
        names(altnames) <- rownames(counts)
    
    grouped <- FALSE
    if (is.null(samples)) {
        if (is.null(colnames(counts)))
            samplenames <- paste("Sample",seq_len(ncol(counts)),sep=" ")
        else
            samplenames <- colnames(counts)
        samples <- list(Samples=samplenames)
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
        boxList <- json <- vector("list",length(samplenames))
        names(boxList) <- names(json) <- samplenames
        for (n in samplenames) {
            boxList[[n]] <- vector("list",length(biotypes))
            names(boxList[[n]]) <- biotypes
            for (b in biotypes)
                boxList[[n]][[b]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(boxList[[n]],plot=FALSE)$stats
            colnames(B) <- biotypes
            oList <- lapply(names(boxList[[n]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },boxList[[n]],B)
            names(oList) <- biotypes
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in seq_len(ncol(B)))
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
                series[[s]]$data <- list(unname(as.list(d[,s])))
                r <- round(d[,s])
                series[[s]]$tooltip=list(
                    pointFormat=paste('<strong>Population: ',
                        length(boxList[[n]][[s]]),'</strong><br/>',
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
                                turboThreshold=10000,
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
                                turboThreshold=10000,
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
                }
            )
        }
        if (out=="json") {
            for (i in seq_along(json))
                #json[[i]] <- .unquote_js_fun(toJSON(json[[i]],
                #   auto_unbox=TRUE,null="null"))
                json[[i]] <- toJSON(json[[i]],auto_unbox=TRUE,null="null")
            return(json)
        }
        else
            return(json)
    }
    else if (by=="biotype") {
        cols <- .getColorScheme(length(samples))
        boxList <- json <- vector("list",length(biotypes))
        names(boxList) <- names(json) <- biotypes
        for (b in biotypes) {
            boxList[[b]] <- vector("list",length(samplenames))
            names(boxList[[b]]) <- samplenames
            for (n in samplenames)
                boxList[[b]][[n]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(boxList[[b]],plot=FALSE)$stats
            colnames(B) <- samplenames
            oList <- lapply(names(boxList[[b]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },boxList[[b]],B)
            names(oList) <- samplenames
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in seq_len(ncol(B)))
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
                if (length(m) > 1)
                    series[[s]]$data <- unname(as.list(d[,m]))
                else
                    series[[s]]$data <- list(unname(as.list(d[,m])))
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
                    list(
                        chart=list(
                            type="boxplot"
                        ),
                        title=list(
                            text=paste("Detection for biotype ",b,
                                " (population: ",lengths(boxList[[b]])[1],
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
                                turboThreshold=10000,
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
                                turboThreshold=10000,
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
                }
            )
        }
        if (out=="json") {
            for (i in seq_along(json))
                #json[[i]] <- .unquote_js_fun(toJSON(json[[i]],
                #   auto_unbox=TRUE,null="null"))
                json[[i]] <- toJSON(json[[i]],auto_unbox=TRUE,null="null")
            return(json)
        }
        else
            return(json)
    }
}

bioDetectionToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    samples <- obj$samples
    status <- obj$status
    plotdata <- obj$user$plotdata
    covars <- obj$user$covars
    
    if (!is.null(samples) && is.list(samples)) {
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
        seriesAbu <- vector("list",3)
        names(seriesAbu) <- c("genome","detectionVSgenome","detectionVSsample")
        seriesAbu$genome <- list()
        seriesAbu$genome$id <- "abu_genome"
        seriesAbu$genome$name <- "% in genome"
        seriesAbu$genome$color <- cols$trans[1]
        seriesAbu$genome$pointPlacement <- -0.2
        seriesAbu$genome$data <- round(as.numeric(plotdata$genome[abu]),3)
        seriesAbu$detectionVSgenome <- list()
        seriesAbu$detectionVSgenome$id <- "abu_detected"
        seriesAbu$detectionVSgenome$name <- "% detected"
        seriesAbu$detectionVSgenome$color <- cols$trans[2]
        seriesAbu$detectionVSgenome$pointPlacement <- 0
        seriesAbu$detectionVSgenome$data <- round(as.numeric(
            plotdata$biotables[[n]][1,abu]),3)
        seriesAbu$detectionVSsample <- list()
        seriesAbu$detectionVSsample$id <- "abu_sample"
        seriesAbu$detectionVSsample$name <- "% in sample"
        seriesAbu$detectionVSsample$color <- cols$trans[3]
        seriesAbu$detectionVSsample$pointPlacement <- 0.2
        seriesAbu$detectionVSsample$data <- round(as.numeric(
            plotdata$biotables[[n]][2,abu]),3)
        seriesNabu <- vector("list",3)
        names(seriesNabu) <- c("genome","detectionVSgenome",
            "detectionVSsample")
        seriesNabu$genome <- list()
        seriesNabu$genome$name <- "% in genome"
        seriesNabu$genome$yAxis <- 1
        seriesNabu$genome$pointStart <- length(abu)
        seriesNabu$genome$linkedTo <- "abu_genome"
        seriesNabu$genome$color <- cols$trans[1]
        seriesNabu$genome$pointPlacement <- -0.2
        seriesNabu$genome$data <- round(as.numeric(plotdata$genome[nabu]),3)
        seriesNabu$detectionVSgenome <- list()
        seriesNabu$detectionVSgenome$name <- "% detected"
        seriesNabu$detectionVSgenome$yAxis <- 1
        seriesNabu$detectionVSgenome$pointStart <- length(abu)
        seriesNabu$detectionVSgenome$linkedTo <- "abu_detected"
        seriesNabu$detectionVSgenome$color <- cols$trans[2]
        seriesNabu$detectionVSgenome$pointPlacement <- 0
        seriesNabu$detectionVSgenome$data <- round(as.numeric(
            plotdata$biotables[[n]][1,nabu]),3)
        seriesNabu$detectionVSsample <- list()
        seriesNabu$detectionVSsample$name <- "% in sample"
        seriesNabu$detectionVSsample$yAxis <- 1
        seriesNabu$detectionVSsample$pointStart <- length(abu)
        seriesNabu$detectionVSsample$linkedTo <- "abu_sample"
        seriesNabu$detectionVSsample$color <- cols$trans[3]
        seriesNabu$detectionVSsample$pointPlacement <- 0.2
        seriesNabu$detectionVSsample$data <- round(as.numeric(
            plotdata$biotables[[n]][2,nabu]),3)
        
        json[[n]] <- switch(jl,
            highcharts = {
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
                    series=c(unname(seriesAbu),unname(seriesNabu))
                )
            }
        )
    }
    
    if (out=="json") {
        for (i in seq_along(json))
            json[[i]] <- toJSON(json[[i]],auto_unbox=TRUE,null="null")
        return(json)
    }
    else if (out=="list")
        return(json)
}

bioSaturationToJSON <- function(obj,by=c("sample","biotype"),
    jl=c("highcharts"),out=c("json","list")) {
    by <- tolower(by[1])
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
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
            M <- plotdata[[n]][,3:ncol(plotdata[[n]]),drop=FALSE]
            
            # To determine the separation
            ord <- sort(M[nrow(M),],decreasing=TRUE,index.return=TRUE)
            abu <- ord$ix[c(1,2)]
            names(abu) <- names(ord$x[c(1,2)])
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
            
            abuSeries <- vector("list",2)
            names(abuSeries) <- names(abu)
            for (s in names(abuSeries)) {
                counter <- counter + 1
                abuSeries[[s]] <- list()
                abuSeries[[s]]$id <- s
                abuSeries[[s]]$name <- s
                abuSeries[[s]]$color <- cols$fill[counter]
                abuSeries[[s]]$data <- makeHighchartsPoints(depth,
                    round(M[,s]))
            }
            
            nabuSeries <- vector("list",length(3:ncol(M)))
            names(nabuSeries) <- names(nabu)
            for (s in names(nabuSeries)) {
                counter <- counter + 1
                nabuSeries[[s]] <- list()
                nabuSeries[[s]]$id <- s
                nabuSeries[[s]]$name <- s
                nabuSeries[[s]]$color <- cols$fill[counter]
                nabuSeries[[s]]$data <- makeHighchartsPoints(depth,
                    round(M[,s]))
            }
            
            json[[n]] <- switch(jl,
                highcharts = {
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
                            unname(abuSeries),unname(nabuSeries))
                    )
                }
            )
        }
        
        if (out=="json") {
            for (i in seq_along(json))
                json[[i]] <- toJSON(json[[i]],auto_unbox=TRUE,null="null")
            return(json)
        }
        else if (out=="list")
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
                }
            )
        }
        if (out=="json") {
            for (i in seq_along(json))
                json[[i]] <- toJSON(json[[i]],auto_unbox=TRUE,null="null")
            return(json)
        }
        else if (out=="list") {
            names(json) <- biotypes
            return(json)
        }
    }
}

readNoiseToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    d <- obj$user
    samples <- obj$samples
    
    # Too many points for a lot of curves of interactive data
    if (nrow(d)>1000) {
        ii <- sort(sample(seq_len(nrow(d)),998))
        ii <- c(1,ii,nrow(d))
        d <- cbind(d[ii,1],d[ii,2:ncol(d)])
    }

    if (is.null(samples)) 
        samples <- seq_len((ncol(d)-1))
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
        series[[n]]$data <- makeHighchartsPoints(d[,1],d[,n],simple=TRUE)
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
            json <- list(
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
        }
    )
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

boxplotToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    b <- obj$plot
    name <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    oList <- obj$user
    
    grouped <- FALSE
    if (is.null(name)) {
        if (is.null(colnames(b$stats)))
            nams <- paste("Sample",seq_len(ncol(b$stats)),sep=" ")
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
        if (length(m) > 1)
            series[[n]]$data <- unname(as.list(d[,m]))
        else
            series[[n]]$data <- list(unname(as.list(d[,m])))
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
            for (i in m) {
                x <- rep(d[1,i],length(oList[[i]]))
                names(x) <- names(oList[[i]])
                outliers[[n]]$data <- c(outliers[[n]]$data,
                    makeHighchartsPoints(x,oList[[i]],unname(altnames)))
            }
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
        }
    )
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else
        return(json)
}

biasPlotToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    counts <- round(nat2log(obj$user$counts),3)
    status <- obj$status
    covar <- obj$user$covar
    covarname <- obj$user$covarname
    samples <- obj$samples
    
    # Too many points for a lot of curves of interactive data
    if (nrow(counts)>2000) {
        #set.seed(seed)
        ii <- sample(seq_len(nrow(counts)),2000)
        counts <- counts[ii,]
        covar <- covar[ii]
    }
    
    # If length bias, not nice to have x-axis at -200k
    if (covarname == "lengthbias")
        minX <- ifelse(max(covar>100),0,"undefined")
    else
        minX <- 0.1

    if (is.null(samples)) {
        if (is.null(colnames(x)))
            samplenames <- paste("Sample",seq_len(ncol(counts)),sep=" ")
        else
            samplenames <- colnames(counts)
        samples <- list(Samples=samplenames)
        cols <- .getColorScheme(length(samples))
    }
    else if (is.list(samples)) { # Is sampleList
        samplenames <- unlist(samples,use.names=FALSE)
        grouped <- TRUE
        cols <- .getColorScheme(length(samplenames))
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
        series[[n]]$data <- lapply(seq_len(length(x)),function(i,x,y) {
            return(c(x[i],y[i])) },round(fit$x,3),round(fit$y,3))
    }
    
    switch(jl,
        highcharts = {
            json <- list(
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
        }
    )
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

filteredToJSON <- function(obj,by=c("chromosome","biotype"),
    jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    by <- tolower(by[1])
    out <- tolower(out[1])
    
    filtered <- obj$user$filtered
    total <- obj$user$total
    cols <- .getColorScheme(2)
    
    if (by=="chromosome") {
        chrTab <- table(as.character(filtered$chromosome))
        chr <- as.numeric(chrTab)
        names(chr) <- names(chrTab)
        chrAllTab <- table(as.character(total$chromosome))
        chrAll <- as.numeric(chrAllTab)
        names(chrAll) <- names(chrAllTab)
        barlabChr <- as.character(chr)        
        perChr <- round(chr/chrAll[names(chr)],3)
        perChr[perChr>1] <- 1
        
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
        series$fraction$data <- unname(perChr)
        
        what <- chr
    }
    else if (by=="biotype") {
        btTab <- table(as.character(filtered$biotype))
        bt <- as.numeric(btTab)
        names(bt) <- names(btTab)
        btAllTab <- table(as.character(total$biotype))
        btAll <- as.numeric(table(as.character(total$biotype)))
        names(btAll) <- names(btAllTab)
        barlabBt <- as.character(bt)
        perBt <- round(bt/btAll[names(bt)],3)
        perBt[perBt>1] <- 1
        
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
        series$fraction$data <- unname(perBt)
        
        what <- bt
    }
    
    json <- switch(jl,
        highcharts = {
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
        }
    )
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

volcanoToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
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
                json <- list(
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
        }
    )
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

scatterToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    x <- obj$x
    y <- obj$y
    altNames <- obj$altnames
    type <- obj$user$covarname
    status <- obj$status
    samples <- obj$samples
    
    if (type=="Mean-Difference") {
        xLab <- "Mean"
        yLab <- "Difference"
    }
    else if (type=="Mean-Variance") {
        xLab <- "Mean"
        yLab <- "Variance"
    }
    else if (type=="X-Y") {
        xLab <- paste("Expression for",samples[1])
        yLab <- paste("Expression for",samples[2])
    }
    
    # Too many points for mostly static, QC plots
    if (length(x)>10000) {
        ii <- sample(seq_len(length(x)),10000)
        x <- x[ii]
        y <- y[ii]
    }
    
    fit <- lowess(x,y)
    
    stMsg <- ""
    if (!is.null(status))
        stMsg <- status
    
    switch(jl,
        highcharts = {
            #if (is.null(altNames))
            #   point.format=paste("<strong>Gene ID: </strong>{point.name}<br>",
            #        "<strong>",xLab,": </strong>{point.x}<br>",
            #        "<strong>",yLab,": </strong>{point.y}",sep="")
            #else
            #    point.format=paste("<strong>Gene name: </strong>",
            #        "{point.alt_name}<br>",
            #        "<strong>Gene ID: </strong>{point.name}<br>",
            #        "<strong>",xLab,": </strong>{point.x}<br>",
            #        "<strong>",yLab,": </strong>{point.y}",sep="")
            point.format=paste("<strong>",xLab,": </strong>{point.x}<br>",
                "<strong>",yLab,": </strong>{point.y}",sep="")
                json <- list(
                    chart=list(
                    type="scatter",
                    zoomType="xy"
                ),
                title=list(
                    text=paste(type,"plot for",stMsg,"samples",samples[1],
                        "and",samples[2])
                ),
                xAxis=list(
                    title=list(
                        text=xLab,
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
                    min=round(min(x),3),
                    max=round(max(x),3)
                ),
                yAxis=list(
                    title=list(
                        useHTML=TRUE,
                        text=yLab,
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
                    min=round(min(y),3),
                    max=round(max(y),3)
                ),
                plotOptions=list(
                    scatter=list(
                        turboThreshold=25000
                    ),
                    line=list(
                        turboThreshold=25000
                    )
                ),
                series=list(
                    list(
                        type="scatter",
                        name="Genes/Transcripts",
                        #color="#00ABEE",
                        color="rgba(0,171,238,0.5)",
                        marker=list(
                            radius=0.5
                        ),
                        #data=makeHighchartsPoints(x,y,unname(altNames)),
                        data=makeHighchartsPoints(x,y,simple=TRUE),
                        tooltip=list(
                            followPointer=FALSE,
                            headerFormat=paste("<span style=",
                                "\"font-size:1.1em;color:{series.color};",
                                "font-weight:bold\">{series.name}<br>",
                                sep=""),
                            pointFormat=point.format
                        )
                    ),
                    list(
                        type="line",
                        name="Trend",
                        marker=list(
                            enabled=FALSE
                        ),
                        color="#E40000",
                        data=lapply(seq_len(length(x)),function(i,x,y) {
                            return(c(x[i],y[i])) 
                        },round(fit$x,3),round(fit$y,3))
                    )
                )
            )
        }
    )
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

rnacompToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- tolower(out[1])
    
    status <- obj$status
    dat <- obj$user$plotdata$data2plot
    samples <- obj$samples
    refColumn <- obj$user$plotdata$refColumn
    
    # Too many points for a lot of curves of interactive data
    if (nrow(dat)>1000) {
        ii <- sort(sample(seq_len(nrow(dat)),998))
        ii <- c(1,ii,nrow(dat))
        dat <- cbind(dat[ii,1],dat[ii,2:ncol(dat)])
    }

    if (is.list(samples))
        samplenames <- unlist(samples)
    
    samplenames <- setdiff(samplenames,refColumn)
    dat <- as.matrix(dat[,samplenames])
    datDens <- apply(dat,2,density,adjust=1.5)
    limY <- c(0,max(vapply(datDens,function (x) max(x$y,na.rm=TRUE),
        numeric(1))))  
    cols <- .getColorScheme(length(samplenames))
    
    # Construct series
    counter <- 0
    densSeries <- vector("list",length(samplenames))
    names(densSeries) <- samplenames
    for (n in names(densSeries)) {
        counter <- counter + 1
        densSeries[[n]] <- list()
        densSeries[[n]]$name=n
        densSeries[[n]]$color=cols$fill[counter]
        densSeries[[n]]$data <- 
            makeHighchartsPoints(datDens[[n]]$x,datDens[[n]]$y)
        densSeries[[n]]$tooltip=list(
            headerFormat=paste("<span style=",
                "\"font-size:1.1em;color:{series.color};",
                "font-weight:bold\">{series.name}<br>",
                sep=""),
            pointFormat=NULL
        )
    }
    counter <- 0
    abSeries <- vector("list",length(samplenames))
    names(abSeries) <- samplenames
    for (n in names(abSeries)) {
        counter <- counter + 1
        abSeries[[n]] <- list()
        abSeries[[n]]$name=paste(n,"median")
        abSeries[[n]]$color=cols$fill[counter]
        abSeries[[n]]$data <- 
            makeHighchartsPoints(rep(median(dat[,n],na.rm=TRUE),nrow(dat)),
                seq(limY[1],limY[2],length.out=nrow(dat)))
        abSeries[[n]]$enableMouseTracking=FALSE
        abSeries[[n]]$dashStyle="Dash"
    }
    
    switch(jl,
        highcharts = {
            json <- list(
                chart=list(
                    type="line",
                    zoomType="xy"
                ),
                title=list(
                    text=paste("RNA-composition for",status,"mapped reads")
                ),
                xAxis=list(
                    title=list(
                        text="log<sub>2</sub>(sample/reference sample)",
                        margin=20,
                        useHTML=TRUE,
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
                    gridLineWidth=1
                ),
                yAxis=list(
                    title=list(
                        text="Density",
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
                    min=limY[1],
                    max=limY[2]
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
                series=c(unname(densSeries),unname(abSeries))
            )
        }
    )
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

biodistToJSON <- function(obj,jl=c("highcharts"),by=c("chromosome","biotype"),
    out=c("json","list")) {
    jl <- tolower(jl[1])
    by <- by[1]
    out <- tolower(out[1])
    
    plotdataChromosome <- obj$user$plotdata$chromosome
    plotdataBiotype <- obj$user$plotdata$biotype
    
    cols <- .getColorScheme()
    
    if (by == "chromosome") {
        # Data series
        series <- vector("list",3)
        names(series) <- c("genome","chromosome_deg","deg_chromosome")
        series$genome <- list()
        series$genome$id <- "genome"
        series$genome$name <- "% chromosome in genome"
        series$genome$color <- cols$trans[1]
        series$genome$pointPlacement <- -0.2
        series$genome$data <- round(as.numeric(plotdataChromosome[1,]),3)
        series$chromosome_deg <- list()
        series$chromosome_deg$id <- "chromosome_deg"
        series$chromosome_deg$name <- "% chromosome in DEG"
        series$chromosome_deg$color <- cols$trans[2]
        series$chromosome_deg$pointPlacement <- 0
        series$chromosome_deg$data <- 
            round(as.numeric(plotdataChromosome[2,]),3)
        series$deg_chromosome <- list()
        series$deg_chromosome$id <- "deg_chromosome"
        series$deg_chromosome$name <- "% DEG in chromosome"
        series$deg_chromosome$color <- cols$trans[3]
        series$deg_chromosome$pointPlacement <- 0.2
        series$deg_chromosome$data <- 
            round(as.numeric(plotdataChromosome[3,]),3)
            
        json <- switch(jl,
            highcharts = {
                list(
                    chart=list(
                        type="column",
                        alignTicks=FALSE
                    ),
                    title=list(
                        text="DEG distribution across chromosomes"
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
                        categories=names(plotdataChromosome),
                        title=list(
                            text="Chromosome",
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
                            min=0,
                            max=ceiling(max(plotdataChromosome[1,])),
                            title=list(
                                text="% of features",
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
                    series=unname(series)
                )
            }
        )
    }
    if (by == "biotype") {
        # For biotype
        abu <- which(plotdataBiotype[1,]>7)
        nabu <- which(plotdataBiotype[1,]<=7)
        
        # Data series
        seriesAbu <- vector("list",3)
        names(seriesAbu) <- c("genome","biotype_deg","deg_biotype")
        seriesAbu$genome <- list()
        seriesAbu$genome$id <- "abu_genome"
        seriesAbu$genome$name <- "% biotype in genome"
        seriesAbu$genome$color <- cols$trans[1]
        seriesAbu$genome$pointPlacement <- -0.2
        seriesAbu$genome$data <- round(as.numeric(plotdataBiotype[1,abu]),3)
        seriesAbu$biotype_deg <- list()
        seriesAbu$biotype_deg$id <- "abu_biotype_deg"
        seriesAbu$biotype_deg$name <- "% biotype in DEG"
        seriesAbu$biotype_deg$color <- cols$trans[2]
        seriesAbu$biotype_deg$pointPlacement <- 0
        seriesAbu$biotype_deg$data <- 
            round(as.numeric(plotdataBiotype[2,abu]),3)
        seriesAbu$deg_biotype <- list()
        seriesAbu$deg_biotype$id <- "abu_deg_biotype"
        seriesAbu$deg_biotype$name <- "% DEG in biotype"
        seriesAbu$deg_biotype$color <- cols$trans[3]
        seriesAbu$deg_biotype$pointPlacement <- 0.2
        seriesAbu$deg_biotype$data <- 
            round(as.numeric(plotdataBiotype[3,abu]),3)
        seriesNabu <- vector("list",3)
        names(seriesNabu) <- c("genome","biotype_deg","deg_biotype")
        seriesNabu$genome <- list()
        seriesNabu$genome$name <- "% biotype in genome"
        seriesNabu$genome$yAxis <- 1
        seriesNabu$genome$pointStart <- length(abu)
        seriesNabu$genome$linkedTo <- "abu_genome"
        seriesNabu$genome$color <- cols$trans[1]
        seriesNabu$genome$pointPlacement <- -0.2
        seriesNabu$genome$data <- round(as.numeric(plotdataBiotype[1,nabu]),3)
        seriesNabu$biotype_deg <- list()
        seriesNabu$biotype_deg$name <- "% biotype in DEG"
        seriesNabu$biotype_deg$yAxis <- 1
        seriesNabu$biotype_deg$pointStart <- length(abu)
        seriesNabu$biotype_deg$linkedTo <- "abu_detected"
        seriesNabu$biotype_deg$color <- cols$trans[2]
        seriesNabu$biotype_deg$pointPlacement <- 0
        seriesNabu$biotype_deg$data <- 
            round(as.numeric(plotdataBiotype[2,nabu]),3)
        seriesNabu$deg_biotype <- list()
        seriesNabu$deg_biotype$name <- "% DEG in biotype"
        seriesNabu$deg_biotype$yAxis <- 1
        seriesNabu$deg_biotype$pointStart <- length(abu)
        seriesNabu$deg_biotype$linkedTo <- "abu_sample"
        seriesNabu$deg_biotype$color <- cols$trans[3]
        seriesNabu$deg_biotype$pointPlacement <- 0.2
        seriesNabu$deg_biotype$data <- 
            round(as.numeric(plotdataBiotype[3,nabu]),3)
        
        json <- switch(jl,
            highcharts = {
                list(
                    chart=list(
                        type="column",
                        alignTicks=FALSE
                    ),
                    title=list(
                        text="DEG distribution across biotypes"
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
                        categories=colnames(plotdataBiotype)[c(abu,nabu)],
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
                                to=ncol(plotdataBiotype)
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
                    series=c(unname(seriesAbu),unname(seriesNabu))
                )
            }
        )
    }
    
    if (out=="json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out=="list")
        return(json)
}

maStatToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- out[1]
    
    a <- obj$x
    f <- obj$y
    p <- obj$user$p
    xlim <- obj$xlim
    ylim <- obj$ylim
    pcut <- obj$pcut
    fcut <- obj$fcut
    altNames <- obj$altnames
    con <- obj$user$con
    #stat <- obj$user$stat
    
    conlab <- strsplit(con,"_vs_")
    conlab <- paste(conlab[[1]][1],"against",conlab[[1]][2])
    #conlab <- paste(conlab[[1]][2],"against",conlab[[1]][1])
    #statMap <- .getStatMap()
    #statlab <- statMap[[stat]]
    
    if (!is.null(p)) {
        upstat <- which(f>=fcut & p<pcut)
        downstat <- which(f<=-fcut & p<pcut)
        up <- which(f>=fcut & p>=pcut)
        down <- which(f<=-fcut & p>=pcut)
        poor <- which(p<pcut & abs(f)<fcut)
        neutral <- setdiff(seq_len(length(a)),
            Reduce("union",list(upstat,downstat,up,down,poor)))
    }
    else {
        upstat <- downstat <- poor <- integer(0)
        up <- which(f>=fcut)
        down <- which(f<=-fcut)
        neutral <- setdiff(seq_len(length(a)),union(up,down))
    }
    
    switch(jl,
        highcharts = {
            series <- list()
            counter <- 0
            if (length(upstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly up-regulated",
                    color="#EE0000",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=a[upstat],
                        y=f[upstat],
                        a=unname(altNames[upstat]),
                        p=if (!is.null(p)) -log10(p[upstat]) else NULL
                    )
                )
            }
            if (length(downstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly down-regulated",
                    color="#00EE00",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=a[downstat],
                        y=f[downstat],
                        a=unname(altNames[downstat]),
                        p=if (!is.null(p)) -log10(p[downstat]) else NULL
                    )
                )
            }
            if (length(up)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Up-regulated",
                    color="#B80000",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=a[up],
                        y=f[up],
                        a=unname(altNames[up]),
                        p=if (!is.null(p)) -log10(p[up]) else NULL
                    )
                )
            }
            if (length(down)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Down-regulated",
                    color="#008000",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=a[down],
                        y=f[down],
                        a=unname(altNames[down]),
                        p=if (!is.null(p)) -log10(p[down]) else NULL
                    )
                )
            }
            if (length(poor)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significant only",
                    color="#FFA600",
                    data=makeHighchartsPoints(
                        x=a[poor],
                        y=f[poor],
                        a=unname(altNames[poor]),
                        p=if (!is.null(p)) -log10(p[poor]) else NULL
                    )
                )
            }
            if (length(neutral)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Neutral",
                    color="#858585",
                    data=makeHighchartsPoints(
                        x=a[neutral],
                        y=f[neutral],
                        a=unname(altNames[neutral]),
                        p=if (!is.null(p)) -log10(p[neutral]) else NULL
                    )
                )
            }
        
            if (is.null(altNames))
                point.format=paste("<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Average expression: </strong>{point.x}<br>",
                    "<strong>Fold change: </strong>{point.y}<br>",
                    "<strong>Significance: </strong>{point.sig}",sep="")
            else
                point.format=paste("<strong>Gene name: </strong>",
                    "{point.alt_name}<br>",
                    "<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Average expression: </strong>{point.x}<br>",
                    "<strong>Fold change: </strong>{point.y}<br>",
                    "<strong>Significance: </strong>{point.sig}",sep="")
                json <- list(
                    chart=list(
                    type="scatter",
                    zoomType="xy"
                ),
                title=list(
                    text=paste("Mean-Difference (MA) plot for",conlab)
                ),
                xAxis=list(
                    title=list(
                        text="Average expression",
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
                        text="Fold change (log<sub>2</sub>)",
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
                        marker=list(
                            radius=2,
                            symbol="circle",
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
                        #events=list(
                        #   legendItemClick=paste("function() {",
                        #       "return false; }")
                        #),
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
                series=series
            )
        }
    )
    if (out == "json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out == "list")
        return(json)
}

dereguloToJSON <- function(obj,jl=c("highcharts"),out=c("json","list")) {
    jl <- tolower(jl[1])
    out <- out[1]
    
    xlim <- obj$xlim
    ylim <- obj$ylim
    pcut <- obj$pcut
    fcut <- obj$fcut
    altNames <- obj$altnames
    fmat <- obj$user$fmat
    pmat <- obj$user$pmat
    
    conlabX <- strsplit(colnames(fmat)[1],"_vs_")
    conlabX <- paste(conlabX[[1]][1],"against",conlabX[[1]][2])
    conlabY <- strsplit(colnames(fmat)[2],"_vs_")
    conlabY <- paste(conlabY[[1]][1],"against",conlabY[[1]][2])
    
    # red
    upupstat <- which(apply(fmat,1,function(x) all(x >= fcut)) &
        apply(pmat,1,function(x) all(x < pcut)))
    # green
    downdownstat <- which(apply(fmat,1,function(x) all(x <= -fcut)) &
        apply(pmat,1,function(x) all(x < pcut)))
    # red3
    upup <- which(apply(fmat,1,function(x) all(x >= fcut)) &
        apply(pmat,1,function(x) any(x >= pcut)))
    # green3
    downdown <- which(apply(fmat,1,function(x) all(x <= -fcut)) &
        apply(pmat,1,function(x) any(x >= pcut)))
    # orange
    updownstat <- which(apply(fmat,1,
        function(x) x[1] >= fcut & x[2] <= -fcut) &
            apply(pmat,1,function(x) all(x < pcut)))
    # blue
    downupstat <- which(apply(fmat,1,
        function(x) x[1] <= -fcut & x[2] >= fcut) &
            apply(pmat,1,function(x) all(x < pcut)))
    # orange3
    updown <- which(apply(fmat,1,function(x) x[1] >= fcut & x[2] <= -fcut) &
        apply(pmat,1,function(x) any(x >= pcut)))
    # blue3
    downup <- which(apply(fmat,1,function(x) x[1] <= -fcut & x[2] >= fcut) &
        apply(pmat,1,function(x) any(x >= pcut)))
    # black
    poor <- which(apply(pmat,1,function(x) all(x < pcut)) &
        apply(fmat,1,function(x) any(abs(x) < fcut)))
    # gray70
    nones <- which(apply(pmat,1,function(x) any(x >= pcut)) &
        apply(fmat,1,function(x) all(abs(x) < fcut)))
    # gray40
    neutral <- setdiff(seq_len(nrow(fmat)),
        Reduce("union",list(upupstat,downdownstat,upup,downdown,updownstat,
            downupstat,updown,downup,poor,nones)))
    
    switch(jl,
        highcharts = {
            series <- list()
            counter <- 0
            if (length(upupstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly both up-regulated",
                    color="#EE0000",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[upupstat,1],
                        y=fmat[upupstat,2],
                        a=unname(altNames[upupstat])
                    )
                )
            }
            if (length(downdownstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly both down-regulated",
                    color="#00EE00",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[downdownstat,1],
                        y=fmat[downdownstat,2],
                        a=unname(altNames[downdownstat])
                    )
                )
            }
            if (length(upup)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Both up-regulated",
                    color="#B80000",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[upup,1],
                        y=fmat[upup,2],
                        a=unname(altNames[upup])
                    )
                )
            }
            if (length(downdown)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Both down-regulated",
                    color="#008000",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[downdown,1],
                        y=fmat[downdown,2],
                        a=unname(altNames[downdown])
                    )
                )
            }
            if (length(updownstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly up-down-regulated",
                    color="#FFB700",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[updownstat,1],
                        y=fmat[updownstat,2],
                        a=unname(altNames[updownstat])
                    )
                )
            }
            if (length(downupstat)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significantly down-up-regulated",
                    color="#0000FF",
                    marker=list(
                        radius=3
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[downupstat,1],
                        y=fmat[downupstat,2],
                        a=unname(altNames[downupstat])
                    )
                )
            }
            if (length(updown)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Up-down-regulated",
                    color="#C68100",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[updown,1],
                        y=fmat[updown,2],
                        a=unname(altNames[updown])
                    )
                )
            }
            if (length(downup)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Down-up-regulated",
                    color="#0000AB",
                    marker=list(
                        radius=2
                    ),
                    data=makeHighchartsPoints(
                        x=fmat[downup,1],
                        y=fmat[downup,2],
                        a=unname(altNames[downup])
                    )
                )
            }
            if (length(poor)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Significant only",
                    color="#1F1F1F",
                    data=makeHighchartsPoints(
                        x=fmat[poor,1],
                        y=fmat[poor,2],
                        a=unname(altNames[poor])
                    )
                )
            }
            if (length(neutral)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="Neutral",
                    color="#00DADE",
                    data=makeHighchartsPoints(
                        x=fmat[neutral,1],
                        y=fmat[neutral,2],
                        a=unname(altNames[neutral])
                    )
                )
            }
            if (length(nones)>0) {
                counter <- counter + 1
                series[[counter]] <- list(
                    name="No regulation",
                    color="#CCCCCC",
                    data=makeHighchartsPoints(
                        x=fmat[nones,1],
                        y=fmat[nones,2],
                        a=unname(altNames[nones])
                    )
                )
            }
            
            # Add fold lines
            counter <- counter + 1
            series[[counter]] <- list(
                name="X up fold threshold",
                color="#000000",
                type="line",
                dashStyle="Dash",
                lineWidth=1,
                marker=list(
                    enabled=FALSE
                ),
                tooltip=list(
                    headerFormat=paste('<strong>{series.name}',
                        '</strong><br/>',sep=""),
                    pointFormat=paste('<strong>Threshold: ',
                        '</strong>{point.x}<br/>',sep="")
                ),
                data=list(round(c(fcut,floor(ylim[1])),3),
                    round(c(fcut,ceiling(ylim[2])),3))
            )
            counter <- counter + 1
            series[[counter]] <- list(
                name="X down fold threshold",
                color="#000000",
                type="line",
                dashStyle="Dash",
                lineWidth=1,
                marker=list(
                    enabled=FALSE
                ),
                tooltip=list(
                    headerFormat=paste('<strong>{series.name}',
                        '</strong><br/>',sep=""),
                    pointFormat=paste('<strong>Threshold: ',
                        '</strong>{point.x}<br/>',sep="")
                ),
                data=list(round(c(-fcut,floor(ylim[1])),3),
                    round(c(-fcut,ceiling(ylim[2])),3))
            )
            counter <- counter + 1
            series[[counter]] <- list(
                name="Y up fold threshold",
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
                        '</strong>{point.y}<br/>',sep="")
                ),
                data=list(round(c(floor(xlim[1]),fcut),3),
                    round(c(ceiling(xlim[2]),fcut),3))
            )
            counter <- counter + 1
            series[[counter]] <- list(
                name="Y down fold threshold",
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
                        '</strong>{point.y}<br/>',sep="")
                ),
                data=list(round(c(floor(xlim[1]),-fcut),3),
                    round(c(ceiling(xlim[2]),-fcut),3))
            )
        
            if (is.null(altNames))
                point.format=paste("<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change X: </strong>{point.x}<br>",
                    "<strong>Fold change Y: </strong>{point.y}<br>",sep="")
            else
                point.format=paste("<strong>Gene name: </strong>",
                    "{point.alt_name}<br>",
                    "<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change X: </strong>{point.x}<br>",
                    "<strong>Fold change Y: </strong>{point.y}<br>",sep="")
                json <- list(
                    chart=list(
                    type="scatter",
                    zoomType="xy"
                ),
                title=list(
                    text=paste("Deregulogram")
                ),
                xAxis=list(
                    title=list(
                        text=paste0("Fold change (log<sub>2</sub>)",conlabX),
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
                        text=paste0("Fold change (log<sub>2</sub>)",conlabY),
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
                        marker=list(
                            radius=2,
                            symbol="circle",
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
                series=series
            )
        }
    )
    if (out == "json")
        return(toJSON(json,auto_unbox=TRUE,null="null"))
    else if (out == "list")
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
            for (i in seq_len(length(starts)))
                substr(js,starts[i]-1,starts[i]-1) <- " "
            ends <- as.numeric(cl[[1]])
            for (i in seq_len(length(starts)))
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
        return(x[seq_len(length(levels(classes)))][design])
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

.getStatMap <- function() {
    return(list(
        deseq="DEseq",
        deseq2="DEseq2",
        edger="edgeR",
        limma="voom",
        nbpseq="NBPSeq",
        noiseq="NOISeq",
        bayseq="baySeq",
        absseq="ABSSeq",
        dss="DSS"
    ))
}
