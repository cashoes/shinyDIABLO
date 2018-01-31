
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinythemes)
library(plotly)

library(network)
library(sna)
library(igraph)
library(intergraph)
library(visNetwork)

library(mixOmics)
library(ggmixOmics)

library(tidyverse)
library(cowplot)
library(GGally)
library(ggnetwork)

source('helpers.R')

server <- function(input,output, session) {

  # Compplot ----
  output$compplot1 <- renderPlotly({
    p <- get(input$selectDataComp, ggmixOmics::ggcompplot(M, comps = c(as.numeric(input$compXComp),
                                                                       as.numeric(input$compYComp))))

    ggplotly(p) %>%
      layout(dragmode = "lasso", height = 800, width = 800)
  }
  #, height = function() {
  #  session$clientData$output_compplot1_width
  #}
  )

  output$hoverComp <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  })

  output$clickComp <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click events appear here (double-click to clear)" else d
  })

  output$brushComp <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })


  # Varplot ----
  plotVar <- mixOmics::plotVar(model1)

  # Get rownames for each Block
  # legacy
  # flowCytometry <- row.names(subset(plotVar, Block == "Flow Cytometry"))
  # luminexCytokine <- row.names(subset(plotVar, Block == "Luminex Cytokine"))
  # metabolomics <- row.names(subset(plotVar, Block == "Metabolomics"))
  # proteomics <- row.names(subset(plotVar, Block == "Proteomics"))
  # transcriptomics <- row.names(subset(plotVar, Block == "Transcriptomics"))

  output$var1 <- renderPlotly({
    p <- get(input$selectDataVar, ggmixOmics::ggvarplot(M))

    ggplotly(p) %>%
      layout(dragmode = "lasso")
  })

  output$varTable <- DT::renderDataTable({
    compN <- paste(c("comp ", input$selectComp), sep="", collapse="")
    table <- as.matrix(model1$loadings[[1]][,compN])
    for(i in 2:nEntries){
      table <- rbind(table, as.matrix(model1$loadings[[i]][,compN]))
    }
    if(input$compare == TRUE){
      tempTable <- as.matrix(model2$loadings$`Flow cytometry`[,compN])
      tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Luminex cytokine'[,compN]))
      tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Metabolomics'[,compN]))
      tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Proteomics'[,compN]))
      tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Transcriptomics'[,compN]))
      table <- cbind(as.matrix(table), as.matrix(tempTable))
    }
    table

    if(input$compare == FALSE)
      colnames(table) <- c('Eigenvector')
    if(input$compare == TRUE)
      colnames(table) <- c('Eigenvector 1', 'Eigenvector 2')

    # table <- as.matrix(table[apply(table[,-1], 1, function(x) !all(x==0)),])
    DT::datatable(data = table,
                  options = list(scrollX = TRUE, scrollY = "275px", autoWidth = FALSE)
    )
  })

  # Biplot ----
  # Make names (TODO)
  samples <- rownames(M$X[[1]])
  output$biplot1 <- renderPlotly({
    p <- get(input$selectDataBi, ggmixOmics::ggbiplot(M, comps = c(as.numeric(input$compXBi),
                                                                   as.numeric(input$compYBi))))
    p$data <- cbind(p$data, samples)

    ggplotly(p) %>%
      layout(dragmode = "lasso")
  })

  # Loadings ----
  output$loadings1 <- renderPlot({mixOmics::plotLoadings(model1, title = "Plot of Loading vectors 1")$graph},
                                 bg = "transparent")
  output$loadings2 <- renderPlot({mixOmics::plotLoadings(model2, title = "Plot of Loading vectors 2")$graph},
                                 bg = "transparent")

  # Circos ----
  output$circos1 <- renderPlot({mixOmics::circosPlot(model1, cutoff = input$cutoff)}, bg = "transparent")
  output$circos2 <- renderPlot({mixOmics::circosPlot(model2, cutoff = input$cutoff)}, bg = "transparent")


  # Heatmap ----
  ncomp <- 2
  corr <- getCorMat(M, method = 'pearson')
  d <- hclust(dist(corr))
  m <- corr[d$order, d$order]
  g1 <- GGally::ggcorr(data = NULL, cor_matrix = m, hjust = 0, size = 0.1, colour = 'white', layout.exp = 0)

  output$heatmap1 <- renderPlot({plot(g1)})

  # Diablo ----
  output$diablo1 <- renderPlot({mixOmics::plotDiablo(model1)$graph}, bg = "transparent")
  output$diablo2 <- renderPlot({mixOmics::plotDiablo(model2)$graph}, bg = "transparent")

  # Network ----

  corThreshold <- reactive(input$threshold) %>% debounce(1000)
  output$network <- renderPlotly({

    # corThreshold <- input$threshold
    # corMat[corMat < corThreshold()] <- 0

    # rownames(corMat) <- make.names(rownames(corMat), unique=TRUE)
    # colnames(corMat) <- make.names(colnames(corMat), unique=TRUE)

    # corMat <- abs(getCorMat(M))
    graph <- igraph::graph_from_adjacency_matrix(corMat, mode = 'lower', weighted = 'weight', diag = F)
    graph <- igraph::delete_edges(graph, E(graph)[E(graph)$weight < corThreshold()])

    # cluster <- igraph::cluster_fast_greedy(graph)

    # # figure out optimal cor threhsold by modularity
    # range <- seq(0.25, 0.95, by = 0.05)
    # purrr::walk(range, g = graph, ~ {
    #   g <- list(...)$g
    #   g <- igraph::delete_edges(g, E(g)[E(g)$weight < .])
    #
    #   s <- sprintf('Cor: %.2f | Num: %03d | Mod: %.2f', ., igraph::clique_num(g), igraph::modularity(igraph::cluster_fast_greedy(g)))
    #   print(s)
    #
    #   # k <- igraph::cluster_fast_greedy(g)
    #   # m <- igraph::modularity(k)
    #   # l <- igraph::groups(k) %>% length()
    #   # list(g = g, k = k, m = m, l = l, thresh = .)
    # })
    #
    # mods <- map(clust, 'm') %>% unlist()
    # nums <- map(clust, 'l') %>% unlist()
    #
    # plot(nums, mods)
    #
    # clust <- clust[[which.max(sign(diff(diff(mods[order(nums)])))[-1])]]
    # graph <- clust$g


    # graph <- graph.adjacency(abs(corMat), weighted = TRUE, mode = "lower")
    # graph <- simplify(graph)

    # Assume ordered
    keeps <- 1:unique(M$ncomp) %>%
      purrr::map(~ mixOmics::selectVar(M, comp = .)) %>%
      purrr::at_depth(2, ~ .x[[1]]) %>%
      purrr::transpose() %>%
      purrr::map(purrr::reduce, union) %>%
      purrr::map(length) %>%
      head(-1)

    V(graph)$group <- rep(names(keeps), keeps)

    # graph information
    output$density <- renderValueBox({
      valueBox(
        value = round(edge_density(graph, loops = FALSE), digits = 3),
        subtitle = "Edge Density",
        icon = icon("anchor")
      )
    })

    output$transitivity <- renderValueBox({
      valueBox(
        value = round(transitivity(graph, type = "global", vids = NULL, weights = NULL, isolates = c("NaN", "zero")), digits = 3),
        subtitle = "Transitivity",
        icon = icon("wifi")
      )
    })

    output$modularity <- renderValueBox({
      valueBox(
        # value = round(modularity(graph, membership(cluster_edge_betweenness(graph))), digits = 3),
        value = round(modularity(graph, membership(cluster_fast_greedy(graph))), digits = 3),
        subtitle = "Modularity",
        icon = icon("gavel")
      )
    })

    # convert plot to ggnetwork
    nodesNedges <- ggnetwork(graph, layout = 'kamadakawai')
    nodesNedges$xend <- as.numeric(nodesNedges$xend)
    nodesNedges$yend <- as.numeric(nodesNedges$yend)
    nodesNedges$y <- as.numeric(nodesNedges$y)
    nodesNedges$x <- as.numeric(nodesNedges$x)

    nodes <- nodesNedges[is.na(nodesNedges$weight), c("x", "y", "group" ,"vertex.names")]
    edges <- nodesNedges[!is.na(nodesNedges$weight), c("x", "y", "vertex.names", "xend", "yend")]

    linkNode <- merge(edges, nodes, by.x = c("xend", "yend"), by.y = c("x","y"))

    nodesNedges <- merge(nodesNedges, linkNode, all.x = TRUE)
    nodesNedges <- nodesNedges[order(nodesNedges[,"na.y"], nodesNedges[, "x"], na.last = FALSE),]

    # reorder
    nodesNedges <- nodesNedges[,c("x", "y", "group", "na.x", "vertex.names", "xend", "yend", "na.y", "weight", "vertex.names.x", "vertex.names.y")]

    # PPI Integration ----
    if(PPIIntegration == TRUE){
      # import PPI data
      data <- PPIList

      # match edges and PPI data
      matches <- matchPPI(nodesNedges[, c("vertex.names.x", "vertex.names.y")], data)
      nodesNedges <- cbind(nodesNedges, matches)
    }


    # plotly interactivity
    toPrint <- NULL
    for(i in 1:nEntries){
      toPrint[[length(toPrint)+1]] <- nodesNedges[nodesNedges$group == dataNames[i],]
    }

    # output$hoverNetM <- renderValueBox({
    #   n <- event_data("plotly_hover")
    #   if (is.null(n)){
    #     value <- "N/A"
    #   }
    #   else
    #     value <- as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
    #   valueBox(
    #     value = value,
    #     subtitle = "Hovered Node"
    #   )
    # })

    # output$clickNetM <- renderValueBox({
    #   n <- event_data("plotly_click")
    #   if (is.null(n)){
    #     value <- "N/A"
    #   }
    #   else
    #     value <- as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
    #   valueBox(
    #     value = value,
    #     subtitle = "Clicked Node"
    #   )
    # })

    output$hoverNet <- renderPrint({
      n <- event_data("plotly_hover")
      if (is.null(n)) "Hover events appear here (unhover to clear)"
      else
        as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
    })

    output$clickNet <- renderPrint({
      n <- event_data("plotly_click")
      if (is.null(n)) "Click events appear here (double-click to clear)"
      else
        as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
    })

    output$brushNet <- renderDataTable({

      n <- event_data("plotly_selected")

      if (is.null(n)) return(NULL)

      n <- n[c("curveNumber", "pointNumber")]
      n[,2] <- n[,2] + 1

      library(sear)
      sear::collections

      print <- nodes %>%
        dplyr::tbl_df() %>%
        dplyr::mutate(group = factor(group, levels = dataNames), curveNumber = as.numeric(group)) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(pointNumber = 1:n()) %>%
        dplyr::inner_join(n, by = c('curveNumber', 'pointNumber')) %>%
        dplyr::pull('vertex.names') %>%
        as.character()

      table <- sear::sear(print) %>%
        dplyr::select(collection, subcollection, geneset, fdr) %>%
        dplyr::arrange(fdr) %>%
        dplyr::slice(1:100) %>%
        as.data.frame()

      # Geneset Enrichment ----
      #   print('Submitting to sear...')
      #   p <- genesetEnrichment(as.vector(print$vertex.names))
      # }
      # if(T){
      # else{
      #   p <- as.vector(print$vertex.names)
      # }

      # DT::datatable(
      #   data = genesetEnrichment(as.vector(print$vertex.names)),
      #   class = 'compact stripe',
      #   # colnames = c('Rank' = 1, 'Collection' = 2, 'Geneset' = 3, 'FDR' = 4),
      #   colnames = c('Collection' = 1, 'Geneset' = 2, 'FDR' = 3),
      #   caption = htmltools::tags$caption(
      #     style = 'caption-side: top; text-align: center;',
      #     'Table 1: ', htmltools::em('Geneset Enrichment Results')
      #   ), options = list(searching = FALSE, autoWidth = TRUE, scrollY = 400, paging = FALSE)
      # )
    }, options = list('lengthMenu' = c(10, 25, 50), 'pageLength' = 10))

    # output$brushNet <- DT::renderDataTable({
    #
    #   n <- event_data("plotly_selected")
    #
    #   print(n)
    #
    #   if (is.null(n)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)"
    #   else {
    #     n <- n[c("curveNumber", "pointNumber")]
    #     n[,2] <- n[,2] + 1
    #
    #     print <- nodes %>%
    #       dplyr::tbl_df() %>%
    #       dplyr::mutate(group = factor(group, levels = dataNames),
    #                     curveNumber = as.numeric(group)) %>%
    #       dplyr::group_by(group) %>%
    #       dplyr::mutate(pointNumber = 1:n()) %>%
    #       as.data.frame() %>%
    #       dplyr::inner_join(., n)
    #
    #     # Geneset Enrichment ----
    #     if(T){
    #       print('Submitting to sear...')
    #       p <- genesetEnrichment(as.vector(print$vertex.names))
    #     }
    #     else{
    #       p <- as.vector(print$vertex.names)
    #     }
    #
    #     p <- DT::datatable(
    #       data = p,
    #       class = 'compact stripe',
    #       # colnames = c('Rank' = 1, 'Collection' = 2, 'Geneset' = 3, 'FDR' = 4),
    #       colnames = c('Collection' = 2, 'Geneset' = 3, 'FDR' = 4),
    #       caption = htmltools::tags$caption(
    #         style = 'caption-side: top; text-align: center;',
    #         'Table 1: ', htmltools::em('Geneset Enrichment Results')
    #       ),
    #       options = list(searching = FALSE, autoWidth = TRUE, scrollY = 400, paging = FALSE)
    #     )
    #     print(str(p))
    #     p
    #   }
    # })

    output$brushNodes <- renderPrint({
      n <- event_data("plotly_selected")

      if (is.null(n)) return(NULL)

      n <- n[c("curveNumber", "pointNumber")]
      n[,2] <- n[,2] + 1

      print <- nodes %>%
        dplyr::tbl_df() %>%
        dplyr::mutate(group = factor(group, levels = dataNames),
                      curveNumber = as.numeric(group)) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(pointNumber = 1:n()) %>%
        as.data.frame() %>%
        dplyr::inner_join(n, c('curveNumber', 'pointNumber'))

      # Geneset Enrichment ----
      p <- as.vector(print$vertex.names)

      p
    })

    # plot graph ----
    # library dependency bug
    library(ggnetwork)

    # for (i in 1:nrow(nodesNedges)){
    #   if(!is.na(nodesNedges[i, "weight"])){
    #     nodesNedges[i, "vertex.names"] <- NA
    #   }
    # }

    plot <- ggplot2::ggplot(nodesNedges, ggplot2::aes(x = x, y= y, xend = xend, yend = yend, text = vertex.names)) +
      # ggnetwork::geom_edges(size = 0.1, color = "grey50") +
      ggnetwork::geom_nodes(ggplot2::aes(fill = group), size = 6, shape = 21, color = 'white') +
      ggthemes::scale_fill_few() +
      # viridis::scale_fill_viridis('', discrete = TRUE) +
      ggplot2::theme_void()

    ggplot <- ggplotly(plot, tooltip = "text") %>%
      layout(dragmode = "lasso")

    ggplot$x$data[[1]]$hoverinfo <- "none"

    ggplot
  })
}

# shinyServer(server)
