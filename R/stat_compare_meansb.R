stat_compare_meansb <- function(mapping = NULL, data = NULL,
         method = NULL, paired = FALSE, method.args = list(), ref.group = NULL,
         comparisons = NULL, hide.ns = FALSE, label.sep = ", ",
         label = NULL, label.x.npc = "left", label.y.npc = "top",
         label.x = NULL, label.y = NULL, vjust = 0, tip.length = 0.03,
         bracket.size = 0.3, step.increase = 0,
         symnum.args = list(),
         geom = "text", position = "identity",  na.rm = FALSE, show.legend = NA,
         inherit.aes = TRUE, ...) {
    
    if(!is.null(comparisons)){
        
        method.info <- .method_info(method)
        method <- method.info$method
        
        method.args <- .add_item(method.args, paired = paired)
        
        pms <- list(...)
        size <- ifelse(is.null(pms$size), 3.88, pms$size)
        color <- ifelse(is.null(pms$color), "black", pms$color)
        
        map_signif_level <- FALSE
        if(is.null(label)) label <- "p.format"
        
        if(.is_p.signif_in_mapping(mapping) | (label %in% "p.signif"))
        {
            map_signif_level <- c("****"=0.0001, "***"=0.001, "**"=0.01,  "*"=0.05, "ns"=Inf)
            if(hide.ns) map_signif_level <- map_signif_level[!names(map_signif_level) %in% "ns"]
            
        }
        
        if(!.is_empty(symnum.args)){
            
            symnum.args.isok <- length(symnum.args$cutpoints == length(symnum.args$symbols))
            if(!symnum.args.isok)
                stop("Incorrect format detected in symnum.args. ",
                     "Check the documentation.")
            map_signif_level <- symnum.args$cutpoints[-1] # the first element is 0 (the minimum p-value)
            names(map_signif_level) <- symnum.args$symbols
            if(hide.ns) map_signif_level <- map_signif_level[!names(map_signif_level) %in% "ns"]
            
        }
        
        if(missing(step.increase)){
            step.increase <- ifelse(is.null(label.y), 0.12, 0)
        }
        ggsignif::geom_signif(comparisons = comparisons, y_position = label.y,
                              test = method, test.args = method.args,
                              step_increase = step.increase, size = bracket.size, textsize = size, color = color,
                              map_signif_level = map_signif_level, tip_length = tip.length,
                              data = data, vjust = vjust)
    }
    
    else{
        mapping <- .update_mapping(mapping, label)
        layer(
            stat = StatCompareMeans, data = data, mapping = mapping, geom = geom,
            position = position, show.legend = show.legend, inherit.aes = inherit.aes,
            params = list(label.x.npc  = label.x.npc , label.y.npc  = label.y.npc,
                          label.x = label.x, label.y = label.y, label.sep = label.sep,
                          method = method, method.args = method.args,
                          paired = paired, ref.group = ref.group,
                          symnum.args = symnum.args,
                          hide.ns = hide.ns, na.rm = na.rm, vjust = vjust,...)
        )
        
    }
    
}
#bytecode: 0x5b05524243f0
#environment: namespace:ggpubr
