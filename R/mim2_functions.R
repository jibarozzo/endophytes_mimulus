# Functions for the analysis of the data from the MIM2 project
# Bolívar Aponte Rolón


# Formatting formula for p value
format.p <- function(p, precision = .001) {
    digits <- -log(precision, base = 10)
    p <- formatC(p, format = 'f', digits = digits)
    if (p < .001) {
        p = paste0('< ', precision)
    }
    if (p >= .001) {
        p = paste0('= ', p)
    }
    sub("", "", p)
}


# Modified functions from mirlyn package #
# The objective of this script expands the the alpha diversity metrics 
# to include shannon, simpson and invsimpson indices. Potentially observed richness.

#Renaming the function to avoid conflicts with the original function
rep_otu_df <- function(x) {
    
    newdat <- lapply(x, otu_table)
    
    finaldat <- matrix(ncol = ncol(newdat[[1]]) * length(newdat), nrow = nrow(newdat[[1]]),
                       dimnames = list(rownames(newdat[[1]]),
                                       paste0(rep(seq_len(length(newdat)), each = ncol(newdat[[1]])), "-",
                                              rep(colnames(newdat[[1]]), times = length(newdat)))))
    
    for (i in seq_along(newdat)) {
        finaldat[, (1 + (i - 1) * ncol(newdat[[1]])):(i * ncol(newdat[[1]]))] <- newdat[[i]]
    }
    
    finaldat
    
}



alpha_rfy_DF <- function(x, diversity = c("shannon", "simpson", "invsimpson")) {
    md <- sample_data(x[[1]])
    
    # Calculate observed richness
    t_otu_table <- t(rep_otu_df(x))
    observed <- specnumber(t_otu_table)
    observed_richness_df <- rownames_to_column(as.data.frame(observed), var = "Unique_ID")
    
    div_df_list <- lapply(diversity, function(index) {
        div_values <- vegan::diversity(t(rep_otu_df(x)), index = index)
        if (!is.matrix(div_values)) {
            div_values <- as.matrix(div_values, ncol = 1, dimnames = list(rownames(div_values), NULL))
        }
        colnames(div_values) <- paste0(index, "_", colnames(div_values))
        return(data.frame(Unique_ID = rownames(div_values), div_values, stringsAsFactors = FALSE))
    })
    
    final <- reduce(div_df_list, function(df1, df2) inner_join(df1, df2, by = "Unique_ID"))
    
    # Merge observed richness with the final dataframe
    final <- merge(final, observed_richness_df, by = "Unique_ID", all.x = TRUE)
    
    final <- cbind(md, final)
    
    return(final)
}



# Example of use
#alpha_rfy_DF(pseq_rrf5, diversity = c("shannon", "simpson", "invsimpson")) # Alpha diversity indices


# Modified version of function tag_facet from "egg" package.
# This function allows me to add a tag fro p-values of anything without losing the facet strip labels.


tag_fat <-function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                   hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}


# Modified function from the "ggpubr" package
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
