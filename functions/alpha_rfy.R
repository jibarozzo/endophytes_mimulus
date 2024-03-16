# Modified functions from mirlyn package #
# The objective of this script expands the the alphha diversity metrics to include shannon, simpson and invsimpson indices. Poteially observed richness.

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