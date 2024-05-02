# Modified version of function tag_facet from "egg" package.
# This function allows me to add a tag fro p-values of anything without losing the facet strip labels.
# Bolívar Aponte Rolón 2024/05/02

tag_fat <-function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                    hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}