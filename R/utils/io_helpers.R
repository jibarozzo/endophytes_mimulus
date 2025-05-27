#' Initialize project paths for the Mimulus endophytes project
#'
#' Creates a list of standardized file paths for the project structure.
#' This function does not create directories, it only defines the paths.
#'
#' @param base_path Optional. Base directory path. If NULL, current working directory is used.
#' @return A list containing paths to key project directories
#' @export
initialize_project_paths <- function(base_path = NULL) {
    # If base_path is not provided, use current working directory
    if (is.null(base_path)) {
        base_path <- getwd()
    }
    
    # Create a list of standard paths used throughout the project
    list(
        # Main data directories
        data = file.path(base_path, "data"),
        clean_data = file.path(base_path, "data", "clean_data"),
        field_data = file.path(base_path, "data", "field_data"),
        filtered = file.path(base_path, "data", "filtered"),
        preprocess = file.path(base_path, "data", "preprocess"),
        
        # Sub-directories for clean data
        asv_tables = file.path(base_path, "data", "clean_data", "ASV_tables"),
        statistics = file.path(base_path, "data", "clean_data", "statistics"),
        taxonomy = file.path(base_path, "data", "clean_data", "taxonomy"),
        
        # Filtered data sub-directories
        filt_8450 = file.path(base_path, "data", "filtered", "filt_8450"),
        filt_8756 = file.path(base_path, "data", "filtered", "filt_8756"),
        
        # Preprocess sub-directories
        pre_8450 = file.path(base_path, "data", "preprocess", "pre_8450"),
        pre_8756 = file.path(base_path, "data", "preprocess", "pre_8756"),
        
        # Output directories
        figures = file.path(base_path, "figures"),
        tables = file.path(base_path, "tables")
    )
}
