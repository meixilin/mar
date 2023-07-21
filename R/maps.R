
ggraster_table <- function(ras) {
    # https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster
    # baseraste<-makebaseraster()
    # Euroclim<-make_Euroclim()
    example <- ras
    test_spdf <- as(example, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    colnames(test_df) <- c("value", "x", "y")
    return(test_df)
}
