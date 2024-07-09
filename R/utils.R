# base::match but with some checks
# TODO: orphaned function not used.
.match <- function(x, xdb) {
    stopifnot(all(x %in% xdb))
    # sort the ids if provided in the wrong orders
    ids = sort(base::match(x, xdb))
    return(ids)
}



################################################################################
#### Function utilities
################################################################################
fc <- function(x) as.character(as.matrix(x))
fn <- function(x) unlist(as.numeric(as.matrix(x)))
transparent <- function(col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) {
        rgb(
            c[1] / 255, c[2] / 255,
            c[3] / 255, alpha
        )
    })
    return(res)
}
normalize <- function(x) {
    if (class(x) != "numeric") {
        message("non-numeric data provided")
        message("attempting to force numeric")
        x <- fn(x)
    }
    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}

# create a table that ggplot can plot using tiles
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
ggraster_plot <- function(ras, mycolor = brewer.pal(8, "Greys")[-1]) {
    # first the raster object in stransformed into a long data table
    rtile <- ggraster_table(ras)
    # then plot using tiles
    p <- ggplot(rtile) +
        geom_tile(aes(x = x, y = y, fill = value)) +
        scale_fill_gradientn(colours = mycolor) +
        # coord_equal()+
        labs(y = "Latitude", x = "Longitude")
    return(p)
}

substrRight <- function(x, lastpos, giveright = T) {
    if (giveright == T) {
        substr(x, nchar(x) - lastpos + 1, nchar(x))
    } else {
        substr(x, 1, nchar(x) - lastpos + 1)
    }
}
