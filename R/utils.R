#' Interpolation of Missing Timesteps in a Raster Stack
#'
#' Interpolates missing time points in a raster stack with layers for different time points. Intended
#' for `terra::SpatRaster` use only.
#'
#' @import terra
#' @import purrr
#' @param raster_stack A `terra::SpatRaster` with at least two layers.
#' @param source_time A numeric vector indicating the timesteps of the `raster_stack`, whatever those
#' may be, e.g., `c(1940, 1945)`.
#' @param target_time A numeric vector indicating the series of timesteps desired in the output, e.g.,
#' `c(1940, 1941, 1942. 1943, 1944, 1945)`.
#' @param time_label `terra::SpatRaster`s do not allow raw numbers as names of raster layers. Therefore
#' a string is needed to place before the timestep number, e.g., "BP", "BCE".
#' @param ... Does nothing. A placeholder for future code improvements.
#' @param method Interpolation method passed on to `terra::approximate`. Default is "linear";
#' alternative is "constant", which implements a step function.
#' @return A `terra::SpatRaster` with as many layers as the length of `target_time`.
#' @export
interpolate_raster <- function(raster_stack, source_time, target_time, time_label, ...,
                               method = "linear") {
  if (!inherits(raster_stack, "SpatRaster")) {
    stop("raster_stack must be a terra::SpatRaster")
  }
  if (!all(is.numeric(source_time), is.numeric(target_time))) {
    stop("source time and target time must be numeric vectors")
  }
  if (length(target_time) < length(source_time)) {
    stop("Target time vector must be longer than source time vector
         \n (otherwise there's no interpolation)")
  }
  if (!is.character(time_label)) {
    stop("time_label must be a string")
  }
  template <- raster_stack[[1]]
  # Create output raster stack
  outputStack <- rast(nlyrs = length(target_time), nrows = nrow(template),
                      ncols = ncol(template),
                      xmin = xmin(template), xmax = xmax(template), ymin = ymin(template),
                      ymax = ymax(template), crs = crs(template), resolution = res(template),
                      vals = NA, names = target_time %>% map_chr(~paste0(time_label, .)))
  for (i in seq_along(source_time)) {
    j <- which(target_time == source_time[i])
    time <- target_time[j]
    inSeq <- which(source_time == time)
    if (length(inSeq) != 0) {
      r <- raster_stack[[inSeq]]
      values(outputStack[[j]]) <- r[]
    }
  }
  # urbanStack <- setZ(urbanStack, z = ts, name = "Years BP")
  interpolated_urban <- approximate(outputStack, method = method)
  return(interpolated_urban)
}
