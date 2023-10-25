#' R6 class representing a disease model of the *Mycoplasma gallisepticum*
#' outbreak in *Haemorhous mexicanus*
#'
#' @description \code{\link[R6:R6Class]{R6}} class representing fixed settings
#' for a spatially-explicit demographic-based SIRI model of disease dynamics. It
#' extends the \code{\link{SimulationModel}} class with parameters for the
#' \code{\link{disease_simulator}} function. It inherits functionality for
#' creating a nested model, whereby a nested template model with fixed
#' parameters is maintained when a model is cloned for various sampled
#' parameters. Also provided are extensions to the methods for checking the
#' consistency and completeness of model parameters.
#'
#' @import R6
#' @importFrom poems SimulationModel
#' @export DiseaseModel

DiseaseModel <- R6Class(
  "DiseaseModel",
  inherit = SimulationModel,
  public = list(
    ## Attributes ##

    # object_generator [inherited]

    #' @field attached A list of dynamically attached attributes (name-value pairs).
    attached = list(),

    ## Methods ##

    # Inherited methods (from GenericClass, GenericModel, SpatialModel, & SimulationModel) #
    #   new_clone(...)
    #   get_attribute_names()
    #   get_attributes(params = list(), ...)
    #   get_attribute(param)
    #   get_attribute_aliases(params = NULL)
    #   set_attributes(params = list(), ...)
    #   set_sample_attributes(params = list(), ...)
    #   is_consistent(params = NULL)
    #   inconsistent_attributes(include_nas = FALSE)
    #   is_complete()
    #   list_completeness()
    #   incomplete_attributes(include_nas = FALSE)

    # Overwritten/overridden methods #

    #' @description
    #' Initialization method sets default aliases and given attributes individually and/or from a list.
    #' @param attribute_aliases A list of alternative alias names for model attributes (form: \code{alias = "attribute"}) to be used with the set and get attributes methods.
    #' @param ... Parameters passed via a \emph{params} list or individually.
    initialize = function(attribute_aliases = NULL, ...) {
      if (!"dispersal" %in% names(c(list(...), list(...)$params))) {
        # set default alias for dispersal
        attribute_aliases <- c(attribute_aliases, list(dispersal_data = "dispersal"))
      }
      if (!"dispersal_source_n_k" %in% names(c(list(...), list(...)$params))) {
        # set default aliases for source n/k
        attribute_aliases <- c(
            attribute_aliases,
            list(
              dispersal_n_k_cutoff = "dispersal_source_n_k$cutoff",
              dispersal_n_k_threshold = "dispersal_source_n_k$threshold"
            )
          )
      }
      if (!"dispersal_target_k" %in% names(c(list(...), list(...)$params))) {
        # set default alias for target k
        attribute_aliases <- c(attribute_aliases,
            list(dispersal_k_threshold = "dispersal_target_k"))
      }
      if (!"dispersal_target_n" %in% names(c(list(...), list(...)$params))) {
        # set default aliases for target n
        attribute_aliases <- c(
            attribute_aliases,
            list(
              dispersal_n_threshold = "dispersal_target_n$threshold",
              dispersal_n_cutoff = "dispersal_target_n$cutoff"
            )
          )
      }
      if (!"dispersal_target_n_k" %in% names(c(list(...), list(...)$params))) {
        # set default aliases for target n/k
        attribute_aliases <- c(
            attribute_aliases,
            list(
              dispersal_target_n_k_threshold = "dispersal_target_n_k$threshold",
              dispersal_target_n_k_cutoff = "dispersal_target_n_k$cutoff"
            )
          )
      }
      super$initialize(attribute_aliases = attribute_aliases, ...)
    },

    #' @description
    #' Returns a boolean to indicate if (optionally selected or all) model attributes (such as dimensions) are consistent.
    #' @param params Optional array of parameter/attribute names.
    #' @return List of booleans (or NAs) to indicate consistency of selected/all attributes.
    list_consistency = function(params = NULL) {
      if (is.null(params)) { # all model attributes
        super$list_consistency()
      } else { # listed attributes
        params <- c(params)
        local_params <- params[which(params %in% c("populations", "initial_abundance", "correlation",
                                                   "carrying_capacity", "density_dependence", "dispersal",
                                                   "dispersal_target_k", "dispersal_target_n", "dispersal_target_n_k",
                                                   "dispersal_source_n_k", "result_stages", "birth",
                                                   "mortality_Sj_winter", "mortality_Sa_winter",
                                                   "mortality_Sa_summer"))]
        consistent_list <- super$list_consistency(params[which(!params %in% local_params)])
        for (param in local_params) {
          param_value <- self$get_attribute(param)
          if (is.null(param_value)) { # ignore incomplete attributes
            consistent_list[[param]] <- NA
          } else {
            consistent_list[[param]] <-
              switch(param,
                     populations = (is.numeric(param_value) && param_value > 0) &&
                       if (!is.null(self$region) && is.numeric(self$region$region_cells)) {
                         (param_value == self$region$region_cells)
                       } else TRUE,
                     initial_abundance =
                       if (is.numeric(self$populations) && is.numeric(self$stages)) {
                         if (any(class(param_value) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
                           if (!is.null(self$region) && self$region$use_raster && !is.null(self$region$region_raster)) {
                             (self$region$raster_is_consistent(param_value) &&
                                self$region$region_cells == self$populations &&
                                raster::nlayers(param_value) %in% c(1, self$stages))
                           } else {
                             NA
                           }
                         } else { # assume matrix or array
                           (is.numeric(param_value) &&
                              nrow(as.matrix(param_value)) == self$populations &&
                              ncol(as.matrix(param_value)) %in% c(1, self$stages))
                         }
                       } else {
                         NA
                       },
                     result_stages =
                       if (is.numeric(param_value) && is.numeric(self$stages)) {
                         (length(param_value) == self$stages)
                       } else {
                         NA
                       },
                     demographic_stochasticity = is.logical(param_value),
                     standard_deviation =
                       if (length(param_value) == 1) {
                         TRUE
                       } else {
                         if (is.numeric(param_value) && is.numeric(self$stages)) {
                           all(dim(as.matrix(param_value)) == self$stages)
                         } else {
                           NA
                         }
                       },
                     correlation = NA,
                     carrying_capacity =
                       if (is.numeric(self$populations) && is.numeric(self$time_steps)) {
                         if (any(class(param_value) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
                           if (!is.null(self$region) && self$region$use_raster && !is.null(self$region$region_raster)) {
                             (self$region$raster_is_consistent(param_value) &&
                                self$region$region_cells == self$populations &&
                                raster::nlayers(param_value) %in% c(1, self$time_steps))
                           } else {
                             NA
                           }
                         } else { # assume matrix or array
                           (is.numeric(param_value) &&
                              nrow(as.matrix(param_value)) == self$populations &&
                              ncol(as.matrix(param_value)) %in% c(1, self$time_steps))
                         }
                       } else {
                         NA
                       },
                     density_dependence = NA,
                     growth_rate_max = , # or
                     dispersal_source_n_k =
                       if (is.list(param_value) && all(c("cutoff", "threshold") %in% names(param_value))) {
                         consistent <- list(cutoff = NA, threshold = NA)
                         for (name in c("cutoff", "threshold")) {
                           if (length(param_value[[name]]) == 1) {
                             consistent[[name]] <- TRUE
                           } else if (length(param_value[[name]]) > 1 && is.numeric(self$populations)) {
                             consistent[[name]] <- (is.numeric(param_value[[name]]) &&
                                                      length(param_value[[name]]) == self$populations)
                           }
                         }
                         all(unlist(consistent))
                       } else {
                         NA
                       },
                     dispersal_target_k =
                       if (length(param_value) == 1) {
                         TRUE
                       } else {
                         if (is.numeric(self$populations)) {
                           (is.numeric(param_value) && length(param_value) == self$populations)
                         } else {
                           NA
                         }
                       },
                     dispersal_target_n =
                       if (is.list(param_value) && all(c("threshold", "cutoff") %in% names(param_value))) {
                         consistent <- list(threshold = NA, cutoff = NA)
                         for (name in c("threshold", "cutoff")) {
                           if (length(param_value[[name]]) == 1) {
                             consistent[[name]] <- TRUE
                           } else if (length(param_value[[name]]) > 1 && is.numeric(self$populations)) {
                             consistent[[name]] <- (is.numeric(param_value[[name]]) &&
                                                      length(param_value[[name]]) == self$populations)
                           }
                         }
                         all(unlist(consistent))
                       } else {
                         NA
                       },
                     dispersal_target_n_k =
                       if (is.list(param_value) && all(c("threshold", "cutoff") %in% names(param_value))) {
                         consistent <- list(threshold = NA, cutoff = NA)
                         for (name in c("threshold", "cutoff")) {
                           if (length(param_value[[name]]) == 1) {
                             consistent[[name]] <- TRUE
                           } else if (length(param_value[[name]]) > 1 && is.numeric(self$populations)) {
                             consistent[[name]] <- (is.numeric(param_value[[name]]) &&
                                                      length(param_value[[name]]) == self$populations)
                           }
                         }
                         all(unlist(consistent))
                       } else {
                         NA
                       },
                     translocation = NA,
                     harvest = NA,
                     mortality = NA,
                     dispersal = NA,
                     abundance_threshold =
                       if (length(param_value) == 1) {
                         TRUE
                       } else {
                         if (is.numeric(self$populations)) {
                           (is.numeric(param_value) && length(param_value) == self$populations)
                         } else {
                           NA
                         }
                       }
              )
          }
        }
        return(consistent_list)
      }
    }
  ),
  private = list(),
  active = list()
)
