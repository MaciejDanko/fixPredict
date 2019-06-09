#' @keywords internal
#' @noRd
msgList <- function(x) {
  tmp <-c('offset() within a formula is not supported yet. Try to fit a model using the gamm4 function and "offset" parameter.',
          'Model matrix (MM) not supported for the predFix function.',
          'Averaging over missig variables has been mot implemented for mgcv package yet.',
          'getMM has been mot implemented for mgcv package yet.',
          'Only one offset can be assigned.',
          'offset() within a formula is not supported yet. Try to fit a model using the gamm4 function and "offset" parameter.',
          'Model matrix (MM) not supported for the predFix function.',
          'Averaging over missig variables has been mot implemented for gamm4 package yet.',
          'getMM has been mot implemented for gamm4 package yet.',
          'offset() within a formula is not supported yet. Try to fit a model using the lme4 function and "offset" parameter.',
          'Model matrix (MM) is missing for the predFix function.',
          'offset() within a formula is not supported yet. Try to use "offset" parameter.',
          'Unknown / unimplemented mixed model class.',
          'Random effects and/or response variables dropped from the newdata.',
          'Unknown level(s) found in the newdata.',
          'Unknown method.',
          'Should not happen.',
          '"as.rate = TRUE" only makes sense for the response scale (set type = "response").',
          'Unknown type.',
          'Not supported formula in the object.',
          'Missing variable(s) in newdata: ',
          '.\n Add missing variables or set average.missing = TRUE.',
          '"as.rate = TRUE" only makes sense when link function is "log".' 
          )
  tmp[x]
}
