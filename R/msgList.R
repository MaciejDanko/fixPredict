#' @keywords internal
#' @noRd
msgList <- function(x) {
  tmp <-c('offset() within a formula is not supported yet. Try to fit a model using the gamm4 function and "offset" parameter.', #1
          'Model matrix (MM) not supported for the predFix function.', #2
          'Averaging over missig variables has been mot implemented for mgcv package yet.', #3
          'getMM has been mot implemented for mgcv package yet.', #4
          'Only one offset can be assigned.', #5
          'offset() within a formula is not supported yet. Try to fit a model using the gamm4 function and "offset" parameter.', #6
          'Model matrix (MM) not supported for the predFix function.', #7
          'Averaging over missig variables has been mot implemented for gamm4 package yet.', #8
          'getMM has been mot implemented for gamm4 package yet.', #9
          'offset() within a formula is not supported yet. Try to fit a model using the lme4 function and "offset" parameter.', #10
          'Model matrix (MM) is missing for the predFix function.', #11
          'offset() within a formula is not supported yet. Try to use "offset" parameter.', #12
          'Unknown / unimplemented mixed model class.', #13
          'Random effects and/or response variables dropped from the newdata.', #14
          'Unknown level(s) found in the newdata.', #15
          'Unknown method.', #16
          'Should not happen.', #17
          '"as.rate = TRUE" only makes sense for the response scale (set type = "response").', #18
          'Unknown type.', #19
          'Not supported formula in the object.', #20
          'Missing variable(s) in newdata: ',#21
          '.\n Add missing variables or set average.missing = TRUE.', #22
          '"as.rate = TRUE" only makes sense when link function is "log".', #23 
          'The fitted random effect formula is currently not supported in the bootstrap computations.', #24
          'Bootstrapping random effects is supported only for merMod objects (lme4 package).', #25
          'The conditional variance cannot be retrived from the model.', #26
          'Bootstrapping random effects is currently not supported for gam objects (gamm4 package).' #27
          )
  tmp[x]
}
