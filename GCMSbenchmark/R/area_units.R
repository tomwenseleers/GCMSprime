# Convert native factor/ion units using only exported factors and acquisition
# times. No matching, reference spectra or simulated areas enter this step.
tic_units <- function(fit, method, files) {
  if(isTRUE(fit$tic_units))return(fit)
  clocks<-stats::setNames(lapply(files,cdf_clock),sample_id(files))
  fit$native_area_kind<-fit$area_kind
  fit$native_profile_kind<-fit$profile_kind
  fit$peaks$native_area<-fit$peaks$area
  fit$peaks$native_height<-fit$peaks$height
  fit$peaks$area_conversion<-rep(NA_real_,nrow(fit$peaks))
  for(j in seq_len(nrow(fit$peaks))) {
    clock<-clocks[[fit$peaks$sample[j]]]
    delta<-stats::median(diff(clock))*60
    if(!is.finite(delta)||delta<=0)stop("Acquisition times must increase.")
    spectrum<-fit$spectra[[j]]
    total<-if(length(spectrum))sum(spectrum[,2L])else 0
    profile<-fit$profiles[[j]]
    profile_factor<-1
    area_factor<-switch(method,
      GCMSprime=1,AMDIS=1,"xcms+CAMERA"=1/delta,MSHub=1/delta,
      "MS-DIAL"={
        profile_factor<-if(total>0&&fit$peaks$height[j]>0)
          total/fit$peaks$height[j]else NA_real_
        profile_factor/delta
      },
      eRah={
        # eRah exports round(1000*S); its Profile is C/max(C), restored
        # to C by the adapter. C %*% S is the original ion reconstruction.
        profile_factor<-total/1000
        step<-if(length(profile)&&nrow(profile)>1L)
          stats::median(diff(profile[,1L]))/stats::median(diff(clock))else 1
        profile_factor*step
      },
      mzmine={
        if(is.null(fit$peaks$quant_mass))stop("mzmine conversion requires the representative m/z.")
        q<-match(round(fit$peaks$quant_mass[j]),round(spectrum[,1L]))
        fraction<-if(!is.na(q)&&total>0)spectrum[q,2L]/total else NA_real_
        # mzmine integrates intensity over retention time in minutes.
        60/delta/fraction
      },stop("Unknown area convention: ",method))
    if(!is.finite(area_factor)||area_factor<=0)
      stop("Cannot determine a total-ion area conversion for ",method," component ",j,".")
    fit$peaks$area_conversion[j]<-area_factor
    fit$peaks$area[j]<-fit$peaks$area[j]*area_factor
    fit$peaks$height[j]<-fit$peaks$height[j]*profile_factor
    if(length(profile)) {
      if(method=="xcms+CAMERA") {
        count<-if(nrow(profile)>1L)sum(stats::approx(profile[,1L],profile[,2L],clock,
          yleft=0,yright=0,ties=mean)$y)else sum(profile[,2L])
        # CAMERA supplies a representative-ion EIC, not a shared-component
        # profile. Scale this explicitly labelled proxy to the grouped TIC area.
        profile_factor<-if(count>0)fit$peaks$area[j]/count else 1
      }
      profile[,2L]<-profile[,2L]*profile_factor
      fit$profiles[[j]]<-profile
    }
    if(total>0) spectrum[,2L]<-spectrum[,2L]/total
    fit$spectra[[j]]<-spectrum
  }
  fit$area_kind<-"integrated counts"
  fit$profile_kind<-if(method=="xcms+CAMERA")"TIC-scaled representative-ion trace"else
    if(method%in%c("mzmine","MSHub"))"unavailable"else "component TIC"
  fit$tic_units<-TRUE
  fit
}
