ri_control <- function(control) {
  defaults<-list(fraction=.01,baseline_window=.5,min_separation=.03,
    min_relative=1e-4,min_molecular_fraction=1e-5,min_molecular_ratio=2,
    min_alkane_fraction=.25,width_quantile=.25)
  if(!is.list(control)||length(control)&&(is.null(names(control))||any(!names(control)%in%names(defaults))))
    stop("Unknown RI calibration control.",call.=FALSE)
  defaults<-utils::modifyList(defaults,control)
  if(any(lengths(defaults)!=1L)||any(!is.finite(unlist(defaults)))||
      defaults$fraction<=0||defaults$fraction>=.5||defaults$baseline_window<=0||
      defaults$min_separation<=0||defaults$min_relative<=0||defaults$min_relative>=1||
      defaults$min_molecular_fraction<=0||defaults$min_molecular_ratio<1||
      defaults$min_alkane_fraction<0||defaults$min_alkane_fraction>1||
      defaults$width_quantile<=0||defaults$width_quantile>1)
    stop("Invalid RI calibration control.",call.=FALSE)
  defaults
}

ri_crossings <- function(time,signal,apex,fraction) {
  level<-signal[apex]*fraction;n<-length(signal)
  left<-apex;right<-apex
  while(left>1L&&signal[left-1L]>level) left<-left-1L
  while(right<n&&signal[right+1L]>level) right<-right+1L
  interp<-function(a,b) if(signal[a]==signal[b]) time[a] else
    time[a]+(time[b]-time[a])*(level-signal[a])/(signal[b]-signal[a])
  start<-if(left>1L) interp(left-1L,left) else time[left]
  end<-if(right<n) interp(right,right+1L) else time[right]
  stats::setNames(as.numeric(c(start,end,end-start)),c("start","end","width"))
}

ri_input <- function(standard,time) {
  if(is.character(standard)) return(GCMS_read_cdf(standard))
  if(is.list(standard)&&!is.null(standard$Y)) {
    if(is.null(time)) time<-standard$time
    standard<-standard$Y
  }
  Y<-as.matrix(standard)
  if(is.null(time)&&!is.null(rownames(Y))) time<-suppressWarnings(as.numeric(rownames(Y)))
  mass<-suppressWarnings(as.numeric(colnames(Y)))
  if(!is.numeric(Y)||nrow(Y)<9L||ncol(Y)<3L||any(!is.finite(Y)|Y<0)||
      length(time)!=nrow(Y)||any(!is.finite(time))||any(diff(time)<=0)||
      length(mass)!=ncol(Y)||any(!is.finite(mass)|mass<0)||anyDuplicated(mass))
    stop("Supply a CDF or a nonnegative scan-by-mass matrix with numeric mass names and increasing times.",call.=FALSE)
  list(Y=Y,time=as.numeric(time),mass=mass)
}

ri_width_curve <- function(time,width,quantile) {
  # A local lower envelope estimates an intensity-independent small-peak duration.
  nbin<-min(5L,max(1L,length(time)%/%3L))
  bins<-if(nbin==1L) rep(1L,length(time)) else cut(seq_along(time),breaks=nbin,labels=FALSE)
  x<-as.numeric(tapply(time,bins,median));y<-as.numeric(tapply(width,bins,
    function(v) stats::quantile(v,quantile,names=FALSE)))
  if(length(x)<4L) return(list(time=x,width=y,spline=NULL))
  fit<-stats::smooth.spline(x,log(y),spar=.6)
  list(time=x,width=exp(stats::predict(fit,x)$y),spline=fit)
}

ri_predict_width <- function(object,time) {
  curve<-object$width_curve
  if(length(curve$time)==1L) return(rep(curve$width,length(time)))
  if(!is.null(curve$spline)) return(exp(stats::predict(curve$spline,
    pmax(min(curve$time),pmin(max(curve$time),time)))$y))
  stats::approx(curve$time,curve$width,xout=time,rule=2)$y
}

#' Calibrate retention indices from an n-alkane standard
#'
#' Identify a homologous alkane series using diagnostic ions and fit a smooth,
#' monotone retention-time calibration. The retention coordinate uses peak onset
#' plus half a local narrow-peak duration, reducing the influence of tailing.
#'
#' @param standard ANDI-MS CDF path, scan-by-mass matrix, or list with Y and time.
#' @param time Increasing scan times in minutes for matrix input.
#' @param carbon Permitted carbon numbers, default 7:40. The molecular ion must
#'   lie within the acquired mass range. Supply the known standard range when available.
#' @param anchors Optional data frame with carbon and retention_time columns,
#'   supplying verified calibration coordinates directly. An optional narrow_width
#'   column supplies local small-peak durations for timing sample profiles.
#' @param spar Smoothing-spline smoothing parameter; NULL selects it by GCV.
#' @param control Named calibration controls, described below.
#' @return An RIcalibration object with accepted anchors, all candidate diagnostics,
#'   fitted spline, local narrow-peak widths and timing convention.
#' @details For an alkane with c carbon atoms the nominal molecular ion is
#'   14c+2. Its excess above the carbon-13 isotope contribution of the adjacent
#'   alkyl ion is used with the 43, 57, 71, 85, 99 and 113 fragment series.
#'   Candidate anchors must have a coherent saturated-hydrocarbon spectrum.
#'   A maximum-weight increasing carbon-number sequence removes duplicate
#'   maxima and inconsistent assignments. Solvent-range C6 is excluded by default.
#'
#'   Peak onset and duration use a relative-height fraction of 0.01 by default.
#'   A local lower envelope of standard peak durations estimates the duration
#'   of a small peak. Both standards and sample profiles use
#'   t_star = onset + narrow_width/2. A smoothing spline maps t_star to 100c;
#'   a monotone cubic interpolation through its fitted anchor values prevents
#'   local reversals. Extrapolation is disabled by default in predict().
#'
#'   control accepts fraction, baseline_window (minutes), min_separation
#'   (minutes), min_relative, min_molecular_fraction, min_molecular_ratio,
#'   min_alkane_fraction and width_quantile. The returned anchor table includes
#'   the molecular-ion evidence and calibration residuals for inspection.
#'   Explicit anchors bypass automatic assignment and use their supplied timing.
#' @references van den Dool, H. and Kratz, P. D. (1963), J. Chromatogr. 11,
#'   463--471, doi:10.1016/S0021-9673(01)80947-X.
#'   <https://webbook.nist.gov/chemistry/gc-ri/>
#' @examples
#' a <- data.frame(carbon=8:15, retention_time=c(3,4,5.2,6.7,8.5,10.6,13,15.7))
#' cal <- calibrate_ri(anchors=a)
#' predict(cal, c(5,10))
#' cal
#' @export
calibrate_ri <- function(standard=NULL,time=NULL,carbon=7:40,anchors=NULL,
                         spar=NULL,control=list()) {
  begun<-proc.time()[["elapsed"]];control<-ri_control(control)
  if(!is.numeric(carbon)||anyNA(carbon)||any(carbon<3|carbon!=floor(carbon))||anyDuplicated(carbon))
    stop("carbon must contain distinct integer carbon numbers.",call.=FALSE)
  candidates<-data.frame();width_curve<-NULL;source<-if(is.character(standard)) standard else NULL
  if(is.null(anchors)) {
    raw<-ri_input(standard,time);Y<-raw$Y;time<-raw$time;mass<-raw$mass
    dt<-median(diff(time));carbon<-sort(carbon[14*carbon+2<=max(mass)&14*carbon+2>=min(mass)])
    if(length(carbon)<4L) stop("At least four alkane molecular ions must be in the acquired mass range.")
    ions<-which(mass%in%c(43,57,71,85,99,113))
    if(length(ions)<3L) stop("The standard lacks alkane diagnostic ion channels.")
    signal<-rowSums(Y[,ions,drop=FALSE])
    span<-max(3L,as.integer(round(control$baseline_window/dt))%/%2L*2L+1L)
    span<-min(span,as.integer(nrow(Y)-1L+(nrow(Y)%%2L)))
    baseline<-stats::runmed(signal,span,endrule="median")
    residual<-pmax(signal-baseline,0)
    smooth<-as.numeric(stats::filter(residual,c(.25,.5,.25),sides=2));smooth[c(1,length(smooth))]<-0
    peaks<-which(diff(sign(diff(smooth)))<0)+1L
    peaks<-peaks[smooth[peaks]>=max(smooth)*control$min_relative]
    keep<-integer()
    for(i in peaks[order(-smooth[peaks])])
      if(!length(keep)||all(abs(time[i]-time[keep])>=control$min_separation)) keep<-c(keep,i)
    keep<-sort(keep)
    rows<-lapply(keep,function(i) {
      radius<-max(8L,as.integer(round(.08/dt)))
      at<-seq.int(max(1L,i-radius),min(nrow(Y),i+radius))
      edge<-unique(c(head(at,3L),tail(at,3L)))
      bg<-apply(Y[edge,,drop=FALSE],2L,median)
      peak_rows<-seq.int(max(1L,i-2L),min(nrow(Y),i+2L))
      spectrum<-colSums(pmax(sweep(Y[peak_rows,,drop=FALSE],2L,bg,"-"),0))
      get<-function(m) {z<-spectrum[match(m,mass)];z[is.na(z)]<-0;z}
      excess<-pmax(0,get(14*carbon+2)-.0107*carbon*get(14*carbon+1))
      o<-order(excess,decreasing=TRUE);j<-o[1L]
      ratio<-excess[j]/max(1,excess[o[min(2,length(o))]])
      total<-sum(spectrum);fragment<-sum(spectrum[ions])/max(1,total)
      geometry<-ri_crossings(time,residual,i,control$fraction)
      data.frame(scan=i,apex=time[i],carbon=carbon[j],molecular_mass=14*carbon[j]+2,
        molecular_fraction=excess[j]/max(1,total),molecular_ratio=ratio,
        alkane_fraction=fragment,signal=smooth[i],start=geometry['start'],
        end=geometry['end'],width=geometry['width'],weight=log1p(excess[j]),row.names=NULL)
    })
    if(!length(rows)) stop("No alkane peaks found; check the standard or supply anchors.")
    candidates<-do.call(rbind,rows)
    good<-with(candidates,molecular_fraction>=control$min_molecular_fraction &
      molecular_ratio>=control$min_molecular_ratio & alkane_fraction>=control$min_alkane_fraction & width>dt)
    eligible<-which(good);best<-numeric(length(eligible));parent<-integer(length(eligible))
    for(k in seq_along(eligible)) {
      prior<-which(seq_along(eligible)<k & candidates$carbon[eligible]<candidates$carbon[eligible[k]])
      if(length(prior)) {a<-prior[which.max(best[prior])];parent[k]<-a;best[k]<-best[a]}
      best[k]<-best[k]+candidates$weight[eligible[k]]
    }
    chosen<-integer()
    if(length(best)) {k<-which.max(best);while(k>0L){chosen<-c(eligible[k],chosen);k<-parent[k]}}
    candidates$accepted<-seq_len(nrow(candidates))%in%chosen
    anchors<-candidates[chosen,,drop=FALSE]
    if(nrow(anchors)<4L) stop("Fewer than four coherent alkane assignments; inspect diagnostics or supply anchors.")
    width_curve<-ri_width_curve(anchors$apex,anchors$width,control$width_quantile)
    temp<-list(width_curve=width_curve)
    anchors$narrow_width<-ri_predict_width(temp,anchors$start)
    anchors$retention_time<-anchors$start+anchors$narrow_width/2
  } else {
    if(!is.data.frame(anchors)||!all(c('carbon','retention_time')%in%names(anchors)))
      stop("anchors requires carbon and retention_time columns.")
    anchors<-anchors[order(anchors$retention_time),,drop=FALSE]
    if('narrow_width'%in%names(anchors)) {
      if(any(!is.finite(anchors$narrow_width)|anchors$narrow_width<=0))
        stop("Explicit narrow_width values must be positive.")
      width_curve<-ri_width_curve(anchors$retention_time,anchors$narrow_width,control$width_quantile)
    }
  }
  if(nrow(anchors)<4L||any(!is.finite(anchors$carbon))||any(anchors$carbon!=floor(anchors$carbon))||
      any(!is.finite(anchors$retention_time))||any(diff(anchors$retention_time)<=0)||any(diff(anchors$carbon)<=0))
    stop("Use at least four increasing alkane times and carbon numbers.")
  spline<-stats::smooth.spline(anchors$retention_time,100*anchors$carbon,spar=spar)
  fitted<-stats::predict(spline,anchors$retention_time)$y
  if(any(diff(fitted)<=0)) stop("The smoothed calibration is not increasing; inspect assignments or adjust spar.")
  anchors$RI<-100*anchors$carbon;anchors$fitted_RI<-fitted;anchors$residual<-anchors$RI-fitted
  structure(list(anchors=anchors,candidates=candidates,spline=spline,width_curve=width_curve,
    fraction=control$fraction,source=source,control=control,
    timing=if(is.null(width_curve)) 'provided' else 'onset_plus_half_narrow_width',
    seconds=proc.time()[['elapsed']]-begun),class='RIcalibration')
}

#' @rdname calibrate_ri
#' @param object A fitted RIcalibration.
#' @param newdata Retention coordinates in minutes.
#' @param extrapolate Use linear continuation outside the calibration range.
#'   FALSE returns NA outside the range.
#' @param ... Additional method arguments.
#' @method predict RIcalibration
#' @export
predict.RIcalibration <- function(object,newdata,extrapolate=FALSE,...) {
  x<-object$anchors$retention_time;y<-object$anchors$fitted_RI
  value<-stats::splinefun(x,y,method='hyman')(newdata)
  left<-which(!is.na(newdata)&newdata<min(x));right<-which(!is.na(newdata)&newdata>max(x))
  if(extrapolate) {
    value[left]<-y[1]+(newdata[left]-x[1])*diff(y[1:2])/diff(x[1:2])
    n<-length(x);value[right]<-y[n]+(newdata[right]-x[n])*diff(tail(y,2))/diff(tail(x,2))
  } else value[c(left,right)]<-NA_real_
  value
}

#' @rdname calibrate_ri
#' @param x An RIcalibration object.
#' @method print RIcalibration
#' @export
print.RIcalibration <- function(x,...) {
  cat(nrow(x$anchors),' alkane anchors; C',min(x$anchors$carbon),'--C',max(x$anchors$carbon),
      '; calibration RMSE ',format(sqrt(mean(x$anchors$residual^2)),digits=3),' RI units\n',sep='')
  invisible(x)
}
