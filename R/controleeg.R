#' Visual-task EEG Dataset - Healthy Controls
#'
#' This dataset contains list of data frames of 45-second EEG recordings (at sampling rate of 128Hz) from 19-channels of 5 healthy controls performing a visual task of counting characters on a flashed image. Channels A1 and A2 located at the earlobes are used for referencing.
#'
#' @format Each data frame contains the following:
#' \describe{
#'   \item{Fp1}{Numeric. Recordings from Channel Fp1.}
#'   \item{Fp2}{Numeric. Recordings from Channel Fp2.}
#'   \item{F3}{Numeric. Recordings from Channel F3.}
#'   \item{F4}{Numeric. Recordings from Channel F4.}
#'   \item{C3}{Numeric. Recordings from Channel C3.}
#'   \item{C4}{Numeric. Recordings from Channel C4.}
#'   \item{P3}{Numeric. Recordings from Channel P3.}
#'   \item{P4}{Numeric. Recordings from Channel P4.}
#'   \item{O1}{Numeric. Recordings from Channel O1.}
#'   \item{O2}{Numeric. Recordings from Channel O2.}
#'   \item{F7}{Numeric. Recordings from Channel F7.}
#'   \item{F8}{Numeric. Recordings from Channel F8.}
#'   \item{T7}{Numeric. Recordings from Channel T7.}
#'   \item{T8}{Numeric. Recordings from Channel T8.}
#'   \item{P7}{Numeric. Recordings from Channel P7.}
#'   \item{P8}{Numeric. Recordings from Channel P8.}
#'   \item{Fz}{Numeric. Recordings from Channel Fz.}
#'   \item{Cz}{Numeric. Recordings from Channel Cz.}
#'   \item{Pz}{Numeric. Recordings from Channel Pz.}
#' }
#'
#' @examples
#' # Recordings from the heathly controls
#' head(controleeg[[1]])
#' head(controleeg[[2]])
#' head(controleeg[[3]])
#' head(controleeg[[4]])
#' head(controleeg[[4]])
#'
#' @source Generated for demonstration purposes.
"controleeg"
