#' This example shows how a user can use the impute.subject() function to impute
#' the visits of a single patient by using the data from another clinical
#' register.

data(patient.data)
data(new.patient)
#' The user must define which features are static/dynamic and
#' continuous/categorical/ordinal.
static.features = c(
  "sex",
  "bmi_premorbid",
  "bmi_diagnosis",
  "fvc_diagnosis",
  "familiality",
  "genetics",
  "ftd",
  "onset_site",
  "onset_age"
)
dynamic.features = c(
  "niv",
  "peg",
  "alsfrs_1",
  "alsfrs_2",
  "alsfrs_3",
  "alsfrs_4",
  "alsfrs_5",
  "alsfrs_6",
  "alsfrs_7",
  "alsfrs_8",
  "alsfrs_9",
  "alsfrs_10",
  "alsfrs_11",
  "alsfrs_12"
)
continuous.features = c("bmi_premorbid",
                        "bmi_diagnosis",
                        "fvc_diagnosis",
                        "onset_age")
categorical.features = c("sex",
                         "familiality",
                         "genetics",
                         "ftd",
                         "onset_site",
                         "niv",
                         "peg")
ordinal.features = c(
  "alsfrs_1",
  "alsfrs_2",
  "alsfrs_3",
  "alsfrs_4",
  "alsfrs_5",
  "alsfrs_6",
  "alsfrs_7",
  "alsfrs_8",
  "alsfrs_9",
  "alsfrs_10",
  "alsfrs_11",
  "alsfrs_12"
)

#' In what follows, the impute.subject() function is used to impute the missing
#' values in the visits of a new patient in a 3 months wide time window.
#' Please note that missing values in the visits outside of this window will not
#' be imputed.
imputed.patient.data <-
  impute.subject(
    subject.to.impute = new.patient,
    # data frame containing two visits with missing data to be imputed
    candidates = patient.data,
    # dataset of patients to be used as candiates for the wkNNMI algorithm
    window_size = 3,
    # how many months of patient data to impute
    K = 5,
    # number of neighbours to consider for the imputation
    static.features = static.features,
    dynamic.features = dynamic.features,
    continuous.features = continuous.features,
    categorical.features = categorical.features,
    ordinal.features = ordinal.features,
    time.feature = "visit_time",
    # the time feature
    sub.id.feature = "subID"
  )
