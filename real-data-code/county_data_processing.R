library(dplyr)


##### Data pre-processing #####

# Load data and remove non-useful features
crime_df <- read.csv('data/crimedata.csv') %>%
  dplyr::select(
    -c(
      # "X.communityname",
      "countyCode",
      "communityCode",
      "LemasSwornFT",
      "LemasSwFTPerPop",
      "LemasSwFTFieldOps",
      "LemasSwFTFieldPerPop",
      "LemasTotalReq",
      "LemasTotReqPerPop",
      "PolicReqPerOffic",
      "PolicPerPop",
      "RacialMatchCommPol",
      "PctPolicWhite",
      "PctPolicBlack",
      "PctPolicHisp",
      OwnOccQrange,
      RentQrange,
      PctPolicAsian,
      PctPolicMinor,
      OfficAssgnDrugUnits,
      NumKindsDrugsSeiz,
      PolicAveOTWorked,
      PolicCars,
      PolicOperBudg,
      LemasPctPolicOnPatr,
      LemasGangUnitDeploy,
      PolicBudgPerPop,
      fold,
      "murders",
      "murdPerPop",
      population,
      "rapes"   ,
      "rapesPerPop"    ,
      "robberies",
      "robbbPerPop"  ,
      "assaults"        ,
      "assaultPerPop",
      "burglaries"  ,
      "burglPerPop" ,
      "larcenies" ,
      "larcPerPop" ,
      "autoTheft",
      "autoTheftPerPop" ,
      "arsons"           ,
      "arsonsPerPop",
      "nonViolPerPop"
    )
  )

# Load county information
county <- read.csv('data/county.csv')
colnames(county) <- c("X.communityname", "County", "state")

# Merge county information with crime data
crime_df <- crime_df %>% left_join(county, by = c("X.communityname","state"))

# Remove missing values
crime_df[crime_df=='?'] <- NA
crime_df <- crime_df %>% na.omit()

# Format data
crime_df$OtherPerCap <- as.numeric(crime_df$OtherPerCap)
crime_df$ViolentCrimesPerPop <-
  as.numeric(crime_df$ViolentCrimesPerPop)

# Scale data
crime_df <-
  crime_df %>% mutate_at(.vars = setdiff(colnames(crime_df), c("X.communityname","state","County")), 
                         .funs = scale)
crime_df$Intercept = 1

# Select counties with more than 3 observations
crime_df$countystate <- paste(crime_df$County, crime_df$state, sep = "_")

counties <- unique(crime_df$countystate)
county_list <- state_list <- NULL
for (i in 1:length(counties)) {
  if (sum(crime_df$countystate == counties[i]) > 3) {
    county_list <- c(county_list, counties[i])
    state_list <- c(state_list, crime_df$state[crime_df$countystate == counties[i]][1])
  }
}

crime_df <- crime_df %>% select(-c(County, state, X.communityname))















