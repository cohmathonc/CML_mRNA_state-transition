
library(survival)

colData<- read.csv("C:\\Users\\uechi\\Downloads\\cml_sampleInfo_dataTable.csv")


# "colData(dat.se)" is the .csv table you have so you should be able to replace that with your table and use this directly
time <- as.numeric(difftime(colData$dod, colData$dob, units = "days"))
status <- ifelse(is.na(colData$dod), 0, 1)
sample_data <-data.frame()
sample_data <- data.frame("time" = time)
sample_data[["status"]] <- status
sample_data[["group"]] <- colData$Group
sample_data[["mouse"]] <- colData$mouse_id

### all mice ### 
# Fit a survival model
surv_object <- with(sample_data, Surv(time, status))

# Create a Kaplan-Meier survival curve
surv_curve <- survfit(surv_object ~ 1)
plot(surv_curve)

### CML mice (group "B") only ###
#b_data <- sample_data[which(sample_data$group=="B"),]
#surv_b <- with(b_data, Surv(time, status))
#b_curve <- survfit(surv_b ~ 1)
#plot(b_curve)

# Create a survival object
surv_object <- Surv(time, status)

# Print the survival object
print(surv_object)