dt1 <- read.csv("Data/all.data.mammals.final.final.csv") %>% as.data.table()
mammal_observation_table <- read.csv("Data/mammals_observation_table.csv") %>% as.data.table()
dt2 <- mammal_observation_table[,.SD,.SDcols=names(dt1)]

key_cols <- c("Project.ID","Unit","Deployment.id","Time_new","Species","Start.date","Timestamp") 
setkeyv(dt1,key_cols)
setkeyv(dt2,key_cols)

rows_only_in_dt1 <- dt1[!dt2]
rows_only_in_dt2 <- dt2[!dt1]
all_differing_rows <- rbindlist(list(rows_only_in_dt1, rows_only_in_dt2))

print(paste("number of rows only in dt1:",nrow(rows_only_in_dt1)))
print(paste("number of rows only in dt2:",nrow(rows_only_in_dt2)))

obs_table <- dt1[dt2, nomatch = 0]

# Select only the columns from dt1
obs_table <- obs_table[, names(dt1), with = FALSE]

# If you want to ensure you have unique rows:
obs_table <- unique(obs_table)

A <- copy(obs_table)
A[,idx:=.I]
A[,grp:=.GRP,by = key(dt1)]

# using the grp and idx indices I noticed that there are multiple observations of the same species at the same deployment and timestamp.
# Currently leave them as is.