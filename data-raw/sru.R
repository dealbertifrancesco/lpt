# data-raw/sru.R
# One-time transformation: trim sru to 13 periods and re-index.
# Run from package root: Rscript data-raw/sru.R

load("data/sru.rda")

sru <- sru[sru$year %in% c(1993:1999, 2014:2019), ]

year_map <- c(setNames(-7:-1, 1993:1999), setNames(0:5, 2014:2019))
sru$year <- as.integer(year_map[as.character(sru$year)])

save(sru, file = "data/sru.rda", compress = "xz")
cat(sprintf("Saved sru: %d rows, %d communes, periods %d to %d\n",
            nrow(sru), length(unique(sru$commune)),
            min(sru$year), max(sru$year)))
