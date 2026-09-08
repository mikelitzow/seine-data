# QA/QC for the age-0 cod / pollock seine CPUE time series
#
# Focused on the current update year, but every structural check also runs
# across all years so problems that slipped in previously get caught.
#
# Run it after add.cpue.R:  Rscript scripts/qaqc_seine_cpue.R
# Nothing is modified -- this only reports.
#
# Background: the final bay/site columns in the time series come from the SITE
# file, not the cpue file, so a typo in the site file lands in the series
# unnoticed. That is how "Japanese " (trailing space) entered in 2024 and
# "Kiul-4" entered in 2025. Checks 12-14 exist to catch that class of error.

suppressMessages(library(tidyverse))

YEAR <- 2026

# ---------------------------------------------------------------- reporting --
results <- list()
check <- function(label, ok, detail = NULL, warn_only = FALSE) {
  status <- if (ok) "PASS" else if (warn_only) "WARN" else "FAIL"
  results[[length(results) + 1]] <<- data.frame(check = label, status = status)
  cat(sprintf("[%-4s] %s\n", status, label))
  if (!ok && !is.null(detail)) {
    out <- capture.output(print(detail))
    cat(paste0("        ", out, collapse = "\n"), "\n")
  }
  invisible(ok)
}
section <- function(x) cat("\n", strrep("-", 70), "\n", x, "\n", strrep("-", 70), "\n", sep = "")

# -------------------------------------------------------------------- input --
dat <- read.csv("./data/age.0_cod_pollock_seine_cpue.csv") %>% select(-1)

# File naming and column layout are not consistent from year to year, so resolve
# each file by trying the known patterns and then normalise the columns the
# checks rely on. 2021 in particular uses cpue_2021.csv (underscore), carries
# month/day/year instead of a date column, and capitalises Species/Length.
pick <- function(...) {
  cand <- c(...)
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

f_cpue <- pick(sprintf("./data/cpue%d.csv", YEAR),
               sprintf("./data/cpue_%d.csv", YEAR))
f_len  <- pick(sprintf("./data/length%d.csv", YEAR),
               sprintf("./data/length_%d.csv", YEAR))
f_site <- pick(sprintf("./data/site%d.csv", YEAR),
               sprintf("./data/site_%d.csv", YEAR))
f_kod  <- pick(sprintf("./data/Kodiak_%d_seine_data_gadid.csv", YEAR),
               sprintf("./data/Kodiak_%d_seine_data_gadid_and_pinks.csv", YEAR),
               sprintf("./data/Kodiak %d seine data - gadid and salmonid.csv", YEAR))

missing_files <- c(cpue = f_cpue, length = f_len, site = f_site, Kodiak = f_kod)
if (any(is.na(missing_files))) {
  cat(sprintf("\nCannot run %d-specific checks -- no file found for: %s\n",
              YEAR, paste(names(missing_files)[is.na(missing_files)], collapse = ", ")))
  cat("The all-years structural checks below still apply.\n")
  YEAR_OK <- FALSE
} else {
  YEAR_OK <- TRUE
}

read_norm <- function(path) {
  x <- read.csv(path)
  blank <- is.na(names(x)) | names(x) == ""
  names(x)[blank] <- paste0("V", seq_len(ncol(x)))[blank]
  # 2021 capitalises these in the length file
  names(x)[names(x) == "Species"] <- "species"
  names(x)[names(x) == "Length"]  <- "length"
  # build a date column where the file only carries month/day/year
  if (!"date" %in% names(x) && !"Date" %in% names(x) &&
      all(c("month", "day", "year") %in% names(x))) {
    x$Date <- paste(x$month, x$day, x$year, sep = "/")
  }
  if (!"Date" %in% names(x) && "date" %in% names(x)) x$Date <- x$date
  if (!"date" %in% names(x) && "Date" %in% names(x)) x$date <- x$Date
  x
}

if (YEAR_OK) {
  cp <- read_norm(f_cpue)
  ln <- read_norm(f_len)
  st <- read_norm(f_site)
  kd <- read_norm(f_kod)

  # blank trailing rows are expected in the raw files; drop them for checking
  cp_r <- cp %>% filter(species != "", !is.na(Station))
  ln_r <- ln %>% filter(species != "", !is.na(Station))
  st_r <- st %>% filter(!is.na(Station))
}

# Mirror the per-year corrections add.cpue.R applies to the raw files, so this
# script validates the pipeline as it actually runs rather than re-reporting
# problems that have already been handled. Keep in step with add.cpue.R.
if (YEAR_OK && YEAR == 2023) {
  # station 531 was used for both Jap-1 (7/30) and Jap-5 (7/31); Jap-5 is 529
  cp_r$Station[cp_r$Station == 531 & cp_r$Site == "Jap-5"] <- 529
}
if (YEAR_OK && YEAR == 2024) {
  # station 589 (Fox-1, 12:10) was a bad tow with no catch record
  st_r$use.for.CPUE[st_r$Station == 589] <- "no"
}

if (YEAR_OK) {

section(sprintf("RAW %d FILES", YEAR))

check(sprintf("cpue%d.csv has rows after dropping blanks", YEAR), nrow(cp_r) > 0)
check(sprintf("length%d.csv has rows after dropping blanks", YEAR), nrow(ln_r) > 0)
check(sprintf("site%d.csv has rows after dropping blanks", YEAR), nrow(st_r) > 0)
cat(sprintf("        blank rows dropped -- cpue: %d, length: %d, site: %d\n",
            nrow(cp) - nrow(cp_r), nrow(ln) - nrow(ln_r), nrow(st) - nrow(st_r)))

# 1. dates must all parse
for (nm in c("cpue", "site", "length", "Kodiak")) {
  v <- switch(nm,
              cpue   = cp_r$date,
              site   = st_r$Date,
              length = ln_r$Date,
              Kodiak = kd$Date)
  bad <- sum(is.na(as.Date(v, format = "%m/%d/%Y")))
  check(sprintf("all dates parse in %s file", nm), bad == 0, sprintf("%d unparseable", bad))
}

# 2. cpue month/day/year columns must agree with the date column
d_a <- as.Date(cp_r$date, format = "%m/%d/%Y")
d_b <- as.Date(paste(cp_r$year, cp_r$month, cp_r$day, sep = "-"))
check("cpue month/day/year agrees with date column", all(d_a == d_b, na.rm = TRUE),
      cp_r[which(d_a != d_b), c("date", "month", "day", "year")])

# 3/4. add.cpue.R drops May samples (julian <= 151) and Chief Cove from every
#      year, so their presence in a raw file is expected, not an error -- 2022
#      legitimately contains both. Report what the pipeline will drop; the
#      all-years section then asserts none of it reached the series.
st_r$julian <- lubridate::yday(as.Date(st_r$Date, format = "%m/%d/%Y"))
n_may <- sum(st_r$julian <= 151, na.rm = TRUE)
n_chief <- sum(trimws(st_r$Bay) == "Chief Cove", na.rm = TRUE)
cat(sprintf("        pipeline will drop -- May samples: %d, Chief Cove rows: %d\n",
            n_may, n_chief))

# 5. use.for.CPUE is populated and valid
check("use.for.CPUE is yes/no for every row", all(st_r$use.for.CPUE %in% c("yes", "no")),
      table(st_r$use.for.CPUE, useNA = "ifany"))

section(sprintf("%d CROSS-FILE CONSISTENCY", YEAR))

st_yes <- st_r %>% filter(use.for.CPUE == "yes")
cp_r$julian <- lubridate::yday(as.Date(cp_r$date, format = "%m/%d/%Y"))

# 6. every retained station has catch records
miss <- setdiff(st_yes$Station, cp_r$Station)
check("every use.for.CPUE=yes station appears in the cpue file", length(miss) == 0,
      st_yes %>% filter(Station %in% miss) %>% select(Station, Site, Bay))

# 7. One station number pointing at two different SITES is the damaging case --
#    it can attach a catch to the wrong place. 2023 used 531 for both Jap-1 and
#    Jap-5; the Jap-5 rows are reassigned to 529 in add.cpue.R.
collide <- cp_r %>%
  distinct(Station, Site) %>%
  count(Station) %>%
  filter(n > 1)
check("each station number maps to a single site in the cpue file",
      nrow(collide) == 0,
      cp_r %>% filter(Station %in% collide$Station) %>%
        distinct(Station, date, Site, Bay) %>% arrange(Station, date))

#    A station spanning two dates is milder -- station numbers are otherwise
#    one-per-set, so it is worth seeing but is not necessarily wrong.
multi_date <- cp_r %>%
  distinct(Station, date) %>%
  count(Station) %>%
  filter(n > 1)
check("each station number appears on a single date in the cpue file",
      nrow(multi_date) == 0,
      cp_r %>% filter(Station %in% multi_date$Station) %>%
        distinct(Station, date, Site, Bay) %>% arrange(Station, date),
      warn_only = TRUE)

# 8. The site<->cpue join is on Station + year + julian, so a station whose
#    gadid records sit at a different julian than its site row would silently
#    become a zero. A station with no gadid records at all is a genuine zero and
#    is not flagged here -- check 6 covers a station missing from the file
#    entirely.
gad <- cp_r %>% filter(species %in% c("Pacific cod", "walleye pollock"))
joinable <- st_yes %>%
  select(Station, Site, Bay, julian_site = julian) %>%
  left_join(gad %>% distinct(Station, julian_cpue = julian), by = "Station") %>%
  group_by(Station, Site, Bay, julian_site) %>%
  summarize(has_gadid = any(!is.na(julian_cpue)),
            joins = any(julian_cpue == julian_site, na.rm = TRUE),
            .groups = "drop")
lost <- joinable %>% filter(has_gadid, !joins)
check("every station with gadid records joins on Station + year + julian",
      nrow(lost) == 0, as.data.frame(lost))

# 8. site NAMES agree across the three files per station.
#    Reported as a warning: the age-1 join is keyed on Station, so a mismatch is
#    no longer silently wrong, but disagreeing vocabularies are still a defect.
nm3 <- st_yes %>%
  distinct(Station, site_file = Site) %>%
  left_join(cp_r %>% distinct(Station, cpue_file = Site), by = "Station") %>%
  left_join(ln_r %>% distinct(Station, length_file = Site), by = "Station") %>%
  filter(site_file != cpue_file | (!is.na(length_file) & site_file != length_file))
check("site names agree across site/cpue/length files", nrow(nm3) == 0, nm3, warn_only = TRUE)

# 9. no duplicate Station x species in cpue (would make pivot_wider build list-columns)
dupsp <- cp_r %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  count(Station, species) %>% filter(n > 1)
check("no duplicate Station x species rows in cpue file", nrow(dupsp) == 0, dupsp)

section(sprintf("%d AGE-1 CUTOFFS", YEAR))

# 10. the cutoff should sit in an empty part of the length distribution.
#     If fish appear near the cut the split is ambiguous and needs a human eye.
cutoff_gap <- function(sp, cut, halo = 20) {
  L <- ln_r %>% filter(species == sp) %>% pull(length)
  near <- sum(L > (cut - halo) & L <= (cut + halo))
  list(n = length(L), near = near,
       tab = table(cut(L, breaks = seq(0, 300, 10))))
}
for (s in list(list("Pacific cod", 120), list("walleye pollock", 100))) {
  r <- cutoff_gap(s[[1]], s[[2]])
  check(sprintf("%s: no fish within 20mm of the %dmm cutoff (n = %d)", s[[1]], s[[2]], r$n),
        r$near == 0, r$tab, warn_only = TRUE)
}

section(sprintf("%d EXTREME CATCHES", YEAR))

# 10b. Very large catches are normally extrapolated from a subsample rather than
#      counted whole, so an extreme value with subsample = "yes" is expected
#      methodology. An extreme value with subsample = "no" is worth a look.
#      Only stations that reach the series are considered -- a big catch at a
#      station already excluded by use.for.CPUE cannot affect the time series.
big <- cp_r %>%
  filter(species %in% c("Pacific cod", "walleye pollock"), CPUE > 1000) %>%
  inner_join(st_r %>% filter(use.for.CPUE == "yes") %>%
               select(Station, subsample), by = "Station") %>%
  select(Station, Site, Bay, species, CPUE, subsample) %>%
  arrange(desc(CPUE))
if (nrow(big) > 0) {
  cat("        catches over 1000, with the subsample flag:\n")
  cat(paste0("        ", capture.output(print(as.data.frame(big), row.names = FALSE)),
             collapse = "\n"), "\n")
}
check("every catch over 1000 came from a subsampled set",
      all(big$subsample == "yes"),
      big %>% filter(subsample != "yes"))

section(sprintf("%d ROW-COUNT RECONCILIATION", YEAR))

# retained stations, after the Chief Cove and May exclusions add.cpue.R applies
n_wgoa <- st_r %>%
  filter(use.for.CPUE == "yes", trimws(Bay) != "Chief Cove", julian > 151) %>%
  nrow()
n_kod  <- nrow(kd)
n_ser  <- sum(dat$year == YEAR)
check(sprintf("series rows for %d == wGOA(%d) + Kodiak(%d) = %d",
              YEAR, n_wgoa, n_kod, n_wgoa + n_kod),
      n_ser == n_wgoa + n_kod, sprintf("series has %d", n_ser))

}  # end of the YEAR-specific checks

section("ALL YEARS -- STRUCTURAL")

# 11. no negative age-0 counts (would mean age-1 was over-subtracted)
neg <- dat %>% filter(cod.age.0 < 0 | pollock.age.0 < 0)
check("no negative age-0 counts in any year", nrow(neg) == 0, neg)

# 12. no leading/trailing whitespace in labels -- this is the "Japanese " class
ws <- dat %>% filter(bay != trimws(bay) | site != trimws(site)) %>% distinct(year, bay, site)
check("no leading/trailing whitespace in bay or site", nrow(ws) == 0, ws)

# 13. bay names that collapse to the same string once case and spacing are
#     normalised -- catches "Cook"/"cook ", "Japanese"/"Japanese "
bay_collapse <- dat %>%
  distinct(bay) %>%
  mutate(key = tolower(trimws(bay))) %>%
  group_by(key) %>% filter(n() > 1) %>% ungroup()
check("no bay names differing only by case/spacing", nrow(bay_collapse) == 0, bay_collapse)

# 14. near-duplicate site names within a bay -- this is the "Kiul-4"/"Kilu-4"
#     and "Kil-1"/"Kilu-1" class of transposition and abbreviation typo.
#     Sites in the same bay that differ only by a number (Kilu-1 vs Kilu-2) or by
#     a compass direction (Eelgrass North vs Eelgrass South) are real, distinct
#     stations, so strip both before deciding whether a pair is suspicious.
strip_noise <- function(x) {
  y <- tolower(x)
  y <- gsub("[0-9]", "", y)
  y <- gsub("\\b(north|south|east|west|n|s|e|w|upper|lower|inner|outer|mid|middle)\\b",
            "", y)
  gsub("[^a-z]", "", y)
}

near_dupes <- dat %>%
  distinct(bay, site) %>%
  group_by(bay) %>%
  group_modify(~ {
    s <- .x$site
    if (length(s) < 2) return(tibble())
    d <- adist(s)
    idx <- which(d > 0 & d <= 2, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    if (nrow(idx) == 0) return(tibble())
    tibble(site_a = s[idx[, 1]], site_b = s[idx[, 2]], edit_dist = d[idx])
  }) %>% ungroup()

suspicious <- near_dupes %>%
  filter(strip_noise(site_a) != strip_noise(site_b))
check("no near-duplicate site names within a bay", nrow(suspicious) == 0,
      as.data.frame(suspicious))

# 15. low-replication sites. add.cpue.R drops every site sampled in fewer than
#     three years (Ugak-1, Sand-1, Kilu-9, Kai-5, BB-5, Kilu-7, Kilu-8, BB-1),
#     so the series should contain none. A failure here means either a site-name
#     typo splitting an established site, or a genuinely new station that needs a
#     keep/drop decision and an entry in drop.sites.
years_per_site <- dat %>% distinct(bay, site, year) %>% count(bay, site, name = "years")
sets_per_site <- dat %>% count(bay, site, name = "sets")
replication <- left_join(years_per_site, sets_per_site, by = c("bay", "site"))

sparse <- replication %>% filter(years <= 2 | sets <= 2) %>%
  arrange(years, sets, bay, site)
check("no site sampled in fewer than three years", nrow(sparse) == 0,
      as.data.frame(sparse))

cat(sprintf("        least-replicated site now spans %d years\n",
            min(replication$years)))

# 16. Two rows for one site on one day are usually legitimate replicate seine
#     sets, so they are only advisory. Rows that are identical in the catch
#     columns as well are far more likely to be a record entered twice.
same_day <- dat %>% count(year, bay, site, julian) %>% filter(n > 1)
identical_rows <- dat %>%
  count(year, bay, site, julian, cod.age.0, pollock.age.0) %>% filter(n > 1)

check("no repeated site x day rows (replicate sets are expected)",
      nrow(same_day) == 0, as.data.frame(same_day), warn_only = TRUE)

# An identical catch is only suspicious when nothing else separates the two
# rows. Site names ARE identical across replicate sets, so they never help. The
# 2006-2020 base file does carry a temperature column that the time series drops
# (add.cpue.R removes it), and two sets with different temperatures are two real
# hauls. Confirmed against the Time column in "Kodiak gadid CPUE 2006-2021.csv":
#   2009 Cook, julian 196   -- all 8 sites seined twice, 5-120 min apart. Real,
#                              and the second pass is now re-dated to julian 197.
#   2010 Cook, julian 239   -- Laminaria East at 1032 and 1215; the 1215 set was
#                              a mis-entered Middle Cove and is now relabelled.
#   2012 Anton Larsen, 195  -- identical 258/21 catch with no temperature on the
#                              second row; now dropped as an accidental repeat.
# This tests the SERIES, not the raw base file, so corrections applied in
# add.cpue.R count as resolved. Temperature is looked up from the base file only
# as evidence. A repeat with no temperature to consult stays flagged.
identical_temps <- read.csv("./data/cpue.data.csv") %>%
  group_by(year, site, julian, cod.age.0, pollock.age.0) %>%
  summarize(temps = paste(temperature, collapse = " | "),
            separable = n_distinct(temperature, na.rm = TRUE) > 1,
            .groups = "drop")

flagged <- identical_rows %>%
  left_join(identical_temps,
            by = c("year", "site", "julian", "cod.age.0", "pollock.age.0")) %>%
  mutate(separable = replace_na(separable, FALSE),
         temps = replace_na(temps, "no temperature recorded"))

unexplained <- flagged %>% filter(!separable) %>% select(-separable)
explained <- flagged %>% filter(separable)

if (nrow(explained) > 0) {
  cat(sprintf("        %d identical-catch pairs are separated by temperature (real replicate sets)\n",
              nrow(explained)))
}
check("identical-catch rows are all separable by temperature",
      nrow(unexplained) == 0, as.data.frame(unexplained))

# 17. julian plausible in every year
check("julian within 130-280 in all years", all(dat$julian >= 130 & dat$julian <= 280),
      dat %>% filter(julian < 130 | julian > 280))

# 18. the two standing exclusions actually held, in every year. These replace
#     the raw-file versions: the raw files may contain both, the series must not.
check("no May samples (julian <= 151) reached the series",
      all(dat$julian > 151), dat %>% filter(julian <= 151))
# The four low-replication bays must never reach the series. Two of them carry
# names other than the one they are known by: Canton Harbor is spelled "Caton
# Harbor" in the raw files, and the Sanak sites sit under bay "NE Harbor".
# Kujulik has never appeared in a cpue/site/length file at all.
excluded_bays <- c("Caton Harbor", "Canton Harbor", "NE Harbor", "Sanak",
                   "Chief Cove", "Kujulik")
check("no Chief Cove / Caton Harbor / NE Harbor / Kujulik rows reached the series",
      !any(dat$bay %in% excluded_bays), dat %>% filter(bay %in% excluded_bays))

section("ALL YEARS -- COVERAGE AND MAGNITUDE")

cat("\nRows and sites sampled per year:\n")
print(as.data.frame(dat %>% group_by(year) %>%
  summarize(rows = n(), bays = n_distinct(bay), sites = n_distinct(site), .groups = "drop")))

cat("\nMean age-0 CPUE per year:\n")
print(as.data.frame(dat %>% group_by(year) %>%
  summarize(cod = round(mean(cod.age.0), 1),
            pollock = round(mean(pollock.age.0), 1), .groups = "drop")))

# 18. flag the update year if its mean is outside the range of all prior years
prior <- dat %>% filter(year < YEAR)
cur   <- dat %>% filter(year == YEAR)
for (sp in c("cod.age.0", "pollock.age.0")) {
  pm <- prior %>% group_by(year) %>% summarize(m = mean(.data[[sp]]), .groups = "drop") %>% pull(m)
  cm <- mean(cur[[sp]])
  check(sprintf("%d mean %s = %.1f vs historical range %.1f - %.1f",
                YEAR, sp, cm, min(pm), max(pm)),
        cm >= min(pm) && cm <= max(pm), warn_only = TRUE)
}

# ------------------------------------------------------------------ summary --
section("SUMMARY")
res <- bind_rows(results)
print(table(res$status))
failed <- res %>% filter(status == "FAIL")
if (nrow(failed) > 0) {
  cat("\nFAILED CHECKS:\n"); print(as.data.frame(failed))
} else {
  cat("\nNo hard failures.\n")
}
warned <- res %>% filter(status == "WARN")
if (nrow(warned) > 0) {
  cat("\nWARNINGS (review, not necessarily errors):\n"); print(as.data.frame(warned))
}
