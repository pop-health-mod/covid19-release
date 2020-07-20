# Plotting cumulative deaths and international comp
# Covid-19
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lemon)


# read in Canadian raw data from John_hopkins data
cad_raw <- read.csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
cad_prov <- cad_raw %>%
  filter(Country.Region == "Canada") %>%
  filter(!(Province.State %in% c("Diamond Princess", "Grand Princess"))) %>%
  select(Province.State) %>%
  unique()

# Population data
pop_raw <- read.csv("Data/WPP2019_TotalPopulationBySex.csv", stringsAsFactors = FALSE)
pop_sub <- subset(pop_raw, MidPeriod == "2019.5")
pop <- data.frame(country = pop_sub$Location,
                  province = pop_sub$Location,
                  pop_total = pop_sub$PopTotal * 1000,
                  stringsAsFactors = FALSE) %>%
  mutate(country = ifelse(country == "United States of America", "US", country)) %>%
  mutate(province = ifelse(province == "United States of America", "US", province)) %>%
  filter(country %in% c("France", "United Kingdom", "US"))

# Population Canada
pop_can <- read.csv("Data/1710000901_databaseLoadingData.csv") %>%
  select(GEO, VALUE) %>%
  filter(GEO != "Canada") %>%
  mutate(country = "Canada") %>%
  rename(province = GEO, pop_total = VALUE)


# Merge population datasets
pop_tot <- bind_rows(pop, pop_can)

dat <- cad_raw %>%
  filter(Country.Region %in% c("Canada", "France", "US", "United Kingdom")) %>%
  select(-Lat, -Long) %>%
  reshape2::melt() %>%
  mutate(variable = gsub("X", "", variable)) %>%
  mutate(variable = as.Date(variable, format = "%m.%d.%y")) %>%
  filter(Province.State %in% c(as.character(cad_prov$Province.State), "")) %>%
  rename(province = Province.State, country = Country.Region, date = variable, 
         deaths = value) %>%
  mutate(coloring = case_when(
    province %in% c("New Brunswick", "Newfoundland and Labrador",
                    "Nova Scotia", "Prince Edward Island") ~ "Atlantic",
    province %in% c("Manitoba", "Saskatchewan", "Alberta") ~ "Prairies",
    province %in% c("Northwest Territories", "Yukon", "Nunavut") ~ "Territories",
    TRUE ~ as.character(province)
  )) %>%
  mutate(province = ifelse(province == "",
                           as.character(country),
                           as.character(province))) %>%
  mutate(labels = ifelse(province %in% c("France", "United Kingdom", "US",
                                         "Quebec", "Ontario", "British Columbia",
                                         "Alberta", "Nova Scotia"),
                             as.character(province),
                             "")) %>%
  mutate(labelscnt = ifelse(province %in% c("France",
                                            "US",
                                            "United Kingdom"),
                            as.character(province),
                            ""))


datfull <- merge(dat, pop_tot, by=c("country", "province")) %>%
  mutate(death_mill = 1e6 * deaths / pop_total)

# We plot the processed data
plot1 <- ggplot() +
  theme_minimal() +
  geom_point(data = datfull %>% filter(date == max(as.Date(date)) & !(province %in% c("France", "US", "United Kingdom"))),
             aes(x = date, y = death_mill, color = coloring, group = province), show.legend = FALSE) +
  geom_point(data = datfull %>% filter(date == max(as.Date(date)) & labelscnt != ""),
             aes(x = date, y = death_mill), color = "darkgrey", show.legend = FALSE) +
  geom_text(data = datfull %>% filter(date == max(as.Date(date))),
             aes(x = date, y = death_mill, label = labels),
            fontface = "bold", size = 2, hjust=0, nudge_x = 2, vjust=0.5) +
  geom_line(data = datfull %>% filter(country == "Canada"),
            aes(x = date, y = death_mill, color = coloring, group = province),
            show.legend = FALSE, size=0.75) +
  geom_line(data = datfull %>% filter(country != "Canada"),
            aes(x = date, y = death_mill, group = province), col = "darkgrey", size=0.75,
            show.legend = FALSE, lty = 2) +
  coord_cartesian(ylim = c(0, max(datfull$death_mill) * 1.01)) +
  scale_x_date(limits=c(as.Date("2020-03-01"), Sys.Date() + 30), expand=c(0,1),
               breaks = "2 weeks", date_labels = "%b %d") +
  scale_color_manual(values = nationalparkcolors::park_palette("ArcticGates")) +
  xlab("Date") + ylab("Cumulative cases per million") +
  theme(text = element_text("Source Sans Pro"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9))
  
plot1

ggsave("cumcases.png", plot = plot1, width = 4, height = 5,
       device = "png",dpi = 320)
