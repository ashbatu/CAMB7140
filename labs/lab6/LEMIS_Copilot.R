# Introduction ----
# The goal of this script is to analyze data from LEMIS, tracking the import of animals and animal products into US ports.
# We'll use a combination of tidyverse and ggplot2 to explore the data and answer a series of questions.

# Load tidyverse and plotly packages ----
library(tidyverse)
library(plotly)

# read in lemis data (lemis_cleaned.tsv) using the readr package ----
lemis <- read_tsv("lemis_cleaned.tsv") # this reads in the lemis data from the file "lemis_cleaned.tsv" and stores it in a variable called 'lemis'.  You can replace "lemis_cleaned.tsv" with the path to your own lemis data file.  Make sure to include the correct file extension (.tsv, .csv, etc.) depending on the format of your data file.

# let's experiment with the prompt above.  What happens if we change .tsv to .csv...or remove these details altogether?


# take a look at the unique levels for any of the variables in the lemis dataset
unique(lemis$description) # this shows the unique levels for the 'description' variable in the lemis dataset.  You can replace 'description' with any other variable name to see the unique levels for that variable.

# Question 1: Use the filter function from dplyr to find the most common (by 'quantity') live mammal taken from the wild for import into the US. ----
lemis.live <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::filter(class == "Mammalia") %>%
  dplyr::filter(source == "Specimens taken from the wild") %>%
  dplyr::arrange(desc(quantity))


# Question 2: Building on your analysis above, filter the dataset to get live animals (use 'generic_name') imported for the purposes of "scientific" or "Biomedical research"
lemis.live.science <- lemis.live %>%
  dplyr::filter(purpose %in% c("Scientific", "Biomedical research")) %>%
  dplyr::arrange(desc(quantity))

# plot this lemis.live data you generated above using ggplot2
lemis.live.science.plot <- ggplot(lemis.live.science, aes(x = generic_name, y = quantity)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Live Mammals Imported for Scientific or Biomedical Research",
       x = "Generic Name",
       y = "Quantity")

# generate the plot above again, but this time using ploty
lemis.live.science.plotly <- ggplotly(lemis.live.science.plot)


# Question 3: From which countries do we import macaques? ----
lemis.macaques <- lemis %>%  
  dplyr::filter(generic_name == "MACAQUE") %>% 
  distinct(country_origin)

# Question 4: From which countries from which we import live bats? ----
lemis.bats <- lemis %>%  
  dplyr::filter(generic_name == "BAT") %>% 
  distinct(country_origin)


# Question 5: For what 'purpose' are these imported bats used? ----
lemis.bats.purpose <- lemis %>%  
  dplyr::filter(generic_name == "BAT") %>% 
  distinct(purpose)

# Question 6: How does the type of bat (use 'specific_name') imported differ between countries (use 'facet_wrap' in your ggplot code)? ----
lemis.bats.specific <- lemis %>%  
  dplyr::filter(generic_name == "BAT") %>% 
  ggplot(aes(x = specific_name, y = quantity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ country_origin) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Types of Bats Imported by Country",
       x = "Specific Name",
       y = "Quantity")


# Question 7: Identify the most expensive (by ‘value’) shipment of "Live specimens (live animals or plants)" to enter the US? ----
lemis.live.expensive <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::arrange(desc(value)) %>%
  head(1)


# Question 8: Identify the most expensive shipment of any kind (live or not)? ----
lemis.expensive <- lemis %>%
  dplyr::arrange(desc(value)) %>%
  head(1)


