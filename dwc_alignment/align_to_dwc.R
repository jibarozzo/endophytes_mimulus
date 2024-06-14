# Script to transform data to Darwin Core for publication at GBIF


# load libaries -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(Biostrings)

# read in data ------------------------------------------------------------

ASV_tab <- read_csv(file = "clean_data/ASV_tables/ASVs_8450_cleaned_10_percent.csv")
samples_tab <- read_xlsx(path = "field_data/Mimulus_CH2_Field_Survey.xlsx", sheet = "Ch2_FieldSurvey")
seq_tab <- readDNAStringSet(filepath = "clean_data/ASV_tables/ASVs_8450_amplicons.fa")
dna_pcr_tab <- read_xlsx(path = "field_data/Mimulus_CH2_Field_Survey.xlsx", sheet = "DNA_PCR_progress")
tax_tab <- read_csv(file = "clean_data/taxonomy/TAXA_8450_cleaned_10_percent.csv")

# reshape data ------------------------------------------------------------

ASV_tab <- ASV_tab %>% 
    dplyr::rename(ASV_name = ...1) %>% 
    pivot_longer(names_to = "sample", values_to = "organismQuantity", cols = -ASV_name) %>% 
    group_by(sample) %>% 
    mutate("sampleSizeValue" = sum(organismQuantity))

samples_tab_with_traits <- samples_tab %>% 
    dplyr::rename(Plant_Sample = Sample) %>% 
    select(Site:Plant_Sample, Unique_ID...15:Notes2) %>% 
    distinct()

samples_tab_no_traits <- samples_tab %>% 
    dplyr::rename(Plant_Sample = Sample) %>% 
    select(Site:Plant_Sample, Unique_ID...15:Sampling_time...17, Habitat) %>% 
    filter(!is.na(Unique_ID...15)) %>% 
    distinct()

seq_tab <-  seq_tab %>%
    as.character(use.names=TRUE) %>%
    tibble::enframe(name = "ASV_name", value = "DNA_sequence")

DNA_tab <- dna_pcr_tab %>% 
    mutate(.keep = "none",
           Unique_ID,
           concentration = `DNA_conc._ng/µl`,
           nucl_acid_ext = Extraction_Method) %>% 
    filter(!is.na(concentration)) %>% 
    distinct()

tax_tab <- tax_tab %>%
    dplyr::rename(ASV_name = ...1,
                  Fungal_Species = Species)
# create DwC tables -------------------------------------------------------

# Join tables
df <- left_join(x = ASV_tab, y = tax_tab, by = "ASV_name") %>% 
    left_join (y = samples_tab_no_traits, by = c("sample" = "Unique_ID...15")) %>% 
    left_join( y = seq_tab, by = "ASV_name") %>% 
    left_join( y = DNA_tab, by = c("sample" = "Unique_ID"))

glimpse(df)

# occurrence table --------------------------------------------------------

occ_tab <- df %>%
    mutate(
        .keep = "none",
        #analogous to transmutate()
        occurrenceID = paste(sample, ASV_name, sep = ":"),
        organismQuantity = organismQuantity,
        organismQuantityType = "DNA sequence reads",
        occurrenceStatus = case_when(organismQuantity > 0 ~ "present", organismQuantity == 0 ~ "absent"),
        eventDate = paste(Sample_Date...16, Sampling_time...17),
        eventID = NA,
        sampleSizeValue = sampleSizeValue,
        sampleSizeUnit = "DNA sequence reads",
        samplingProtocol = NA,
        eventType = "Sample",
        decimalLatitude = Latitude,
        decimalLongitude = Longitude,
        coordinatePrecision = NA,
        coordinateUncertaintyInMeters = NA,
        locationID = NA,
        locality = "Yosemite National Park, CA, USA",
        countryCode = "US",
        continent = "North America",
        minimumElevationInMeters = Elevation_m,
        maximumElevantionInMeters = Elevation_m,
        basisOfRecord = "MaterialSample",
        recordedBy = "Bolívar Aponte Rolón",
        materialSampleID = NA,
        associatedSequences = NA,
        identificationRemarks = NA,
        identificationReferences = NA,
        verbatimIdentification = paste0(Kingdom, Phylum, Class, Order, Family, Genus, Species.x),
        scientificName = NA,
        scientificNameID = NA,
        taxonID = NA,
        taxonRank = NA
    ) %>%
    glimpse()

                      
           Plant_occ_core <- df %>%
               mutate(.keep = "none", #analogous to transmutate()
                      occurrenceID = sample, #This is not sufficiently unique, placeholder for now
                      individualCount = 1,
                      occurrenceStatus = "present",
                      eventDate = paste(Sample_Date...16, Sampling_time...17),
                      eventID = NA,
                      sampleSizeValue = NA,
                      sampleSizeUnit = NA,
                      samplingProtocol = NA,
                      eventType = "Sample",
                      decimalLatitude = Latitude,
                      decimalLongitude = Longitude,
                      coordinatePrecision = NA,
                      coordinateUncertaintyInMeters = NA,
                      locationID = NA,
                      locality = "Yosemite National Park, CA, USA",
                      countryCode = "US",
                      continent = "North America",
                      minimumElevationInMeters = Elevation_m,
                      maximumElevantionInMeters = Elevation_m,
                      basisOfRecord = "MaterialSample",
                      recordedBy = "Bolívar Aponte Rolón",
                      materialSampleID = NA,
                      associatedSequences = NA,
                      identificationRemarks = NA,
                      identificationReferences = NA,
                      verbatimIdentification = Plant_Species,
                      scientificName = NA,
                      scientificNameID = NA,
                      taxonID = NA,
                      taxonRank = NA
               ) %>%
               glimpse()

           #all occurrences
           all_occ <- left_join(Plant_occ_core, Fungal_occ_core)
           
# extended MeasurementOrFact extension ------------------------------------


           eMOF <- samples %>%
               pivot_longer(names_to = measurementType,
                            values_to = measurementValue,
                            cols = Leaf_thickness,ACI,Leaf_toughness, Leaf_lobe_index) %>%
               mutate(.keep = "none",
                      occurrenceID = Unique_ID...15,
                      measurementType,
                      measurementValue
               ) %>%
               glimpse()
           #extended Measurement Or Fact extension
           eMOF <- samples %>%
               pivot_longer(cols = c(Leaf_thickness,
                                     ACI,
                                     Leaf_toughness,
                                     Leaf_lobe_index),
                            names_to = 'measurementType',
                            values_to = 'measurementValue',
               ) %>%
               mutate(.keep = "none",
                      occurrenceID = Unique_ID...15,
                      measurementType,
                      measurementTypeID = NA,
                      measurementValue,
                      measurementValueID = NA,
                      measurmentUnit = NA,
                      measurementUnitID = NA
               ) %>%
               glimpse()

# DNA derived Data extension ----------------------------------------------


           DNA <- df %>%
               mutate(.keep = "none",
                      occurrenceID = paste(sample, ASV_name, sep = ":"),
                      DNA_sequence,
                      env_broad_scale = NA,
                      env_medium = NA,
                      pcr_primer_forward = NA,
                      pcr_primer_reverse = NA,
                      pcr_primer_name_forward = NA,
                      pcr_primer_name_reverse = NA,
                      pcr_primer_reference = NA,
                      target_gene = NA,
                      target_subfragment = NA)
           