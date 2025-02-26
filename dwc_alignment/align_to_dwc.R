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

# reshape data and rename a few columns for convenience -------------------

ASV_tab <- ASV_tab %>% 
    dplyr::rename(ASV_name = ...1) %>% 
    pivot_longer(names_to = "Plant_Sample", values_to = "organismQuantity", cols = -ASV_name) %>%
    group_by(Plant_Sample) %>% 
    mutate("sampleSizeValue" = sum(organismQuantity))

samples_tab_with_traits <- samples_tab %>% 
    dplyr::rename(Plant_Sample = Sample,
                  Leaf_Sample = Unique_ID...15) %>% 
    select(Site:Plant_Sample, Leaf_Sample:Notes2) %>% 
    distinct()

samples_tab_no_traits <- samples_tab %>% 
    dplyr::rename(Plant_Sample = Sample,
                  Leaf_Sample = Unique_ID...15) %>% 
    select(Site:Plant_Sample, Leaf_Sample:Sampling_time...17, Habitat) %>% 
    filter(!is.na(Leaf_Sample)) %>% 
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
    left_join (y = samples_tab_no_traits, by = "Leaf_Sample") %>% 
    left_join( y = seq_tab, by = "ASV_name") %>% 
    left_join( y = DNA_tab, by = c("Leaf_Sample" = "Unique_ID"))

glimpse(df)

# occurrence table --------------------------------------------------------

# Find definitions here: https://dwc.tdwg.org/terms/

occ_tab_ASV <- df %>%
    mutate(
        .keep = "none", #analogous to transmutate()
        occurrenceID = paste(Leaf_Sample, ASV_name, sep = ":"),
        organismQuantity = organismQuantity,
        organismQuantityType = "DNA sequence reads",
        occurrenceStatus = case_when(organismQuantity > 0 ~ "present", organismQuantity == 0 ~ "absent"),
        verbatimEventDate = paste("Sample_Date: ", Sample_Date...16, "Sample_time: ", Sampling_time...17),
        
        # Jump through some hoops to format the date-time as ISO8601 string in UTC
        eventDate = case_when(is.na(Sampling_time...17) & !is.na(Sample_Date...16) ~ Sample_Date...16, 
                              !is.na(Sampling_time...17) & !is.na(Sample_Date...16) ~ paste(Sample_Date...16, Sampling_time...17)) %>% 
            lubridate::parse_date_time(orders = c("mdy", "mdy %I:%M:%S %p"), tz = "US/Pacific") %>% 
            lubridate::with_tz(time, tzone = "UTC") %>% 
            as.character() %>% 
            stringr::str_replace(pattern = "$", replacement = "Z") %>% 
            stringr::str_replace(pattern = "\\s", replacement = "T") %>% 
            stringr::str_replace(pattern = "T07:00:00Z", replacement = ""),
        
        eventID = NA, # Maybe something like "MIM2:LEF:ASV:202120220407" for example. This is MIM2 project, LEF = LEAF, ASV = ASV, 202120220407 = unique identifier (2021 and 2022 between months 04 and 07)
        sampleSizeValue = sampleSizeValue,
        sampleSizeUnit = "DNA sequence reads",
        samplingProtocol = "Aponte Rolón, B. (2023). High‐Molecular‐Weight SPRI‐aided DNA extraction from Mimulus (Phrymaceae) leaf tissue.
        https://dx.doi.org/10.17504/protocols.io.bp2l6xn8rlqe/v2",
        eventType = "Sample",
        decimalLatitude = Latitude,
        decimalLongitude = Longitude,
        coordinatePrecision = "0.000278 (nearest second)",
        coordinateUncertaintyInMeters = NA,
        geodeticDatum = 'WGS84', #double check this
        locationID = Site,
        locality = NA,  # Must be derived from lat/lon, could just say Sierra Nevada, CA 
        # ^ Maybe this can be a polygon?
        
        # 3 Locations
        # Yosemite National Park, CA, USA
        # Stanislaus National Forest, CA, USA
        # Sierra National Forest, CA, USA
        
        countryCode = "US",
        continent = "North America",
        minimumElevationInMeters = Elevation_m,
        maximumElevantionInMeters = Elevation_m,
        basisOfRecord = "MaterialSample",
        recordedBy = "Bolívar Aponte Rolón",
        materialSampleID = NA, # This is gonna be either Plant_Sample or Unique_ID (from the initial data entry unique identifier)
        associatedSequences = NA, # BioProject or SRA accession number
        identificationRemarks = NA,
        identificationReferences = NA, #Jepsn guide
        verbatimIdentification = paste0(Kingdom, Phylum, Class, Order, Family, Genus, Fungal_Species),
        
        #Will need to be aligned with GBIF taxonomic backbone
        scientificName = NA, 
        scientificNameID = NA,
        taxonID = NA,
        taxonRank = NA
    )

                      
occ_tab_plant <- df %>%
    mutate(
        .keep = "none",
        occurrenceID = paste(Site, Plant_Sample, sep = ":"),#This is not sufficiently unique, placeholder for now
        individualCount = 1, #not sure if I'm interpreting this correctly
        occurrenceStatus = "present",
        
        # Jump through some hoops to format the date-time as ISO8601 string in UTC
        eventDate = case_when(is.na(Sampling_time...17) & !is.na(Sample_Date...16) ~ Sample_Date...16, 
                              !is.na(Sampling_time...17) & !is.na(Sample_Date...16) ~ paste(Sample_Date...16, Sampling_time...17)) %>% 
            lubridate::parse_date_time(orders = c("mdy", "mdy %I:%M:%S %p"), tz = "US/Pacific") %>% 
            lubridate::with_tz(time, tzone = "UTC") %>% 
            as.character() %>% 
            stringr::str_replace(pattern = "$", replacement = "Z") %>% 
            stringr::str_replace(pattern = "\\s", replacement = "T") %>% 
            stringr::str_replace(pattern = "T07:00:00Z", replacement = ""),
        
        eventID = # Maybe something like "MIM2:LEF:ASV:202120220407" for example. This is MIM2 project, LEF = LEAF, ASV = ASV, 202120220407 = unique identifier (2021 and 2022 between months 04 and 07),
        sampleSizeValue = 400,
        sampleSizeUnit = milligram,
        samplingProtocol = NA,
        eventType = "Sample",
        decimalLatitude = Latitude,
        decimalLongitude = Longitude,
        coordinatePrecision = NA,
        coordinateUncertaintyInMeters = NA,
        locationID = Site,
        locality = NA,  # Must be derived from lat/lon, could just say Sierra Nevada, CA 
        # ^ Maybe this can be a polygon?
        
        # 3 Locations
        # Yosemite National Park, CA, USA
        # Stanislaus National Forest, CA, USA
        # Sierra National Forest, CA, USA
        
        countryCode = "US",
        continent = "North America",
        minimumElevationInMeters = Elevation_m,
        maximumElevantionInMeters = Elevation_m,
        basisOfRecord = "MaterialSample",
        recordedBy = "Bolívar Aponte Rolón",
        materialSampleID = NA, # This is gonna be either Plant_Sample or Unique_ID (from the initial data entry unique identifier)
        associatedSequences = NA,
        identificationRemarks = NA,
        identificationReferences = NA,
        verbatimIdentification = Species,
        
        #Will need to be aligned with GBIF taxonomic backbone
        scientificName = NA, 
        scientificNameID = NA,
        taxonID = NA,
        taxonRank = NA
    )

#all occurrences
all_occ <- rbind(occ_tab_plant, occ_tab_ASV)

# extended MeasurementOrFact extension ------------------------------------

# Find definitions here: https://rs.gbif.org/extension/obis/extended_measurement_or_fact.xml

           eMOF <- samples_tab_with_traits %>%
               pivot_longer(names_to = "measurementType",
                            values_to = "measurementValue",
                            cols = c("Leaf_thickness","ACI","Leaf_toughness", "Leaf_lobe_index")) %>% #Look for LMA

                   mutate(.keep = "none",
                      occurrenceID = Leaf_Sample,
                      measurementType,
                      measurementValue,
                      
                      # Get identifiers for these 
                      measurementTypeID = NA,
                      measurementValueID = NA,
                      measurmentUnit = NA,
                      measurementUnitID = NA
                      
               )
           
# DNA derived Data extension ----------------------------------------------

# Find definitions here: https://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2022-02-23.xml

           DNA <- df %>%
               mutate(.keep = "none",
                      occurrenceID = paste(Leaf_Sample, ASV_name, sep = ":"),
                      DNA_sequence,
                      env_broad_scale = 'alpine biome[ENVO:01001835]',
                      env_medium = 'meadow soil[ENVO:00005761] | granite ENVO_01000356', #make this ASV specific
                      env_local_scale = 'plant-associated environment[ENVO:01001001]', #check if there is something more specific
                      pcr_primer_forward = 'CACTCTTTCCCTACACGACGCTCTTCCGATCTCTTGGTCATTTAGAGGAAGTAA',
                      pcr_primer_reverse = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCTGCGTTCTTCATCGATGC',
                      pcr_primer_name_forward = 'ITS1F',
                      pcr_primer_name_reverse = 'ITS2',
                      pcr_primer_reference = NA, #insert DOIs for references, pipe-delimited e.g. DOI#1|DOI#2
                      target_gene = 'ITS',
                      target_subfragment = 'ITS1',
                      lib_layout = 'paired',
                      nucl_acid_ext = NA, # extraction protocol,
                      nucl_acid_amp = NA, # amplification protocol 
                      samp_size = NA, #mass of sample that was collected
                      seq_meth = 'Illumina MiSeq...', #DOuble check this
                      mid = NA,
                      adapters = NA,
                      neg_cont_type = NA, #negative control
                      pos_cont_type = NA # positive control
                      )

# Information to fill NA values -------------------------------------------

           