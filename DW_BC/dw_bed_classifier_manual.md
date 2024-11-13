# Deep-Water Bed Classifier

2024-11-13

**Authors:** Budai, S., Colombera, L., McArthur, A. and Patacci, M.

contact: soma.budai@unipv.it

------------------------------------------------------------------------

This document explain how the `dw_bed_classifier.R` works and how it
uses deep-water sedimentary facies data to discriminate and classify
lithological (sand/gravel) bed types of deep-water successions produced
by sediment gravity flows.

#### 1. Loading required packages

``` r
library(dplyr)
library(stringr)
```

In order to run, the script requires the `dplyr` and `stringr` package
for the manipulation dataframes and strings.

#### 2. Loading dummy dataset

The script utilizes one mandatary dataset containing facies data and two
other optional datasets containing data on deep-water depositional
systems and architectural elements linked to the deposits. In the
illustrative example presented here, the script is supplied with three
dummy datasets containing partially redacted data extracted from the
DMAKS database. If the script is opened separately not as part of the R
project provided path to the files should be changed. Errors can be
avoided by pasting the whole path to the files!

``` r
facies_data_dummy <- read.csv('DW_BC/dummy_dataset/facies_data_dummy.csv', header = TRUE)
system_data_dummy <- read.csv('DW_BC/dummy_dataset/system_data_dummy.csv', header = TRUE)
element_data_dummy <-  read.csv('DW_BC/dummy_dataset/element_data_dummy.csv', header = TRUE)
```

#### 3. Function

The script contains a function used to assign facies data to beds
classified according to bed types.

#### 3.1 Function parameters

The function output is governed by a series of attributes including the
names of the input datasets: `facies_data`, `system_data` and
`element_data`. A series of Boolean parameters are used to specify
whether the bed classification should include vertical facies changes
and/or lamination data (see paper) and whether the depositional-system
and architectural-element information should be included in the output
data or not (`include_system_info`, `include_elemenet_info`).

``` r
bed_classifier <- function (facies_data = NULL, 
  include_vchanges = FALSE, 
  include_lam = FALSE, 
  system_data = NULL, 
  include_system_info = FALSE, 
  element_data = NULL, 
  include_element_info = FALSE) {...}
```

#### 3.2 Identifying deep-water beds

The initial part of the script carries out the filtering of the facies
data. It then identifies distinct beds based on the presence of bed
bouding surfaces (indicated by the ‘trans_from_below_event_b’ field) in
the coded sedimentary log and assigns an ID to each bed.

``` r
#identify beds
facies_data_bedboundary <- facies_data %>% 
  filter(!is.na(trans_from_below_event_b)) %>% 
  filter(facies_type %in% c('_G','_S','gS','sG','(g)S', 'mS', 'sM')) %>%
  filter(thickness_type %in% c('true (maximum)','true (not maximum)','apparent')) 
  
facies_data_bednr <- facies_data_bedboundary  %>% 
  mutate(new_bed = case_when(data_order!= lag(data_order)+1 | 
    trans_from_below_event_b == "Y" ~ 1, TRUE ~ 0)) %>% 
  mutate(artificial_bed_id = cumsum(new_bed))
  
has_y <- facies_data_bednr %>% 
  filter(trans_from_below_event_b == "Y") %>% 
  group_by(artificial_bed_id)  %>% 
  count(trans_from_below_event_b) %>% 
  dplyr::select(artificial_bed_id,n)
  
beds_raw <- left_join(facies_data_bednr, has_y, by = "artificial_bed_id") %>% 
  filter(!is.na(n))
```

#### 3.3 Calculate bed thickness

For each bed, its thickness and start and end ‘facies_ID’ is determined.

``` r
#calculate bed thickness
  
beds_raw <- beds_raw %>% 
  group_by(artificial_bed_id) %>% 
  mutate(bed_thickness = sum(thickness)) %>% 
  mutate(start_facies = min(facies_ID)) %>%
  mutate(end_facies = max(facies_ID)) %>% 
  ungroup()
```

#### 3.4 Calculate grain size fraction

For each deep-water bed, the thickness of gravel, sand and sandy-muddy
facies is calculated and additionally expressed as a fraction relative
to the bed thickness. If the fraction of sandy-muddy facies is higher
than zero, the bed will fall into the ‘SM’ grain-size category; if the
fraction of gravel is 1 then the bed is classified as ‘G’, if the
faction is 0 than the bed is classified as ‘S’. Fractions between 1 and
0 result in the assignment of ‘sG’ and ‘gS’ categories.

``` r
#grain size
  
beds_gs <- beds_raw %>%
  group_by(artificial_bed_id) %>%
  mutate(gravel_thickness = sum(ifelse(facies_type %in% c('sG','_G'), thickness, 0 ))) %>%
  mutate(gsand_thickness = sum(ifelse(facies_type %in% c('gS','(g)S'), thickness, 0 ))) %>%
  mutate(sand_thickness = sum(ifelse(facies_type %in% c('_S'), thickness, 0 ))) %>%
  mutate(heterolithic_thickness = sum(ifelse(facies_type %in% c('mS', 'sM'), 
      thickness, 0 ))) %>%
  mutate(gravel_fraction = round((gravel_thickness/bed_thickness),2)) %>%
  mutate(gsand_fraction = round((gsand_thickness/bed_thickness),2)) %>%
  mutate(sand_fraction = round((sand_thickness/bed_thickness),2)) %>%
  mutate(heterolithic_fraction = round((heterolithic_thickness/bed_thickness),2)) %>%
  ungroup() %>%
  mutate(gs_category = case_when(heterolithic_fraction > 0 ~ 'SM',
                                   gravel_fraction == 1 ~ 'G', 
                                   gravel_fraction == 0 ~ 'S', 
                                   gravel_fraction < 1 & gravel_fraction >= 0.5 ~ 'sG', 
                                   gravel_fraction < 0.5 & gravel_fraction > 0 ~ 'gS')) %>%
    distinct(artificial_bed_id, .keep_all = TRUE)
```

#### 3.5 Determine sharp vertical facies transitions within beds

The next part of the function only runs if the vertical facies changes
were selected to be included in the bed-type classification and in the
output dataframe. Numerical values are assigned to each facies in the
beds based on a combination of their facies and sand grain-sizes in a
fining order. These values range from 17 for gravel to 1 for very fine
sandy mud. Datasets for which the sand grain-size was not provided are
filtered out. For each bed, the script calculates the difference in
these assigned values between each facies and the overlying one. Then
the minimum and maximum of these differences for each bed, termed here
as ‘grain-size trend scores’, are determined. If the minimum and maximum
values are equal to zero or if the bed is composed of only one facies,
then the bed is classed as showing no vertical trend (category N). If
the maximum value of the calculated differences is negative, then a
fining-upward trend (F) is assigned to the bed; if the minimum value is
positive then the bed is classified as coarsening upward (C). If both
coarsening and fining trends are observed among the facies units
composing the bed at hand (i.e., if the minimum value is negative and
the maximum value is positive) without any preference in order, then the
bed falls into the ’B’ category, which flags the presence of both
grain-size trends.

``` r
if (include_vchanges == TRUE) {
    
    
  combined_gs_values <- data.frame(
        combined_gs = c('G',  'S',    'Svery coarse', 'Scoarse',  'Smedium',
                        'Sfine','Svery fine', 'MS','MSvery coarse','MScoarse',
                        'MSmedium','MSfine','MSvery fine','SM','SMvery coarse', 
                        'SMcoarse', 'SMmedium','SMfine','SMvery fine'),
        gs_value = c(17, 14,  16, 15, 14, 13, 12, 8,  10, 9,  8,  7,  6,  3,  
                      5,  4,  3,  2, 1)
  )
    
  no_sand_gs <- beds_raw %>% 
    dplyr::select(artificial_bed_id, facies_type, grain_size_sand) %>%
    mutate(simple_gs = case_when(facies_type %in% c('_G', 'sG')~ 'G', 
                                 facies_type %in% c('(g)S','_S','gS') ~ 'S', 
                                 facies_type %in% c('mS') ~ 'MS', 
                                 facies_type %in% c('sM') ~ 'SM')) %>%
    mutate(has_sand_gs = case_when(simple_gs %in% c('S','MS') & is.na(grain_size_sand)~'N',
                                   TRUE ~ 'Y')) %>%
    filter(has_sand_gs == 'N') %>% 
    distinct(artificial_bed_id) %>% 
    pull(artificial_bed_id) %>% as.list()
    
  combined_vchanges <- beds_raw %>% 
      mutate(simple_gs = case_when(facies_type %in% c('_G', 'sG')~ 'G', 
                                   facies_type %in% c('(g)S','_S','gS') ~ 'S', 
                                   facies_type %in% c('mS') ~ 'MS', 
                                   facies_type %in% c('sM') ~ 'SM')) %>%
      mutate(general_structure = case_when(is.na(general_structure) ~ 'nd', 
                                           TRUE ~ general_structure )) %>% 
      mutate(lamination_type = case_when(is.na(lamination_type) ~ 'nd', 
                                         TRUE ~ lamination_type )) %>%
      mutate(combined_gs = case_when(simple_gs == 'G' ~ simple_gs,
                                     is.na(grain_size_sand) ~ simple_gs, 
                                     !is.na(grain_size_sand) ~ paste(simple_gs,
                                            grain_size_sand, sep = ''))) %>%
      dplyr::select(facies_ID,log_ID,data_order,artificial_bed_id,trans_from_below,
                    trans_from_below_event_b,
                    general_structure, lamination_type, combined_gs) %>%
      left_join(combined_gs_values, by = 'combined_gs') %>%
      mutate(v_diff = case_when(lead(artificial_bed_id) != artificial_bed_id & lag(artificial_bed_id) != artificial_bed_id ~ 100,
                                is.na(lead(artificial_bed_id)) & lag(artificial_bed_id) != artificial_bed_id ~ 100,
                                is.na(lead(artificial_bed_id)) & lag(artificial_bed_id) == artificial_bed_id ~ 200,
                                lead(artificial_bed_id) != artificial_bed_id ~ 200, 
                                lead(trans_from_below == 'gradational') & (lead(general_structure) == general_structure | lead(lamination_type) == lamination_type)  ~ 0, 
                                TRUE ~ lead(gs_value)-gs_value )) %>%
      filter(v_diff != 200) %>%
      group_by(artificial_bed_id) %>%
      mutate(v_min = min(v_diff)) %>%
      mutate(v_max = max(v_diff)) %>%
      ungroup() %>%
      mutate(combined_trend = case_when(is.na(v_max) ~ 'N', v_max == 100 ~ 'N', v_min == 0 & v_max == 0 ~ 'N' ,v_max <= 0 & v_min < 0 ~ 'F', 
                                        v_min < 0 & v_max > 0 ~ 'B', v_min >= 0 & v_max > 0 ~ 'C')) %>%
      dplyr::select(artificial_bed_id,combined_trend) %>%
      distinct(artificial_bed_id, .keep_all = TRUE)
    
}
```

#### 3.5 Calculate lamination thickness fraction

The next step of the function only runs if information on the lamination
is selected to be included in the bed-type classification and in the
output dataframe. The script calculates the cumulative thickness of the
laminated or stratified portion (e.g. planar-parallel lamination, non
planar-parallel, ripple cross-lamination, mesoform-scale
cross-stratification and wavy cross-stratification) of each bed and
compares it to the bed thickness. It also filters out beds recorded in
logs where information on lamination of the deposits was missing. Based
on these thickness ratios, three categories can be assigned: beds that
lack any kind of lamination (m), bed that are fully made of laminated
facies (l) and beds that are only in part made of laminated facies (x).

``` r
if (include_lam == TRUE) {
    
    lamination_thickness <- beds_raw %>% 
      filter(general_structure == 'laminated')%>% 
      group_by(artificial_bed_id) %>% 
      mutate(laminated_thickness = sum(thickness)) %>% 
      ungroup() %>%
      dplyr::select(artificial_bed_id,laminated_thickness) %>% 
      distinct(artificial_bed_id, .keep_all = TRUE)
    
    log_lamination_list <-  facies_data %>% filter(!is.na(log_facies_suit_proportions_type)) %>% 
      filter(str_detect(log_facies_suit_proportions_type,'general structures')) %>% distinct(log_ID) %>% pull(log_ID) %>% as.list()
    
    # subset_lamination_list <-  system_data_dmaks %>% filter(!is.na(facies_suit_proportions_type)) %>% 
    #   filter(str_detect(facies_suit_proportions_type,'general structures')) %>% pull(subset_ID) %>% as.list()
    
    lamination_ratio <- left_join(beds_gs,lamination_thickness, by = 'artificial_bed_id') %>% 
      mutate(laminated_thickness = case_when(is.na(laminated_thickness)~ 0, TRUE ~ laminated_thickness)) %>%
      mutate(laminated_fraction = round((laminated_thickness/bed_thickness),2)) %>% 
      mutate(laminated_fraction = case_when(!log_ID %in% log_lamination_list ~ 100, TRUE ~ laminated_fraction )) %>%
      mutate(lam_category = case_when(laminated_fraction == 100 ~ 'nd',
                                      laminated_fraction == 1 ~ 'l', 
                                      laminated_fraction == 0 ~ 'm', 
                                      TRUE ~ 'x')) %>%
      dplyr::select(artificial_bed_id, laminated_fraction, lam_category)
    
    
}
```

#### 3.6 Create bed data table

Based on the values of Boolean attributes `include_vchanges` and
`include_lam`, the script creates a joint dataframe containing
information on the stored beds.

``` r
bed_table <- beds_gs %>% mutate(code = gs_category)
  
  
if (include_vchanges == TRUE & include_lam == FALSE) {
    
    bed_table <- left_join(bed_table,combined_vchanges, by = 'artificial_bed_id') %>%
      filter(!artificial_bed_id %in% no_sand_gs) %>%
      mutate(code = paste(code, combined_trend, sep = '-')) %>%
      mutate(code2 = combined_trend) %>%
      dplyr::select(artificial_bed_id, code, code2,
                    bed_thickness, gravel_fraction, heterolithic_fraction,
                    log_ID, element_ID, subset_ID, case_study_ID, system_ID)
    
} 
  
if (include_lam == TRUE & include_vchanges == FALSE) {
    
    bed_table <- left_join(bed_table,lamination_ratio, by = 'artificial_bed_id')%>%
      filter(lam_category != 'nd') %>%
      mutate(code = paste(code, lam_category, sep = '-')) %>%
      mutate(code2 = lam_category) %>%
      dplyr::select(artificial_bed_id, code, code2,
                    bed_thickness, gravel_fraction, heterolithic_fraction,
                    log_ID, element_ID, subset_ID, case_study_ID, system_ID)
}  
  
if (include_vchanges == TRUE & include_lam == TRUE) {
    
    bed_table <- left_join(bed_table,combined_vchanges, by = 'artificial_bed_id') %>% 
      left_join(lamination_ratio, by = 'artificial_bed_id')%>%
      filter(!artificial_bed_id %in% no_sand_gs) %>%
      filter(lam_category != 'nd') %>%
      mutate(code = paste(code, combined_trend, lam_category, sep = '-')) %>%
      mutate(code2 = paste(combined_trend, lam_category, sep = '-')) %>%
      dplyr::select(artificial_bed_id, code, code2,
                    bed_thickness, gravel_fraction, heterolithic_fraction,
                    gs_category, combined_trend, lam_category,
                    log_ID, element_ID, subset_ID, case_study_ID, system_ID)
    
}
```

#### 3.7 Inclusion of other data

Based on the values of Boolean attributes `include_element_info` and
`include_system_info`, the script creates a joint dataframe containing
information on stored beds and on the deep-marine system to which the
beds are linked (age, global climate etc.) and/or on they parent
architectural elements (e.g. type, thickness, width etc.). The script
returns a dataframe that can be further modified in R.

``` r
if (include_element_info == TRUE & !is.null(element_data)) {
    bed_table <- left_join(bed_table, element_data, by = 'element_ID')
}
  
if (include_system_info == TRUE & !is.null(system_data)) {
    bed_table <- left_join(bed_table, system_data, by = c('system_ID', 'subset_ID',              'case_study_ID'))
}
  
return(bed_table)
```

#### 4. Example run

``` r
bed_data <- bed_classifier(facies_data = facies_data_dummy, include_vchanges = TRUE, include_lam = TRUE, 
                                                 system_data = system_data_dummy, include_system_info = TRUE, 
                                                 element_data = element_data_dummy, include_element_info = TRUE)
```

Example output (does not contain every column):

| artificial_bed_id | code | bed_thickness | gravel_fraction | laminated_fraction | log_ID | element_ID |
|----|----|----|----|----|----|----|
| 1 | S-N-m | 0.28 | 0 | 0 | 3 | 50 |
| 2 | S-N-m | 0.1 | 0 | 0 | 3 | 50 |
| 3 | S-N-m | 0.12 | 0 | 0 | 4 | 67 |
| 4 | S-N-x | 0.07 | 0 | 0.71 | 4 | 67 |
| 5 | S-N-x | 0.11 | 0 | 0.36 | 4 | 67 |
| 6 | S-F-m | 0.12 | 0 | 0.83 | 5 | 52 |
