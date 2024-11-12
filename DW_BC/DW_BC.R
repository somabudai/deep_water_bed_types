library(RMariaDB)
library(ggplot2)
library(dplyr)
library(stringr)

#READ DUMMY DATA

facies_data_dummy <- read.csv('facies_data_dummy.csv', header = TRUE)
system_data_dummy <- read.csv('system_data_dummy.csv', header = TRUE)
element_data_dummy <-  read.csv('element_data_dummy.csv', header = TRUE)


#FUNCTION

#' Bed classifier
#'
#' @param facies_data Input dataframe that contains facies data collected from
#'                    sedimentary logs. 
#' @param include_vchanges A boolean input value which indicates of the
#'                    output dataframe should include the grain size trend
#'                    in the bed-type classification.
#' @param include_lam A boolean input value which indicates if the
#'                    output dataframe should include the lamination proportion
#'                    in the bed-type classification.
#' @param system_data Optional input dataframe that contains data on the deep-water system
#'                    in which the analysed facies was deposited.
#' @param include_system_info A boolean input value which indicates if the
#'                    output dataframe should include information on the host system.
#' @param element_data Optional input dataframe that contains data on the deep-water elements
#'                    that host the analysed beds.
#' @param include_element_info A boolean input value which indicates if the
#'                    output dataframe should include information parent element.
#'
#' @return A dataframe that contains the classified beds with optional information
#'        on the host system and/or the parent archtiectural element.

#' @examples
#' 
#' bed_classifier(facies_data = facies_data_dummy, include_vchanges = TRUE, include_lam = TRUE, 
#' system_data = system_data_dummy, include_system_info = TRUE, 
#' element_data = element_data_ddummy, include_element_info = TRUE)
#' 
bed_classifier <- function (facies_data = NULL, include_vchanges = FALSE, include_lam = FALSE, 
                            system_data = NULL, include_system_info = FALSE, 
                            element_data = NULL, include_element_info = FALSE) {
  
  if (is.null(facies_data)) {
    stop('Input facies data is missing!')
  }
  
  #identify beds
  facies_data_bedboundary <- facies_data %>% 
    filter(!is.na(trans_from_below_event_b)) %>% 
    filter(facies_type %in% c('_G','_S','gS','sG','(g)S', 'mS', 'sM')) %>%
    filter(thickness_type %in% c('true (maximum)','true (not maximum)','apparent')) 
  
  facies_data_bednr <- facies_data_bedboundary  %>% 
    mutate(new_bed = case_when(data_order != lag(data_order)+1 | trans_from_below_event_b == "Y" ~ 1, TRUE ~ 0))%>% 
    mutate(artificial_bed_id = cumsum(new_bed))
  
  has_y <- facies_data_bednr %>% 
    filter(trans_from_below_event_b == "Y") %>% 
    group_by(artificial_bed_id)  %>% 
    count(trans_from_below_event_b) %>% 
    dplyr::select(artificial_bed_id,n)
  
  beds_raw <- left_join(facies_data_bednr, has_y, by = "artificial_bed_id") %>% filter(!is.na(n))
  
  #calculate bed thickness
  
  beds_raw <- beds_raw %>% 
    group_by(artificial_bed_id) %>% 
    mutate(bed_thickness = sum(thickness)) %>% 
    mutate(start_facies = min(facies_ID)) %>%
    mutate(end_facies = max(facies_ID)) %>% 
    ungroup()
  
  #grain size
  
  beds_gs <- beds_raw %>%
    group_by(artificial_bed_id) %>%
    mutate(gravel_thickness = sum(ifelse(facies_type %in% c('sG','_G'), thickness, 0 ))) %>%
    mutate(gsand_thickness = sum(ifelse(facies_type %in% c('gS','(g)S'), thickness, 0 ))) %>%
    mutate(sand_thickness = sum(ifelse(facies_type %in% c('_S'), thickness, 0 ))) %>%
    mutate(heterolithic_thickness = sum(ifelse(facies_type %in% c('mS', 'sM'), thickness, 0 ))) %>%
    mutate(gravel_fraction = round((gravel_thickness/bed_thickness),2)) %>%
    mutate(gsand_fraction = round((gsand_thickness/bed_thickness),2)) %>%
    mutate(sand_fraction = round((sand_thickness/bed_thickness),2)) %>%
    mutate(heterolithic_fraction = round((heterolithic_thickness/bed_thickness),2)) %>%
    ungroup() %>%
    mutate(gs_category = case_when(heterolithic_fraction > 0 ~ 'SM',
                                   gravel_fraction == 1 ~ 'G', 
                                   gravel_fraction == '0' ~ 'S', 
                                   gravel_fraction < 1 & gravel_fraction >= 0.5 ~ 'sG', 
                                   gravel_fraction < 0.5 & gravel_fraction > 0 ~ 'gS')) %>%
    distinct(artificial_bed_id, .keep_all = TRUE)
  
  #vertical changes
  if (include_vchanges == TRUE) {
    
    
    combined_gs_values <- data.frame(
      combined_gs = c('G',	'S',	'Svery coarse',	'Scoarse',	'Smedium',	'Sfine',	'Svery fine',	
                      'MS',	'MSvery coarse',	'MScoarse',	'MSmedium',	'MSfine',	'MSvery fine',	
                      'SM',	'SMvery coarse',	'SMcoarse',	'SMmedium',	'SMfine',	'SMvery fine'),
      gs_value = c(17, 14,	16,	15,	14,	13,	12,	
                   8,	10,	9,	8,	7,	6,	3,	
                   5,	4,	3,	2,	1)
    )
    
    no_sand_gs <- beds_raw %>% dplyr::select(artificial_bed_id, facies_type, grain_size_sand) %>%
      mutate(simple_gs = case_when(facies_type %in% c('_G', 'sG')~ 'G', facies_type %in% c('(g)S','_S','gS') ~ 'S', 
                                   facies_type %in% c('mS') ~ 'MS', facies_type %in% c('sM') ~ 'SM')) %>%
      mutate(has_sand_gs = case_when(simple_gs %in% c('S','MS') & is.na(grain_size_sand) ~ 'N', TRUE ~ 'Y')) %>%
      filter(has_sand_gs == 'N') %>% distinct(artificial_bed_id) %>% pull(artificial_bed_id) %>% as.list()
    
    combined_vchanges <- beds_raw %>% 
      mutate(simple_gs = case_when(facies_type %in% c('_G', 'sG')~ 'G', facies_type %in% c('(g)S','_S','gS') ~ 'S', 
                                   facies_type %in% c('mS') ~ 'MS', facies_type %in% c('sM') ~ 'SM')) %>%
      mutate(general_structure = case_when(is.na(general_structure) ~ 'nd', TRUE ~ general_structure )) %>% 
      mutate(lamination_type = case_when(is.na(lamination_type) ~ 'nd', TRUE ~ lamination_type )) %>%
      mutate(combined_gs = case_when(simple_gs == 'G' ~ simple_gs,
                                     is.na(grain_size_sand) ~ simple_gs, 
                                     !is.na(grain_size_sand) ~ paste(simple_gs,grain_size_sand, sep = ''))) %>%
      dplyr::select(facies_ID,log_ID,data_order,artificial_bed_id,trans_from_below,trans_from_below_event_b, 
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
      utate(code2 = lam_category) %>%
      dplyr::select(artificial_bed_id, code, code2,
                    bed_thickness, gravel_fraction, heterolithic_fraction, laminated_fraction,
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
                    bed_thickness, gravel_fraction, heterolithic_fraction, laminated_fraction,
                    gs_category, combined_trend, lam_category,
                    log_ID, element_ID, subset_ID, case_study_ID, system_ID)
    
  }
  
  
  if (include_element_info == TRUE & !is.null(element_data)) {
    bed_table <- left_join(bed_table, element_data, by = 'element_ID')
  }
  
  if (include_system_info == TRUE & !is.null(system_data)) {
    bed_table <- left_join(bed_table, system_data, by = c('system_ID', 'subset_ID', 'case_study_ID'))
  }
  
  return(bed_table)
  
}



classification_test <- bed_classifier(facies_data = facies_data_dummy, include_vchanges = TRUE, include_lam = TRUE, 
                                                 system_data = system_data_dummy, include_system_info = TRUE, 
                                                 element_data = element_data_dummy, include_element_info = TRUE)
