library(data.table)
library(magrittr)

# all the data!

ohie_data <- list.files("data/interim/OHIE/", full.names = T) %>%
  lapply(fread) %>%
  Reduce(function(x, y) merge(x, y, by = "person_id", all.x = T), .)

# rename variables based on prepare_data.do
name_tbl <- fread(text =
"	ohp_all_mo_matchn_30sep2009	ohp_all_mo_admin
	ohp_all_mo_firstn_30sep2009	ohp_all_mo_survey
	ohp_all_ever_firstn_30sep2009	ohp_all_ever_survey
	ohp_all_ever_firstn_survey0m	ohp_all_ever_survey0m
	ohp_all_ever_firstn_survey6m	ohp_all_ever_survey6m
	ohp_all_end_30sep2009	ohp_all_end_admin
	ohp_all_mo_firstn_survey0m	ohp_all_mo_survey0m
	ohp_all_mo_firstn_survey6m	ohp_all_mo_survey6m
	ohp_std_ever_matchn_30sep2009	ohp_std_ever_admin
	ohp_std_ever_firstn_30sep2009	ohp_std_ever_survey
	tanf_ever_matchn_30sep2009	postn_tanf_bin
	tanf_ever_prenotify07	prenany_tanf_bin
	tanf_ever_firstn_survey12m	postn_survey12m_tanf_bin
	tanf_ever_presurvey12m	pren_survey12m_tanf_bin
	tanf_tot_hh_firstn_survey12m	postn_survey12m_tanf_hh_amt
	tanf_tot_hh_30sep2009	postn_tanf_hh_amt
	tanf_tot_hh_prenotify07	prenany_tanf_hh_amt
	tanf_tot_hh_presurvey12m	pren_survey12m_tanf_hh_amt
	snap_ever_matchn_30sep2009	postn_snap_bin
	snap_ever_presurvey12m	pren_survey12m_snap_bin
	snap_tot_hh_30sep2009	postn_snap_hh_amt
	snap_tot_hh_prenotify07	prenany_snap_hh_amt
	snap_tot_hh_presurvey12m	pren_survey12m_snap_hh_amt
	snap_ever_prenotify07	prenany_snap_bin
	snap_ever_firstn_survey12m	postn_survey12m_snap_bin
	snap_tot_hh_firstn_survey12m	postn_survey12m_snap_hh_amt
	zip_msa_list	zip_msa
	wave_survey0m	draw
	wave_survey12m	draw_survey_12m ",
fill = T
)[,.(OldNames = V2, NewNames = V3)]

# rename data
setnames(ohie_data, name_tbl$OldNames, name_tbl$NewNames)

# write to large csv file
fwrite(ohie_data, "data/interim/OHIE/cleaned_ohie_data.csv")

