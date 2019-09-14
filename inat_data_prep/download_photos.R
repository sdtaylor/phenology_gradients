library(tidyverse)
source('config.R')


download_images = function(metadata, save_dir){
  metadata$local_image_file = NA
  for(i in 1:nrow(metadata)){
    local_filename = paste0(save_dir, metadata$inat_id[i], '.jpg')
    download_url = metadata$image_url[i]
    # Get the original full size photo
    download_url = stringr::str_replace(download_url,'medium','original')
    attempt1 = try(download.file(download_url, local_filename))
    if(class(attempt1)=='try-error'){
      # If it doesn't work try again with the medium res link
      download_url = metadata$image_url[i]
      attempt2 = try(download.file(download_url, local_filename))
    }
    metadata$local_image_file[i] = local_filename
  }
  return(metadata)
}

#Drop images that aren't hosted on the inaturalist domain
drop_non_inat_photos = function(metadata){
  metadata = metadata[grepl('inaturalist', metadata$image_url),]
  return(metadata)
}


#######################################

rudbeckia = read_csv(paste0(data_folder, 'inat_rudbeckia_hirta.csv')) %>%
  filter(quality_grade=='research') %>%
  rename(inat_id = id)

rudbeckia = download_images(rudbeckia, save_dir = paste0(data_folder, 'photos/rudbeckia/'))

maianthemum = read_csv(paste0(data_folder, 'inat_maianthemum_canadense.csv')) %>%
  filter(quality_grade=='research') %>%
  rename(inat_id = id)

maianthemum = download_images(maianthemum, save_dir = paste0(data_folder, 'photos/maianthemum/'))


all_metadata = rudbeckia %>%
  bind_rows(maianthemum) %>%
  select(inat_id, date = observed_on, latitude, longitude, positional_accuracy, scientific_name)


write_csv(all_metadata, paste0(data_folder, 'inaturalist_metadata.csv'))
