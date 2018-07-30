#' @export
sc5mc.get_indexes_file <- function(version='bravo'){
	system.file(glue('config/indexes/{version}.yaml'), package='sc5mc')
}

sc5mc.get_config_file <- function(){
	system.file(glue('config/config.yaml'), package='sc5mc')
}

#' @export
sc5mc.dump_config_file <- function(filename){
	file.copy(sc5mc.get_config_file(), filename)
}

#' @export
sc5mc.dump_minimal_config <- function(filename){
	file.copy(system.file(glue('config/minimal_config.yaml'), package='sc5mc'), filename)
}

sc5mc.get_defaults_file <- function(){
	system.file(glue('config/defaults.yaml'), package='sc5mc')
}

# random_id <- function(){
# 	glue('{ids::proquint(n_words=1)}_{ids::random_id(bytes=2)}')
# }

#' @export
sc5mc.init_merge_pipeline <- function(config_files, workdir=getwd(), description=' ', log_file=NULL){
	if (!is.null(log_file)){
		logging::addHandler(logging::writeToFile, file=log_file)
		logging::removeHandler('basic.stdout')
		on.exit(logging::removeHandler('logging::writeToFile'))
	}
	all_conf <- map(config_files, yaml::yaml.load_file)
	plate_id <- all_conf[[1]]$plate_id	
	plate_workdir <- glue('{workdir}/{plate_id}')
	system(glue('mkdir -p {plate_workdir}'))

	conf <- all_conf[[1]]
	conf$workdir <- plate_workdir
	conf$fastq_dir <- map(all_conf, 'fastq_dir') %>% do.call(c, .)
	conf$seq_id <- map(all_conf, 'seq_id') %>% do.call(c, .)


	conf$pipeline_steps <- 'bam2smat'
	conf$description <- description

	if (has_name(all_conf[[1]], 'annotations')){
		conf$annotations <- all_conf[[1]]$annotations
	}

	readr::write_lines(yaml::as.yaml(conf), glue('{plate_workdir}/config.yaml'))			

	
	workdirs <- map_chr(all_conf, 'workdir')	

	logging::loginfo('Creating symbolic links for fastqs at %s/fastq', plate_workdir)
	raw_fastqs <- tibble(illu_dir = map(workdirs, ~ list.files(glue('{.x}/fastq/'), full.names=TRUE)) %>% do.call(c, .) ) %>% 
		purrrlyr::by_row(~ list.files(glue('{.x$illu_dir}/raw'), full.names=TRUE, pattern='.*\\.fastq\\.gz$'), .to='fn') %>%  
		unnest(fn) %>% 
		mutate(base = basename(fn), illu_dir=basename(illu_dir)) %>% 
		group_by(base) %>% 
		mutate(n = n(), i = 1:n()) %>% 
		mutate(base1 = ifelse(n > 1, paste0(i, '_', base), base)) %>% 
		ungroup() %>% 
		select(illu_dir, fn, base=base1) %>% 
		mutate(dest = glue('{plate_workdir}/fastq/{illu_dir}/raw/{base}'))

	walk(unique(raw_fastqs$illu_dir), ~ system(glue('mkdir -p {plate_workdir}/fastq/{.x}/raw/split')))
	a <- plyr::adply(raw_fastqs, 1, function(.x) system(glue('ln -sf {.x$fn} {.x$dest}')), .parallel = TRUE)
	

	fastqs <- tibble(illu_dir = map(workdirs, ~ list.files(glue('{.x}/fastq/'), full.names=TRUE)) %>% do.call(c, .) ) %>% 
		purrrlyr::by_row(~ list.files(glue('{.x$illu_dir}/raw/split'), full.names=TRUE, pattern='.*\\.fastq\\.gz$'), .to='fn') %>% 
		unnest(fn) %>% 
		mutate(base = basename(fn), illu_dir=basename(illu_dir)) %>% 
		group_by(base) %>% 
		mutate(n = n(), i = 1:n()) %>% 
		mutate(base1 = ifelse(n > 1, paste0(i, '_', base), base)) %>% 
		ungroup() %>% 
		select(illu_dir, fn, base=base1) %>% 
		mutate(dest = glue('{plate_workdir}/fastq/{illu_dir}/raw/split/{base}'))

	
	a <- plyr::adply(fastqs, 1, function(.x) system(glue('ln -sf {.x$fn} {.x$dest}')), .parallel = TRUE)	

	system(glue('mkdir -p {plate_workdir}/bam'))
	logging::loginfo('Creating symbolic links for bams at %s/bam', plate_workdir)

	bams <- tibble(fn = map(workdirs, ~ list.files(glue('{.x}/bam/'), pattern='*\\.bam$', full.names=TRUE)) %>% do.call(c, .) ) %>% 
		mutate(base = basename(fn)) %>% 
		group_by(base) %>% 
		mutate(n = n(), i = 1:n()) %>% 
		mutate(base1 = ifelse(n > 1, paste0(i, '_', base), base)) %>% 
		ungroup() %>% 
		select(fn, base=base1) %>% 
		mutate(dest = glue('{plate_workdir}/bam/{base}'))

	a <- plyr::adply(bams, 1, function(.x) system(glue('ln -sf {.x$fn} {.x$dest}')), .parallel = TRUE)
	
	cell_metadata_fn <- glue('{workdirs[1]}/cell_metadata.csv')
	if (file.exists(cell_metadata_fn)){
		file.copy(cell_metadata_fn, glue('{plate_workdir}/cell_metadata.csv'))
	}
	
	indexes_fn <- glue('{workdirs[1]}/indexes.yaml')	

	file.copy(indexes_fn, glue('{plate_workdir}/indexes.yaml'))	

	logging::loginfo('Created pipeline directories and files at %s', plate_workdir)
	logging::loginfo('Please fill config.yaml before running sc5mc.run_pipeline')	

}

#' @export
sc5mc.init_pipeline <- function(plate_id=NULL, workdir=getwd(), indexes_file = sc5mc.get_indexes_file(), config_file=NULL, seq_id=NULL, description=' ', raw_fastq_dir=NULL, log_file=NULL){	
	if (!is.null(log_file)){
		logging::addHandler(logging::writeToFile, file=log_file)
		logging::removeHandler('basic.stdout')
		on.exit(logging::removeHandler('logging::writeToFile'))
	}
	default_conf <- yaml::yaml.load_file( sc5mc.get_config_file())
	if (is.null(config_file)){		
		if (is.null(plate_id)){
			stop('Please provide plate_id')
		}
		conf <- default_conf
	}  else {
		conf <- yaml::yaml.load_file(config_file)
		if (!has_name(conf, 'plate_id')){
			if (!is.null(plate_id)){
				conf$plate_id <- plate_id
			} else {
				stop('Please provide plate_id in the config file or as a parameter')		
			}
		}
		conf <- plyr::defaults(conf, default_conf)
		plate_id <- conf$plate_id
	}

	plate_workdir <- glue('{workdir}/{plate_id}')
	system(glue('mkdir -p {plate_workdir}/fastq'))
	file.copy(indexes_file, glue('{plate_workdir}/indexes.yaml'))
	
	conf$plate_id <- plate_id
	conf$workdir <- plate_workdir
	conf$fastq_dir <- glue('{plate_workdir}/fastq')
	conf$raw_fastq_dir <- raw_fastq_dir

	if (!is.null(seq_id)){
		conf$seq_id <- seq_id	
	}
	
	conf$description <- description
	readr::write_lines(yaml::as.yaml(conf), glue('{plate_workdir}/config.yaml'))		
	logging::loginfo('Created pipeline directories and files at %s', plate_workdir)
	logging::loginfo('Please fill config.yaml before running sc5mc.run_pipeline')
}

#' @export
sc5mc.symlink_fastq <- function(raw_fastq_dir, config_file){
	conf <- yaml::yaml.load_file(config_file)
		
	illu_indexes <-  purrr::map(conf$exp_indexes, ~ .[[2]]) %>% simplify() %>% unique()
	illu_indexes <- illu_indexes[!is.na(illu_indexes) & illu_indexes != 'NA']
	workdir <- conf$workdir
	fastq_dir <- conf$fastq_dir
	orig_wd <- getwd()
	purrr::walk(illu_indexes, ~ {
			setwd(orig_wd)
			illu_fastq_dir <- glue('{fastq_dir}/{.x}/raw')
			dir.create(illu_fastq_dir, recursive=TRUE)			
			fns <- list.files(glue('{raw_fastq_dir}/{.x}/raw'), pattern='*.gz$', full.names=TRUE)
			fns <- c(fns, list.files(glue('{raw_fastq_dir}/{.x}'), pattern='*.gz$', full.names=TRUE))
			setwd(illu_fastq_dir)
			walk2(fns, basename(fns), function(from, to) file.symlink(from, to))			
		})	
	setwd(orig_wd)
}

yaml2indexes <- function(fn, indexes_file, all_indexes_file=system.file('config/indexes/indexes.csv', package = 'sc5mc'), remove_missing=TRUE, index1_len=8, index2_len=8){
    indexes <- fread(all_indexes_file) %>% as.tibble()
    indexes_yml <- yaml::yaml.load_file(indexes_file)
    conf_yml <- yaml::yaml.load_file(fn)
    yml <- c(conf_yml, indexes_yml)

    plate_id <- yml$plate_id

    all_bcds <- readr::read_csv(yml$plate) %>% gather('column', 'index', -row)
    bcds <- all_bcds

    for (i in 1:length(yml$exp_indexes)){
        idx <- names(yml$exp_indexes)[i]
        vars <- yml$exp_indexes[[i]] 

        bcds$index <- gsub(glue('_{yml$pbat2_vars[[idx]]}_'), glue('_{vars[1]}_'), bcds$index)

        if (is.list(vars[2])){
        	illu_indexes <- vars[2][[1]]
        	stopifnot(length(yml$illumina_vars[[idx]]) == length(illu_indexes))
        	for (j in 1:length(illu_indexes)){
        		bcds$index <- gsub(glue('_{yml$illumina_vars[[idx]][j]}$'), glue('_{illu_indexes[j]}'), bcds$index)		
        	}        	
        } else {
        	bcds$index <- gsub(glue('_{yml$illumina_vars[[idx]]}$'), glue('_{vars[2]}'), bcds$index)	
        }
        
    }
    
    bcds <- bcds %>% left_join(all_bcds %>% mutate(orig = TRUE), by = c("row", "column", "index")) %>% mutate(index = ifelse(is.na(orig), index, NA)) %>% select(-orig) # set indexes that are not present to NA
    
    bcds <- bcds %>% separate(index, c('index1', 'index2', 'illumina_index')) %>% mutate(plate_pos = paste0(column, row), empty = plate_pos %in% yml$empty) 
  
    bcds <- bcds %>% mutate(index2 = ifelse(!is.na(index2), paste0('revComp', index2), index2))

    
    bcds <- bcds %>% 
        rename(index = index1) %>%
        left_join(indexes, by = "index") %>% 
        rename(index1 = index, index1.seq = seq, index=index2) %>% 
        left_join(indexes, by = "index") %>% 
        rename(index2 = index, index2.seq = seq) %>%
        mutate(index1.seq = substring(index1.seq, 1, index1_len), index2.seq = substring(index2.seq, 1, index2_len))

    bcds <- bcds %>% select(row, column, plate_pos, illumina_index, index1, index1.seq, index2, index2.seq, empty)

    
    bcds <- bcds %>% left_join(purrr::imap_dfr(yml$exp_indexes, ~ tibble(batch_id = .y, index2=paste0('revComp', .x[1]), illumina_index=.x[2]) ) %>% unnest(illumina_index), by = c("illumina_index", "index2"))
    
    if (remove_missing){
    	bcds <- bcds %>% mutate(illumina_index = ifelse(illumina_index == 'NA', NA, illumina_index))
    	bcds <- bcds %>% filter(!is.na(illumina_index), !is.na(index1), !is.na(index2))
    }

    bcds <- bcds %>% arrange(column, row) %>% mutate(cell_id = paste0(plate_id, '.', 1:n()))
	bcds <- bcds %>% select(cell_id, batch_id, plate_pos, column, row, empty, illumina_index, index1, index1.seq, index2, index2.seq)

    return(bcds)
}

#' @export
sc5mc.run_pipeline <- function(config_file=NULL, workdir=NULL, indexes_file=NULL, log_file=NULL, raw_fastq_dir=NULL, defaults_file=sc5mc.get_defaults_file(), overwrite=FALSE, regions=NULL){
	if (!is.null(log_file)){
		logging::addHandler(logging::writeToFile, file=log_file)
		logging::removeHandler('basic.stdout')
		on.exit(logging::removeHandler('logging::writeToFile'))
	}

	if (is.null(config_file) && is.null(workdir)){
		loginfo('Please provide either config file or workdir')
	}
	if (is.null(config_file)){
		config_file <- glue('{workdir}/config.yaml')
	}
	conf <- yaml::yaml.load_file(config_file)

	workdir <- conf$workdir	

	if (is.null(indexes_file)){
		indexes_file <- glue('{workdir}/indexes.yaml')
	}

	if (file.exists(glue('{workdir}/defaults.yaml'))){
		loginfo('taking defaults file from "%s/defaults.yaml"', workdir)
		defaults_file <- glue('{workdir}/defaults.yaml')
	} else {
		file.copy(defaults_file, glue('{workdir}/defaults.yaml'))
	}

	cell_metadata_fn <- glue('{workdir}/cell_metadata.csv')
	cell_metadata <- yaml2indexes(config_file, indexes_file, remove_missing=TRUE)
	
	if (!is.null(conf$annotations)){
		cell_metadata <- cbind(cell_metadata, as_tibble(conf$annotations) )
	}
	loginfo("Writing cell metadata to %s", cell_metadata_fn)
	readr::write_csv(cell_metadata, cell_metadata_fn)

	raw_fastq_dir <- raw_fastq_dir %||% conf$raw_fastq_dir
	if (!is.null(raw_fastq_dir)){
		loginfo("Creating symlinks from %s to %s", raw_fastq_dir, conf$fastq_dir)
		sc5mc.symlink_fastq(raw_fastq_dir, config_file)
	}

	pipeline_dir <- glue('{workdir}/pipeline')
	system(glue('mkdir -p {pipeline_dir}'))
	defaults <- yaml::yaml.load_file(defaults_file)
	loginfo(pipeline_dir)
	gpatterns_indexes_file <- glue('{pipeline_dir}/gpatterns_config.csv')
	gpatterns_config_file <- glue('{pipeline_dir}/gpatterns_config.yaml')
	gpatterns_log_file <- glue('{pipeline_dir}/gpatterns_log')

	exp_conf <-  list(plyr::defaults(list(parallel=TRUE, config_file=gpatterns_indexes_file), conf))	
	names(exp_conf) <- conf$plate_id	
	defaults[['experiments']] <- exp_conf
	
	cell_metadata %>% mutate(experiment = conf$experiment, lib = cell_id) %>% readr::write_csv(gpatterns_indexes_file)
	readr::write_lines(yaml::as.yaml(defaults), gpatterns_config_file)
	
	file.remove(glue('{workdir}/pipeline/finished_mapping')) # always try to map unmapped files

	if (length(overwrite) == 1 && overwrite == TRUE){
		sc5mc.clean_pipeline(config_file, steps=c('demultiplexing', 'mapping', 'bam2smat'))
	} else if (any(c('demultiplexing', 'mapping', 'bam2smat') %in% overwrite)){
		sc5mc.clean_pipeline(config_file, steps=overwrite)
		overwrite <- FALSE
	}
	
	gpatterns::gpatterns.pipeline(gpatterns_config_file, log_file=gpatterns_log_file, run_dir=pipeline_dir, overwrite=overwrite)	

	smat_dir <- glue('sparse_matrices/{conf$plate_id}')
	file.symlink(smat_dir, glue('{workdir}/sc_data'))
	stats_file <- glue('sparse_matrices/{conf$plate_id}/smat_stats.tsv')
	file.symlink(stats_file, glue('{workdir}/qc_stats.tsv'))
	
	red_message("generating QC report")	
	sc5mc.pipeline_qc(workdir, raw_fastq_dir = raw_fastq_dir, regions=regions)	
	red_message("pipeline finished")
	blue_message('load the data using the following command:')
	blue_message('"scmat <- smat.load(\'{workdir}/sc_data/smat\')")')
}

#' @export
sc5mc.clean_pipeline <- function(config_file, steps){
	conf <- yaml::yaml.load_file(config_file)
	workdir <- conf$workdir	
	if ('demultiplexing' %in% steps){
		illumina_indexes <- map_chr(conf$exp_indexes, ~.x[2]) %>% unique()
		walk(illumina_indexes, ~ {
			dir <- glue('{workdir}/fastq/{.x}/raw/split')
			loginfo('removing %s', dir)
			if (dir_exists(dir)){
				dir_delete(dir)
			}			
		})
		file.remove(glue('{workdir}/pipeline/finished_demultiplexing'))
		steps <- c(steps, 'mapping', 'bam2smat')
	}
	if ('mapping' %in% steps){
		loginfo('removing %s', glue('{workdir}/bam'))	
		if (dir_exists(glue('{workdir}/bam'))){
			dir_delete(glue('{workdir}/bam'))
		}		
		file.remove(glue('{workdir}/pipeline/finished_mapping'))		
		steps <- c(steps, 'bam2smat')
	}
	if ('bam2smat' %in% steps){
		loginfo('removing %s', glue('{workdir}/sparse_matrices'))
		if (dir_exists(glue('{workdir}/sparse_matrices'))){
			dir_delete(glue('{workdir}/sparse_matrices'))
		}		
		loginfo('removing %s', glue('{workdir}/sc_data'))
		file.remove(glue('{workdir}/sc_data'))
		file.remove(glue('{workdir}/pipeline/finished_bam2smat'))
	}	

}

#' @export
sc5mc.pipeline_qc <- function(workdir, ofn = NULL, raw_fastq_dir=NULL, subtitle=NULL, regions=NULL){
	config_file <- glue('{workdir}/config.yaml')
	conf <- yaml::yaml.load_file(config_file)

	defaults <- yaml::yaml.load_file(glue('{workdir}/defaults.yaml'))
	genome_conf <- gpatterns:::apply_genome_conf(conf, defaults)

	smat_dir <- glue('{workdir}/sparse_matrices/{conf$plate_id}/smat')

	smat <- smat.load(smat_dir)
	
	if (has_name(genome_conf, 'groot')){
		gsetroot(genome_conf$groot)
		opt <- options(gmax.data.size=1e9)
		on.exit(options(opt))
		temp_cgdb <- glue('{workdir}/cgdb')
		dir.create(temp_cgdb)
		on.exit(fs::dir_delete(temp_cgdb))
		cgdb_init(temp_cgdb, intervals=smat$intervs, overwrite = TRUE)
		db <- cgdb_load(temp_cgdb)
		db <- cgdb_add_plate(db, smat, plate_name=conf$plate, verbose=FALSE)
		if (all(has_name(smat$cell_metadata, c('cell_source', 'cell_type', 'treatment')))){
			db <- db %>% group_by_cells(cell_source, cell_type, treatment)
		}
		on.exit(freemem(db))		
	} else {
		db <- NULL
	}
	  
	raw_fastq_dir <- raw_fastq_dir %||% conf$raw_fastq_dir
    
    date <- gsub('.+_20', '', conf$seq_id)   

    if (is.null(subtitle)){
    	if (all(has_name(smat$cell_metadata, c('cell_source', 'cell_type', 'treatment')))){
    		cells_annot <-  smat$cell_metadata %>% unite('name', cell_source, cell_type, treatment, sep=', ')  %>% slice(1)	%>% pull(name)
    	} else {
    		cells_annot <- NULL
    	}        

        subtitle <- glue('{conf$experiment}, {cells_annot} ({date})')                 
    }

    ofn <- ofn %||% glue('{workdir}/qc.png')

    blue_message('QC report file: {ofn}')

    if (!is.null(regions) && is.character(regions)){
    	regions <- gintervals.load(regions)
    }

    message(ofn)    
    smat <- sc5mc.qc_plot(smat, ofn=ofn, raw_reads_dir=raw_fastq_dir, subtitle=subtitle, db = db, regions=regions)    

    # if (!is.null(regions)){    	
    	# blue_message('capture QC report file: {workdir}/qc_capture.png')
		# smat <- sc5mc.qc_plot(smat, ofn=glue('{workdir}/qc_capture.png'), raw_reads_dir=raw_fastq_dir, subtitle=subtitle, db = db, regions=regions, capture_stats=TRUE)	
    # }

    rm(smat)	
    
}


#' @export
sc5mc.smat_per_experiment <- function(config, groot=NULL, prefix=NULL, workdir=tempdir(), use_sge = TRUE, name='', description='', run_commands=TRUE, log_prefix=TRUE, keep_tidy_cpgs=TRUE, load_existing=FALSE, cell_metadata=NULL, single_cell=TRUE, step_file=NULL, ...){	
	orig_config <- config

	if (has_name(config, 'experiment_type')){
		config <- config %>% filter(experiment_type == 'single_cell')
		if (nrow(config) == 0){
			logging::loginfo('no single cell experiments')
			return(config)
		}
	}

	if (has_name(config, 'workdir')){
		workdir <- config$workdir
	}

	if (!is.null(cell_metadata)){
		cell_metadata <- glue(cell_metadata)
	}

	if (has_name(config, 'bams')){
		config <- config %>% select(experiment, lib, bams, workdir, groot) %>% unnest(bams) %>% rename(bam = bam_file)	
	} else {
		config <- config %>% select(experiment, lib, bam=bam_file, workdir, groot) 
	}	
	
	exp_cfg <- config %>% distinct(experiment, .keep_all=TRUE) %>% mutate(prefix=glue(prefix), name=glue(name), description=glue(description))

	run_per_exp <- function(exp_name){		
		logging::loginfo('running for experiment: %s', exp_name)

		cmd_cfg <- exp_cfg %>% filter(experiment == exp_name) 
		exp_config <- config %>% filter(experiment == exp_name)
		
		gpatterns:::do.call_tibble(sc5mc:::smat.from_bams, cmd_cfg, list(metadata=exp_config, use_sge=use_sge, keep_tidy_cpgs=keep_tidy_cpgs, load_existing=load_existing, cell_metadata=cell_metadata, single_cell=single_cell))
	}
	
	if (run_commands){
		res <- map(unique(exp_cfg$experiment), ~ run_per_exp(.x))	
	}	
	
	config <- orig_config %>% left_join(exp_cfg %>% select(experiment, smat_prefix=prefix, smat_name=name, smat_description=description), by='experiment')
	return(config)
}