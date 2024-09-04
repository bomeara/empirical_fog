run_datasets <- function(geiger_data) {
	dataset_names <- names(geiger_data)
	results <- list()
	for (dataset_name in dataset_names) {
		print(dataset_name)
		dataset <- geiger_data[[dataset_name]]
		phy <- dataset$phy
		phy$node.label <- rep(1, ape::Nnode(phy))
		dat <- dataset$dat
		df <- data.frame(names(dat), rep(1, length(dat)), dat)
		rownames(df) <- NULL
		models <- c("BM1", "OU1")
		fog_possibilities <- c("none", "estimate")
		for (model in models) {
			for (fog_choice in fog_possibilities) {
				ouwie_run <- OUwie::OUwie(phy=phy, data=df, model = model, mserr = fog_choice, algorithm="three.point")
				ouwie_run$chosen_model <- model
				ouwie_run$fog <- fog_choice
				ouwie_run$dataset_name <- dataset_name
				results[[length(results)+1]] <- ouwie_run
			}
		}
	}
	return(results)			
}

# After O'Meara et al. 2006 eq 1
expected_disparity_excluding_sigma_sq <- function(phy) {
	C_matrix <- ape::vcv(phy)
	ntax <- ape::Ntip(phy)
	result <- (1/ntax)*(sum(diag(C_matrix))) - (1/(ntax^2))*sum(C_matrix)
	return(result)
}


run_datasets_geiger <- function(geiger_data) {
	dataset_names <- names(geiger_data)
	results <- list()
	for (dataset_name in dataset_names) {
		print(dataset_name)
		dataset <- geiger_data[[dataset_name]]
		phy <- dataset$phy
		dat <- dataset$dat
		models <- c("BM", "OU", "EB", "kappa", "delta", "white")
		fog_possibilities <- c("none", "estimate")	
		for (model in models) {
			for (fog_choice in fog_possibilities) {
				if(model=="white" & fog_choice=="estimate") {
					next
				}
				print(paste0("Running ", model, " with ", fog_choice, " for ", dataset_name))
				geiger_run <- geiger::fitContinuous(phy, dat, model=model, SE=ifelse(fog_choice=="none", 0, NA))	
				geiger_run$chosen_model <- model
				geiger_run$fog <- fog_choice
				geiger_run$dataset_name <- dataset_name
				geiger_run$other_param <- NA
				
				if(model=="BM") {
					geiger_run$expected_tip_variance_along_tree <- geiger_run$opt$sigsq*expected_disparity_excluding_sigma_sq(phy)
				}
				if(model=="OU") {
					transformed_tree <- rescale(phy, model="OU", geiger_run$opt$alpha)
					geiger_run$other_param <- geiger_run$opt$alpha
					geiger_run$expected_tip_variance_along_tree <- geiger_run$opt$sigsq*expected_disparity_excluding_sigma_sq(transformed_tree)
				}
				if(model=="EB") {
					transformed_tree <- rescale(phy, model="EB", geiger_run$opt$a)
					geiger_run$other_param <- geiger_run$opt$a
					geiger_run$expected_tip_variance_along_tree <- geiger_run$opt$sigsq*expected_disparity_excluding_sigma_sq(transformed_tree)
				}
				if(model=="kappa") {
					transformed_tree <- rescale(phy, model="kappa", geiger_run$opt$kappa)
					geiger_run$other_param <- geiger_run$opt$kappa
					geiger_run$expected_tip_variance_along_tree <- geiger_run$opt$sigsq*expected_disparity_excluding_sigma_sq(transformed_tree)
				}
				if(model=="delta") {
					transformed_tree <- rescale(phy, model="delta", geiger_run$opt$delta)
					geiger_run$other_param <- geiger_run$opt$delta
					geiger_run$expected_tip_variance_along_tree <- geiger_run$opt$sigsq*expected_disparity_excluding_sigma_sq(transformed_tree)
				}
				if(model=="white") {
					geiger_run$expected_tip_variance_along_tree <- 0
					# so as to better align with the other models
					geiger_run$opt$SE <- geiger_run$opt$sigsq
					geiger_run$opt$sigsq <- 0
				}
				
				results[[length(results)+1]] <- geiger_run
			}
		}
	}
	return(results)
}

get_data <- function() {
	all_data <- list()
	data(caniformia, package="geiger")
	all_data[['caniformia']] <- caniformia
	data(carnivores, package="geiger")
	carnivore_data <- carnivores$dat[, 'mean']
	names(carnivore_data) <- rownames(carnivores$dat)
	carnivores$dat <- carnivore_data
	all_data[['carnivores']] <- carnivores
	data(geospiza, package="geiger")
	geospiza_cleaned<- geiger::treedata(geospiza$phy, geospiza$dat, warnings=FALSE)
	geospiza_cleaned$dat <- geospiza_cleaned$dat[,'wingL']
	all_data[['geospiza']] <- geospiza_cleaned
	data(primates, package="geiger")
	all_data[['primates']] <- primates
	return(all_data)
}

get_Alencar_et_al_data <- function(all_data) {
	traits <- as.data.frame(readxl::read_xlsx("data/AlencarEtAl/Table_S2_JAN2024.xlsx"))
	phy <- ape::read.tree("data/AlencarEtAl/PleuroTreeToWork.tre")
	traits$tempK <- traits$temp + 273.15
	traits$log_tempK <- log(traits$tempK)
	traits$log_svl <- log(traits$female.svl)
	traits$log_topo_complexity <- log(traits$topo.complex)
	traits$log_precip <- log(traits$prec)
	focal_traits <- c('log_tempK', 'log_svl', 'log_topo_complexity', 'log_precip')
	for (trait_name in focal_traits) {
		focal_trait <- traits[[trait_name]]
		names(focal_trait) <- traits$Species
		focal_trait <- focal_trait[!is.na(focal_trait)]
		cleaned <- geiger::treedata(phy, focal_trait, warnings=FALSE)
		all_data[[paste0("AlencarEtAl_", trait_name)]] <- list(phy=cleaned$phy, dat=cleaned$dat[,1])
	}
	return(all_data)
}

summarize_ouwie_results <- function(all_results) {
	all_df <- data.frame()
	for (i in sequence(length(all_results))) {
		ouwie_run <- all_results[[i]]
		ouwie_df <- data.frame(dataset=ouwie_run$dataset_name, model=ouwie_run$chosen_model, fog=ouwie_run$fog, data_mean = mean(ouwie_run$data[,2]), AICc=ouwie_run$AICc, alpha=ouwie_run$solution['alpha', 1], sigma.sq=ouwie_run$solution['sigma.sq', 1], theta1=NA, theta2=NA, sigma.sq.me=0, ntax=ape::Ntip(ouwie_run$phy), treeheight=max(ape::branching.times(ouwie_run$phy)), loglik=ouwie_run$loglik)
		try(ouwie_df$theta1 <- ouwie_run$theta[1,1], silent=TRUE)
		try(ouwie_df$theta2 <- ouwie_run$theta[2,1], silent=TRUE)
		if(!is.null(ouwie_run$sigma.sq.me)) {
			ouwie_df$sigma.sq.me <- ouwie_run$sigma.sq.me
		}
		ouwie_df$tip_fog_as_sd <- sqrt(ouwie_df$sigma.sq.me)
		ouwie_df$tip_fog_as_fraction_trait_mean <- ouwie_df$tip_fog_as_sd  / ouwie_df$data_mean
		all_df <- rbind(all_df, ouwie_df)	
	}
	all_df <- all_df |> arrange(dataset, AICc)
	all_df$program <- "ouwie"
	return(all_df)
}

summarize_geiger_results <- function(all_geiger_results) {
	all_df <- data.frame()
	for (i in sequence(length(all_geiger_results))) {
		geiger_run <- all_geiger_results[[i]]
		geiger_df <- data.frame(dataset=geiger_run$dataset_name, model=geiger_run$chosen_model, fog=geiger_run$fog,  AICc=geiger_run$opt$aicc, loglik=geiger_run$opt$lnL, sigma.sq=geiger_run$opt$sigsq, sigma.sq.me=0, expected_tip_variance_along_tree=NA, other_param=geiger_run$other_param)

		if(!is.null(geiger_run$opt$SE)) {
			geiger_df$sigma.sq.me <- geiger_run$opt$SE
		}
		if(!is.null(geiger_run$expected_tip_variance_along_tree)) {
			geiger_df$expected_tip_variance_along_tree <- geiger_run$expected_tip_variance_along_tree
		}
		geiger_df$tip_fog_as_sd <- sqrt(geiger_df$sigma.sq.me)
		all_df <- rbind(all_df, geiger_df)	
	}
	all_df <- all_df |> arrange(dataset, AICc)
	all_df$program <- "geiger"
	return(all_df)
}

organize_final_df <- function(all_df) {
	datasets <- unique(all_df$dataset)
	final_df <- data.frame()
	for (focal_dataset in datasets) {
		focal_df <- all_df[all_df$dataset==focal_dataset,]
		focal_df$delta_AICc <- focal_df$AICc-min(focal_df$AICc)
		focal_df$data_mean <- mean(focal_df$data_mean, na.rm=TRUE)
		focal_df$ntax <- mean(focal_df$ntax, na.rm=TRUE)
		focal_df$treeheight <- mean(focal_df$treeheight, na.rm=TRUE)
		final_df <- rbind(final_df, focal_df)
	}
	final_df$tip_fog_as_fraction_trait_mean <- final_df$tip_fog_as_sd / final_df$data_mean
	final_df$fog[final_df$model=="white"] <- "only"

	final_df$tip_variance_divided_by_tree_plus_tip_variance <- final_df$sigma.sq.me / (final_df$expected_tip_variance_along_tree + final_df$sigma.sq.me)
	final_df <- final_df |> arrange(dataset, AICc) |> select(dataset, model, fog, program, delta_AICc, tip_fog_as_fraction_trait_mean, tip_variance_divided_by_tree_plus_tip_variance, sigma.sq, sigma.sq.me, other_param, ntax, treeheight)
	return(final_df)
}