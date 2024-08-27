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
				ouwie_run <- OUwie::OUwie(phy=phy, data=df, model = model, mserr = fog_choice)
				ouwie_run$chosen_model <- model
				ouwie_run$fog <- fog_choice
				ouwie_run$dataset_name <- dataset_name
				results[[length(results)+1]] <- ouwie_run
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
	geospiza_cleaned<- geiger::treedata(geospiza$phy, geospiza$dat)
	geospiza_cleaned$dat <- geospiza_cleaned$dat[,'wingL']
	all_data[['geospiza']] <- geospiza_cleaned
	data(primates, package="geiger")
	all_data[['primates']] <- primates
	return(all_data)
}

summarize_results <- function(all_results) {
	all_df <- data.frame()
	for (i in sequence(length(all_results))) {
		ouwie_run <- all_results[[i]]
		ouwie_df <- data.frame(dataset=ouwie_run$dataset_name, model=ouwie_run$chosen_model, fog=ouwie_run$fog, data_mean = mean(ouwie_run$data[,2]), AICc=ouwie_run$AICc, alpha=ouwie_run$solution['alpha', 1], sigma.sq=ouwie_run$solution['sigma.sq', 1], theta1=NA, theta2=NA, sigma.sq.me=0)
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
	return(all_df)
}