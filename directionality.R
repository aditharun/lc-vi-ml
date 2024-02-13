

getShap <- function(merged, class1varname, model_seq_xgb){

	label.merge <- merged %>% mutate(status = ifelse(status == class1varname, 1, 0)) %>% pull(status)

	fparams <- model_seq_xgb$finalModel$tuneValue %>% as_tibble() 

	ft_cols <- attr(model_seq_xgb$terms, "term.labels") %>% as.character() %>% gsub("`", "", .)

	merged <- merged %>% select(c(ft_cols, "status"))

	bst <- xgboost::xgboost(as.matrix( (merged %>% select(-status))),label.merge , nrounds = fparams$nrounds, eta = fparams$eta, max_depth = fparams$max_depth, subsample = fparams$subsample, gamma = fparams$gamma, colsample_bytree = fparams$colsample_bytree, min_child_weight = fparams$min_child_weight, objective = "binary:logistic", nthread = 2, verbose = 0)

	shap.df <- xgboost::xgb.plot.shap(data = as.matrix((merged %>% mutate(status = ifelse(status == class1varname, 1, 0)) %>% select(-status))), model = bst, top_n = 90, plot=FALSE)$shap_contr %>% as_tibble() %>% mutate_all(~ abs(.)) %>% colMeans()  %>% as.data.frame() 


	shap.df <- shap.df %>% mutate(var = rownames(shap.df)) %>% as_tibble() %>% magrittr::set_colnames(c("mean_abs_shap", "var")) %>% arrange(desc(mean_abs_shap)) %>% mutate(id = 1:n()) 


	true_names <- papervarnames_file %>% read_csv() %>% select(2,3)

	df_shap <- xgboost::xgb.plot.shap(data = as.matrix((merged %>% mutate(status = ifelse(status == class1varname, 1, 0)) %>% select(-status))), model = bst, top_n = 90, plot=FALSE)$shap_contr %>% as_tibble() 

	dff <- df_shap %>% mutate_all(~ abs(.)) %>% colMeans() %>% as.data.frame() %>% arrange(desc(.)) 

	dff$varnames <- rownames(dff)

	dff <- dff %>% as_tibble() %>% left_join(true_names, by=c("varnames"="new_names")) 

	colname.order <- dff$proper_names

	#+ geom_quasirandom(dodge.width = .65, size = 2.25)

	plt <- df_shap %>% mutate(idx = 1:n(), status = merged$status) %>% mutate(status = ifelse(status == "LC", "LC", "PVS")) %>% pivot_longer(-c(idx, status)) %>% left_join(true_names, by=c("name"="new_names")) %>% mutate(proper_names = factor(proper_names, levels = colname.order)) %>% group_by(proper_names, status) %>% summarize(m = mean(value)) %>% ggplot(aes(x=proper_names, y=m, color=status, fill=status)) + geom_col() + coord_flip() + theme(panel.grid = element_blank(), panel.grid.major.y = element_line(color = "grey70", size = 0.1), panel.background = element_rect(color = "black", fill = "transparent")) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.3) + xlab("") + ylab("Mean Shapley Value") + theme(legend.title = element_blank())

	return(list(shap.df, plt))

}


xgbModel <- xgb_stats(merged, "status", class1_variablename, cv_number, cv_repeat)

	results_xgb_shap <- xgbModel$stats 
	varimp_store_xgb_shap <- list(getShap(merged, class1_variablename, xgbModel$model_seq))

	startvarlen <- floor(nrow(varimp_store_xgb_shap[[1]][[1]])/5)*5 - 5

	while (startvarlen >= 1){

		fts <- c(varimp_store_xgb_shap[[1]][[1]] %>% dplyr::slice(1: startvarlen) %>% pull(var), "status")

		fts <- fts %>% gsub("`", "", .)

		gs <- xgb_stats((merged %>% select(fts)), "status", class1_variablename, cv_number, cv_repeat)

		results_xgb_shap <- rbind(results_xgb_shap, gs$stats)

		if (startvarlen == 1){
			tmp.varimp <- data.frame(mean_abs_shap = 1, var = fts[fts!="status"], id=1)
		} else{
			tmp.varimp <- getShap(merged, class1_variablename, gs$model_seq)
		}

		varimp_store_xgb_shap[[length(varimp_store_xgb_shap) + 1]] <- tmp.varimp

		if (startvarlen <= 20){

			startvarlen <- startvarlen - 1
		} else{
				startvarlen <- startvarlen - 5
		}

	}

#