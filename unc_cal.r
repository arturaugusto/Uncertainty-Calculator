main <- function(object){
	require("numDeriv")
	# Remove the spare row
	object$table_data <- head(object$table_data, -1)
	w_s <- function (ui, df, ci = rep(1, length(ui)), uc = sqrt(sum((ci*ui)^2))){
		(uc^4)/sum(((ci*ui)^4)/df)
	}

	distributions <- list(Rect. = sqrt(3), Norm. = 2, Triang. = sqrt(6), U = sqrt(2))

	out_table_data <- data.frame( VI = numeric(), VC = numeric(), e = numeric(), MPE =  numeric(), U = numeric(), k = numeric(), veff = numeric() )

	details <- c()

	# Iterate over data
	data_nrows <- NROW(object$table_data)
	for(i in seq(data_nrows)){
		out <- tryCatch({
			# quantities means
			x <- c()
			# quantities uncertanties
			u <- c()
			# max permissible erro
			mpe <- c()
			description <- c()
			nu <- c()
			var_names <- c()
			# Create env for each row
			row_env <- new.env()
			# Set influence quantities values to new env
			if( length(object$value$influence_quantities) ){
				for(iq_name in unlist(object$value$influence_quantities["name"]) ){
					if(iq_name %in% object$table_data[i,] ){
						iq_value <- as.numeric(object$table_data[i,][iq_name][[1]])
						assign(iq_name, iq_value, envir = row_env)
					}
				}
			}

			data_row <- (object$table_data)[i,]
			n_row_vars <- NROW(object$value$variables)
			for(j in seq(n_row_vars)){
				# Procede if its not invisible
				var_name <- object$value$variables$name[[j]]
				if(!"Invisible" %in% object$value$variables$kind[[j]]){

					# Match the title for replications of this variable using regex
					match <- ls( object$table_data[1,], pattern = paste("^(", var_name, "){1}(\\s){1}([0-9])+$", sep="") )

					# Select input data, and convert to numeric
					#inputs <- lapply(data_row[,match], function(x){as.numeric(x)})
					inputs <- as.numeric(data_row[,match])

					# Filter only numeric and convert to vector
					numeric_input <- Filter(function(x) !is.na(x) && is.numeric(x), inputs)

					# Number of values
					n <- length(numeric_input)

					# Readout for quantit is the mean
					readout <- mean(numeric_input)
					assign(var_name, readout, envir = row_env)

					# Standard desviation
					s <- sd(numeric_input)

					# Uncertanty associated. For this kind of uncertanty, the distribution is "normal"
					s_r <- s/sqrt(n)

					# Store data on lists
					x <- c(x, readout)
					u <- c(u, s_r)
					description <- c(description, paste("Readout", " (", var_name,")", sep="") )
					var_names <- c(var_names, var_name)
					nu <- c(nu, n-1)
				}
			}
			# Iterate again, now with all enviroment var set
			for(j in seq(n_row_vars)){
				var_name <- object$value$variables$name[[j]]
				# Copy env
				u_eval_env <- new.env()
				for(n in ls(row_env, all.names=TRUE)) assign(n, get(n, row_env), u_eval_env)
				
				# Find which is the index on lookup table for current var and table_data line
				lookup <- object$lookup[which((object$lookup$row_index == (i-1)) & (object$lookup$var == var_name)),]
				
				# Note that the data on lookup starts on 0, so we need to add 1
				snippet <- object$asset_snippets$snippets[lookup$snippet_index+1,]
				range <- snippet$value$ranges[[1]][lookup$range_index+1,]
				assign("range", range, envir = u_eval_env)
				# Assign some aliases
				assign("range_start", range$limits$start, envir = u_eval_env)
				assign("range_end", range$limits$end, envir = u_eval_env)
				assign("full_scale", range$limits$fullscale, envir = u_eval_env)

				# Assign some helper variables
				var_name_index <- which(object$value$variables$name == var_name)
				is_UUT <- "UUT" == object$value$variables$kind[var_name_index]
				print("Measurement" == range$kind)
				assign("is_UUT", is_UUT, envir = u_eval_env)
				assign("is_meas", "Measurement" == range$kind, envir = u_eval_env)
				assign("is_source", "Source" == range$kind, envir = u_eval_env)
				assign("is_fixed", "Fixed" == range$kind, envir = u_eval_env)

				# Define readout
				if(length(ls(u_eval_env, pattern=paste("^", var_name, "$", sep="")))){
					readout <- get(var_name, envir = u_eval_env)
					assign("readout", readout, envir = u_eval_env)
				}

				# Evaluate function that CAN overwrite or set the readout
				if(!is.null(range$nominal_value)){
					with(u_eval_env, eval(parse(text=range$nominal_value)))
				}
				
				# Define var value
				if(length(ls(u_eval_env, pattern=paste("^readout$", sep="")))){
					readout <- get("readout", envir = u_eval_env)
					assign(var_name, readout, envir = row_env)
				}

				# Check for reclassifications if reclass
				# index exists and var isnt of UUT kind
				if(!is.na(lookup$reklass_index) && ("UUT" != object$value$variables[1,]$kind) ){
					uncertanty_set <- range$reclassifications[[1]][lookup$reklass_index+1,]$uncertainties
					# Evaluate function for correction
					with(u_eval_env, eval(parse(text=range$reclassifications[[1]][lookup$reklass_index+1,]$correction)))
				}else{
					uncertanty_set <- range$uncertainties
				}
				uncertanty_set <- uncertanty_set[[1]]
				
				# Assign to var on ROW_env, the readout with corrections
				if( length(ls(u_eval_env, pattern="^correct_readout$")) ){
					assign(var_name, u_eval_env$correct_readout, envir = row_env)
				}
				# Iterate over uncertanties
				n_uncertainties <- nrow(uncertanty_set)
				for(k in seq(n_uncertainties)){
					# First, set to 0 variables on env
					assign("u", 0, u_eval_env)
					assign("mpe", 0, u_eval_env)
					with(u_eval_env, eval(parse(text=uncertanty_set[k,]$formula)))
					# Store data on lists
					x <- c(x, 0)
					u <- c(u, get("u", envir = u_eval_env)/as.numeric(distributions[uncertanty_set[k,]$distribution]) )
					mpe <- c(mpe, get("mpe", envir = u_eval_env) )
					description <- c(description, paste(uncertanty_set[k,]$description, " (", var_name,")", sep="") )
					var_names <- c(var_names, var_name)
					nu <- c(nu, 9999)

				}
			}
			expr <- parse(text=object$value$formula)
			coefs <- c(rep(1,length(var_names)))
			# Iterate again, now with all enviroment var set
			for(j in seq(n_row_vars)){
				var_name <- object$value$variables$name[[j]]
				arg_list <- as.list(row_env)

				# Define derivation point as the current var value
				deriv_point <- as.numeric(arg_list[var_name])

				# Remove derivation point to define the function to be derivated
				arg_list[var_name] <- NULL

				fdx <- function(x){
					expr <- parse(text=object$value$formula)
					arg_list[var_name] <- x
					for(k in expr){
						expr_sub <- do.call("substitute", list(k, arg_list))
						y <- eval(expr_sub)
					}
					return(y)
				}

				coef <- grad(fdx, deriv_point)
				
				for(coef_index in which(var_names == var_name)){
					coefs[coef_index] = coef
				}
			}

			veff <- round(w_s(u,nu))
			k <- round(qt(0.977,df=veff), 2)
			uc <- sqrt(sum((u*coefs)^2))
			MPE <- sum(mpe)
			U <- uc*k

			uut_index <- which(object$value$variables$kind == "UUT")
			var_name <- object$value$variables[uut_index,][["name"]]
			
			VI_first_sample <- as.character(object$table_data[paste(var_name, "1")][[1]][1])

			VI <- get(var_name, row_env)

			with(row_env, eval(parse(text=object$value$formula)))
			e <- NA
			if(length(ls(row_env, pattern="^e$"))){
				e <- get("e", row_env)
			}

			VC <- as.numeric(VI) - as.numeric(e)

			new_row <- sapply(names(out_table_data), function(x){return(get(x))})
			out_table_data[i,] <- new_row

			row_details <- list(
				list(
					list(
						type = "raw",
						value = "Uncertainties components:"
						),
					list(
						type = "table",
						value = data.frame(var_name = var_names, x = x, u = u, coef = coefs, nu = nu, row.names=description)
						)
					)
				)
			details <- c(details, row_details)
			
		},
		error = function(c){ 
			row_details <- list(
				list(
					list(
						type = "errors", value = c$message
						)
					) 
				)
			# Access to global, set details with the errors
			details <<- c(details, row_details)
		},
		#warning = function(c){ print(c$message) },
		#message = function(c){ print(c$message) },
		finally = {})
		# Resume on error
		if(inherits(out, "error")){
			next
		} 
	}
	return(list(data = out_table_data, details = details))
}
