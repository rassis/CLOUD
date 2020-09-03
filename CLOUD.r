library(dplyr)
library(readr)
library(mvtnorm) # For multivariate Gaussian distribution
library(moments) # For generating moments
library(keras)
library(tensorflow)

GenerateTrainingData <- function(m, Nk, singlecopy_filename, duplicate_filename, training_prefix) {
	# Could make following into function parameters
	propPvC = 0.5  # For subfunctionalization
	logThetaMin = -4
	logThetaMax = 4
	logAlphaMin = 0
	logAlphaMax = 3
	logSigmaSqMin = -2
	logSigmaSqMax = 3

	scenarios <- c("cons", "neoparent", "neochild", "sub", "spec")
	tissues <- seq(1, m)

	training_data <- matrix(0, nrow = 5*Nk, ncol = 2 + 8*m) # m tissues with 3 expressionv alues and 5 parameters each, tPC, and class label

	set_of_tPC <- as_tibble(read.table(duplicate_filename, header = TRUE))  %>%
	  mutate(tPC = TPC / TPCA) %>%
	  select(tPC)
	set_of_tPC <- set_of_tPC$tPC

	for(sindex in 1:length(scenarios)) {
		scenario <- scenarios[sindex]


		for(i in 1:Nk) {
			tPCA <- 1
			tPC <- set_of_tPC[sample(1:length(set_of_tPC), 1)]

			alpha <- 10^runif(m, min = logAlphaMin, max = logAlphaMax)
			sigmaSq <- 10^runif(m, min = logSigmaSqMin, max = logSigmaSqMax)
			
			thetaP <- c()
			thetaC <- c()
			thetaA <- c()
			
			if(scenario == "cons") { # CONSERVATION
				thetaA <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaP <- thetaA
				thetaC <- thetaA
			}
			else if(scenario == "neoparent") { # NEO(P)
				thetaA <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaP <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaC <- thetaA
			}
			else if(scenario == "neochild") { # NEO(C)
				thetaA <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaC <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaP <- thetaA
			}
			else if(scenario == "sub") { # SUBFUNCTIONALIZATION        
				thetaP <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaC <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaA <- propPvC * thetaP + (1 - propPvC) * thetaC
			}
			else if(scenario == "spec") { # SPECIALIZATION
				thetaA <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaP <- runif(m, min = logThetaMin, max = logThetaMax)
				thetaC <- runif(m, min = logThetaMin, max = logThetaMax)
			}    

			expression_vec <- c()
			mu <- rep(0, 3)
			CovMat <- matrix(0, nrow = 3, ncol = 3)
			for(j in 1:m) {
				mu[1] = (1 - exp(-alpha[j] * tPC)) * thetaP[j] + exp(-alpha[j] * tPC) * thetaA[j]
				mu[2] = (1 - exp(-alpha[j] * tPC)) * thetaC[j] + exp(-alpha[j] * tPC) * thetaA[j]
				mu[3] = thetaA[j]
				CovMat[1, 1] = sigmaSq[j] / (2 * alpha[j])
				CovMat[2, 2] = CovMat[1, 1]
				CovMat[3, 3] = CovMat[1, 1]
				CovMat[1, 2] = exp(-2 * alpha[j] * tPC) * sigmaSq[j] / (2 * alpha[j])
				CovMat[2, 1] = CovMat[1, 2]
				CovMat[1, 3] = exp(-2 * alpha[j] * tPCA) * sigmaSq[j] / (2 * alpha[j])
				CovMat[3, 1] = CovMat[1, 3]
				CovMat[2, 3] = CovMat[1, 3]
				CovMat[3, 2] = CovMat[1, 3]
				expression_vec <- c(expression_vec, rmvnorm(1, mean = mu, sigma = CovMat))
			}
			
			rowIndex <- (sindex - 1)*Nk + i # get the row of the training dataset
			training_data[rowIndex, 1] = sindex
			training_data[rowIndex, 2] = tPC
			for(j in 1:(3*m)) {
			  training_data[rowIndex, 2 + j] = expression_vec[j]
			}
			for(j in 1:m) {
			  training_data[rowIndex, 2 + 3*m + j] = thetaP[j]
			}
			for(j in 1:m) {
			  training_data[rowIndex, 2 + 4*m + j] = thetaC[j]
			}
			for(j in 1:m) {
			  training_data[rowIndex, 2 + 5*m + j] = thetaA[j]
			}
			for(j in 1:m) {
			  training_data[rowIndex, 2 + 6*m + j] = log10(alpha[j])
			}
			for(j in 1:m) {
			  training_data[rowIndex, 2 + 7*m + j] = log10(sigmaSq[j])
			}
		}
	}

	column_labels <- c("Class", "tPC")
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("eP", j, sep = ""), paste("eC", j, sep = ""), paste("eA", j, sep = ""))
	}
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("ThetaP", j, sep = ""))
	}
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("ThetaC", j, sep = ""))
	}
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("ThetaA", j, sep = ""))
	}
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("LogAlpha", j, sep = ""))
	}
	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("LogSigmaSq", j, sep = ""))
	}


	colnames(training_data) <- column_labels

	write.table(training_data, file = paste(training_prefix, ".data", sep = ""), row.names = FALSE)
	
	GenerateFeatures(m, singlecopy_filename, paste(training_prefix, ".data", sep = ""), paste(training_prefix, ".features", sep = ""))
	GenerateClassifierResponse(paste(training_prefix, ".data", sep = ""), paste(training_prefix, ".classes", sep = ""))
	GeneratePredictorResponse(paste(training_prefix, ".data", sep = ""), paste(training_prefix, ".responses", sep = ""))
}

GenerateFeatures <- function(m, singlecopy_filename, input_filename, feature_filename) {
  minexp = 1e-4
  errorexp = 1e-5
  
  single <- read.table(singlecopy_filename, header = TRUE)
  for(i in 1:nrow(single)) {
    for(j in 1:ncol(single)) {
      single[i,j] <- log10(single[i,j] + minexp + rnorm(1, 0, errorexp)) # Transform data to log scale, accounting for expression of 0
    }
  }
  eS1S2dist <- sqrt( rowSums((single[,1:m] - single[,(m+1):(2*m)])^2)  ) # Euclidean distance between single-copy expression profiles
  maxS1S2dist <- max(abs(eS1S2dist))
  eS1S2cor <- c()
  for(i in 1:nrow(single)) {
    eS1S2cor[i] = cor(as.numeric(single[i,1:m]), as.numeric(single[i,(m+1):(2*m)]), method = "pearson") 
  }
  rm(single)

  features <- as_tibble( as.matrix(read.table(input_filename, header = TRUE)) ) %>% 
    select(tPC, starts_with("eP"), starts_with("eC"), starts_with("eA")) 
  
  eSumPC <- select(features, starts_with("eP")) + select(features, starts_with("eC"))
  colnames(eSumPC) <- paste("eSumPC", seq(1, m), sep = "")
  
  features <- bind_cols(features, as_tibble(eSumPC))
  
  rm(eSumPC)
  
  eDistPvC <- sqrt( rowSums( ( select(features, starts_with("eP")) - select(features, starts_with("eC")) )^2 ) ) # Euclidean distance between P and C
  eDistPvA <- sqrt( rowSums( ( select(features, starts_with("eP")) - select(features, starts_with("eA")) )^2 ) ) # Euclidean distance between P and A
  eDistCvA <- sqrt( rowSums( ( select(features, starts_with("eC")) - select(features, starts_with("eA")) )^2 ) ) # Euclidean distance between C and A
  eDistPCvA <- sqrt( rowSums( ( select(features, starts_with("eSumPC")) - select(features, starts_with("eA")) )^2 ) ) # Euclidean distance between PC and A
  
  features <- mutate(features, DistPvC = eDistPvC, DistPvA = eDistPvA, DistCvA = eDistCvA, DistPCvA = eDistPCvA)
  
  rm(eDistPvC)
  rm(eDistPvA)
  rm(eDistCvA)
  rm(eDistPCvA)
  
  features <- features %>%
    mutate(PBS_P = (DistPvC + DistPvA - DistCvA)/2) %>%
    mutate(PBS_C = (DistPvC + DistCvA - DistPvA)/2) %>%
    mutate(PBS_At = (DistPvA + DistCvA - DistPvC)/2) %>%
    rowwise() %>%
    mutate(RankDistPvC = mean(eS1S2dist < DistPvC)) %>%
    mutate(RankDistPvA = mean(eS1S2dist < DistPvA)) %>%
    mutate(RankDistCvA = mean(eS1S2dist < DistCvA)) %>%
    mutate(RankDistPCvA = mean(eS1S2dist < DistPCvA)) %>%
    mutate(DistPvC_m1 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m1 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m1 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m1 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m2 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m2 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m2 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m2 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m3 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m3 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m3 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m3 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m4 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m4 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m4 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m4 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m5 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m5 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m5 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m5 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m6 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m6 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m6 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m6 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m7 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m7 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m7 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m7 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvC_m8 = moment((eS1S2dist - DistPvC) / maxS1S2dist, order = 8, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPvA_m8 = moment((eS1S2dist - DistPvA) / maxS1S2dist, order = 8, central = FALSE, absolute = FALSE)) %>%
    mutate(DistCvA_m8 = moment((eS1S2dist - DistCvA) / maxS1S2dist, order = 8, central = FALSE, absolute = FALSE)) %>%
    mutate(DistPCvA_m8 = moment((eS1S2dist - DistPCvA) / maxS1S2dist, order = 8, central = FALSE, absolute = FALSE)) %>%
    ungroup()

  
  eCorPvC <- c()
  eCorPvA <- c()
  eCorCvA <- c()
  eCorPCvA <- c()

  # Pearson correlations
  for(i in 1:nrow(features)) {
    eCorPvC[i] <- cor( as.numeric(select(features[i,], starts_with("eP"))), as.numeric(select(features[i,], starts_with("eC"))), method = "pearson" )
    eCorPvA[i] <- cor( as.numeric(select(features[i,], starts_with("eP"))), as.numeric(select(features[i,], starts_with("eA"))), method = "pearson" )
    eCorCvA[i] <- cor( as.numeric(select(features[i,], starts_with("eC"))), as.numeric(select(features[i,], starts_with("eA"))), method = "pearson" )
    eCorPCvA[i] <- cor( as.numeric(select(features[i,], starts_with("eSumPC"))), as.numeric(select(features[i,], starts_with("eA"))), method = "pearson" )
  } 
  
  features <- mutate(features, CorPvC = eCorPvC, CorPvA = eCorPvA, CorCvA = eCorCvA, CorPCvA = eCorPCvA)
  
  rm(eCorPvC)
  rm(eCorPvA)
  rm(eCorCvA)
  rm(eCorPCvA)
  
  features <- features %>%
    rowwise() %>%
    mutate(RankCorPvC = mean(eS1S2cor < CorPvC)) %>%
    mutate(RankCorPvA = mean(eS1S2cor < CorPvA)) %>%
    mutate(RankCorCvA = mean(eS1S2cor < CorCvA)) %>%
    mutate(RankCorPCvA = mean(eS1S2cor < CorPCvA)) %>%
    mutate(CorPvC_m1 = moment(eS1S2cor - CorPvC, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m1 = moment(eS1S2cor - CorPvA, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m1 = moment(eS1S2cor - CorCvA, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m1 = moment(eS1S2cor - CorPCvA, order = 1, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m2 = moment(eS1S2cor - CorPvC, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m2 = moment(eS1S2cor - CorPvA, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m2 = moment(eS1S2cor - CorCvA, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m2 = moment(eS1S2cor - CorPCvA, order = 2, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m3 = moment(eS1S2cor - CorPvC, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m3 = moment(eS1S2cor - CorPvA, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m3 = moment(eS1S2cor - CorCvA, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m3 = moment(eS1S2cor - CorPCvA, order = 3, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m4 = moment(eS1S2cor - CorPvC, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m4 = moment(eS1S2cor - CorPvA, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m4 = moment(eS1S2cor - CorCvA, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m4 = moment(eS1S2cor - CorPCvA, order = 4, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m5 = moment(eS1S2cor - CorPvC, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m5 = moment(eS1S2cor - CorPvA, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m5 = moment(eS1S2cor - CorCvA, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m5 = moment(eS1S2cor - CorPCvA, order = 5, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m6 = moment(eS1S2cor - CorPvC, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m6 = moment(eS1S2cor - CorPvA, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m6 = moment(eS1S2cor - CorCvA, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m6 = moment(eS1S2cor - CorPCvA, order = 6, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m7 = moment(eS1S2cor - CorPvC, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m7 = moment(eS1S2cor - CorPvA, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m7 = moment(eS1S2cor - CorCvA, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPCvA_m7 = moment(eS1S2cor - CorPCvA, order = 7, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvC_m8 = moment(eS1S2cor - CorPvC, order = 8, central = FALSE, absolute = FALSE)) %>%
    mutate(CorPvA_m8 = moment(eS1S2cor - CorPvA, order = 8, central = FALSE, absolute = FALSE)) %>%
    mutate(CorCvA_m8 = moment(eS1S2cor - CorCvA, order = 8, central = FALSE, absolute = FALSE)) %>%  
    mutate(COrPCvA_m8 = moment(eS1S2cor - CorPCvA, order = 8, central = FALSE, absolute = FALSE)) %>%
    ungroup() %>%
    as.matrix()
  
  write.table(features, file = feature_filename, row.names = FALSE)
}

GenerateClassifierResponse <- function(input_filename, response_filename) {
  training_data <- as_tibble( as.matrix(read.table(input_filename, header = TRUE)) ) %>% 
    select(Class) %>%
    as.matrix()
  
  response <- matrix(0, nrow = nrow(training_data), ncol = 5)
  for(i in 1:nrow(training_data)) {
    response[i, training_data[i,1]] = 1
  }
  colnames(response) <- c("cons", "neoparent", "neochild", "sub", "spec")
  
  write.table(response, file = response_filename, row.names = FALSE)
}

GeneratePredictorResponse <- function(input_filename, response_filename) {
  training_data <- as_tibble( as.matrix(read.table(input_filename, header = TRUE)) ) %>% 
    select(starts_with("theta"), starts_with("LogAlpha"), starts_with("LogSigma")) %>%
    as.matrix()
  
  write.table(training_data, file = response_filename, row.names = FALSE)
}


ClassifierCV <- function(m, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix) {
  CV <- 5
  lambdas <- 10^seq(log_lambda_min, log_lambda_max, length = num_lambda)
  gammas <- seq(gamma_min, gamma_max, length = num_gamma)
  
  X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Y <- as.matrix(read.table(paste(training_prefix, ".classes", sep = ""), header = TRUE))
   
  # standardize the input for training
  X_means <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  for(j in 1:ncol(X)) {
	  X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }

  # CV-fold cross validation
  val_loss <- array(0, dim = c(length(gammas), length(lambdas)))

  # Randomly choose balanced training/validation sets per fold
  foldid_cons <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_neoparent <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_neochild <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_sub <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_spec <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid <- c(foldid_cons, foldid_neoparent, foldid_neochild, foldid_sub, foldid_spec)
  rm(foldid_cons)
  rm(foldid_neoparent)
  rm(foldid_neochild)
  rm(foldid_sub)
  rm(foldid_spec)
  
  
  # Perform K-fold CV, where K = CV
  for(curr_fold in 1:CV) {
	  Xval <- X[foldid == curr_fold, ]
	  Xtrain <- X[foldid != curr_fold, ]
	  Yval <- Y[foldid == curr_fold, ]
	  Ytrain <- Y[foldid != curr_fold, ]

	  # standardize the input for train and val based on train
	  temp_means <- colMeans(Xtrain)
	  temp_sds <- apply(Xtrain, 2, sd)
	  for(j in 1:ncol(Xtrain)) {
		  Xtrain[,j] = (Xtrain[,j] - temp_means[j]) / temp_sds[j]
		  Xval[,j] = (Xval[,j] - temp_means[j]) / temp_sds[j]
	  }
	  
	  for(i in 1:length(gammas)) {
		  for(j in 1:length(lambdas)) {
		  	model <- keras_model_sequential()
			  model %>%
				  layer_dense(units = 256, 
					  	        activation = 'relu',
						  	      kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
							        input_shape = c(ncol(Xtrain))) %>%
				  layer_dense(units = 128, 
					  		      activation = 'relu',
						  	      kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
				  layer_dense(units = 5,
							        activation = 'softmax',
							        kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))

			  model %>% 
			    compile(loss = 'categorical_crossentropy',
					        optimizer = optimizer_adam(),
					        metrics = c('categorical_crossentropy'))
					  
			  history <- model %>% 
			    fit(Xtrain, Ytrain,
				      epochs = num_epochs,
				      batch_size = batchsize,
				      validation_data = list(Xval, Yval),
				      verbose = 0) # verbose = 0  ensures it is silent
					  
			  val_loss[i,j] <- val_loss[i,j] + as.data.frame(history) %>%
				  filter(data == "validation", metric == "categorical_crossentropy") %>%
				  select(value) %>%
				  min()
		  }
	  }
	
	  rm(Xval)
	  rm(Xtrain)
	  rm(Yval)
  	rm(Ytrain)
  }

  val_loss <- val_loss / CV
	
  gamma_opt = gammas[ which(val_loss == min(val_loss), arr.ind = TRUE)[1] ]
  lambda_opt = lambdas[ which(val_loss == min(val_loss), arr.ind = TRUE)[2] ]
  cv_results <- matrix(0, 1, 3)
  cv_results[,1] = min(val_loss)
  cv_results[,2] = gamma_opt
  cv_results[,3] = lambda_opt
  colnames(cv_results) <- c("Loss", "Gamma", "Lambda")
  
  write.table(cv_results, file = paste(training_prefix, ".classifier_cv", sep = ""), row.names = FALSE)
  
  ClassifierFit(m, batchsize, num_epochs, training_prefix)
}

PredictorCV <- function(m, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix) {
  CV <- 5
  lambdas <- 10^seq(log_lambda_min, log_lambda_max, length = num_lambda)
  gammas <- seq(gamma_min, gamma_max, length = num_gamma)
  
  X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Y <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE))
  
  # standardize the input for training
  X_means <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  Y_means <- colMeans(Y)
  Y_sds <- apply(Y, 2, sd)
  for(j in 1:ncol(X)) {
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  for(j in 1:ncol(Y)) {
    Y[,j] = (Y[,j] - Y_means[j]) / Y_sds[j]
  }
  
  # CV-fold cross validation
  val_loss <- array(0, dim = c(length(gammas), length(lambdas)))
  
  # Randomly choose balanced training/validation sets per fold
  foldid_cons <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_neoparent <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_neochild <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_sub <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid_spec <- sample(rep(seq(CV), length = nrow(X)/5))
  foldid <- c(foldid_cons, foldid_neoparent, foldid_neochild, foldid_sub, foldid_spec)
  rm(foldid_cons)
  rm(foldid_neoparent)
  rm(foldid_neochild)
  rm(foldid_sub)
  rm(foldid_spec)
  
  
  # Perform K-fold CV, where K = CV
  for(curr_fold in 1:CV) {
    Xval <- X[foldid == curr_fold, ]
    Xtrain <- X[foldid != curr_fold, ]
    Yval <- Y[foldid == curr_fold, ]
    Ytrain <- Y[foldid != curr_fold, ]
    
    # standardize the input for train and val based on train
    temp_Xmeans <- colMeans(Xtrain)
    temp_Xsds <- apply(Xtrain, 2, sd)
    temp_Ymeans <- colMeans(Ytrain)
    temp_Ysds <- apply(Ytrain, 2, sd)
    for(j in 1:ncol(Xtrain)) {
      Xtrain[,j] = (Xtrain[,j] - temp_Xmeans[j]) / temp_Xsds[j]
      Xval[,j] = (Xval[,j] - temp_Xmeans[j]) / temp_Xsds[j]
    }
    for(j in 1:ncol(Ytrain)) {
      Ytrain[,j] = (Ytrain[,j] - temp_Ymeans[j]) / temp_Ysds[j]
      Yval[,j] = (Yval[,j] - temp_Ymeans[j]) / temp_Ysds[j]
    }
    
    for(i in 1:length(gammas)) {
      for(j in 1:length(lambdas)) {
        model <- keras_model_sequential()
        model %>%
          layer_dense(units = 256, 
                      activation = 'relu',
                      kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
                      input_shape = c(ncol(Xtrain))) %>%
          layer_dense(units = 128, 
                      activation = 'relu',
                      kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
          layer_dense(units = ncol(Ytrain),
                      activation = 'linear',
                      kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
        
        model %>% 
          compile(loss = 'mean_squared_error',
                  optimizer = optimizer_adam(),
                  metrics = c('mean_squared_error'))
        
        history <- model %>% 
          fit(Xtrain, Ytrain,
              epochs = num_epochs,
              batch_size = batchsize,
              validation_data = list(Xval, Yval),
              verbose = 0) # verbose = 0  ensures it is silent
        
        val_loss[i,j] <- val_loss[i,j] + as.data.frame(history) %>%
          filter(data == "validation", metric == "mean_squared_error") %>%
          select(value) %>%
          min()
      }
    }
    
    rm(Xval)
    rm(Xtrain)
    rm(Yval)
    rm(Ytrain)
  }
  
  val_loss <- val_loss / CV
  
  gamma_opt = gammas[ which(val_loss == min(val_loss), arr.ind = TRUE)[1] ]
  lambda_opt = lambdas[ which(val_loss == min(val_loss), arr.ind = TRUE)[2] ]
  cv_results <- matrix(0, 1, 3)
  cv_results[,1] = min(val_loss)
  cv_results[,2] = gamma_opt
  cv_results[,3] = lambda_opt
  colnames(cv_results) <- c("Loss", "Gamma", "Lambda")
  
  write.table(cv_results, file = paste(training_prefix, ".predictor_cv", sep = ""), row.names = FALSE)
  
  PredictorFit(m, batchsize, num_epochs, training_prefix)
}

ClassifierFit <- function(m, batchsize, num_epochs, training_prefix) {
  cv_results <- read.table(paste(training_prefix, ".classifier_cv", sep = ""), header = TRUE)
  gamma_opt <- cv_results$Gamma
  lambda_opt <- cv_results$Lambda
  
  X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Y <- as.matrix(read.table(paste(training_prefix, ".classes", sep = ""), header = TRUE))
  
  # standardize the input for training
  X_means <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  for(j in 1:ncol(X)) {
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  
  std_params <- data.frame("Xmeans" = X_means, "Xsds" = X_sds)
  write.table(std_params, paste(training_prefix, ".X_stdparams", sep = ""), row.names = FALSE)
  rm(std_params)
  
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 256, 
                activation = 'relu',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
                input_shape = c(ncol(X))) %>%
    layer_dense(units = 128, 
                activation = 'relu',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
    
    layer_dense(units = 5,
                activation = 'softmax',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
  
  model %>% 
    compile(loss = 'categorical_crossentropy',
            optimizer = optimizer_adam(),
            metrics = c('categorical_crossentropy', 'accuracy'))
  
  history <- model %>% 
    fit(X, Y,
        epochs = num_epochs,
        batch_size = batchsize,
        verbose = 0) # verbose = 0  ensures it is silent
  
  model %>% save_model_hdf5(paste(training_prefix, ".classifier.hdf5", sep = ""))
}

PredictorFit <- function(m, batchsize, num_epochs, training_prefix) {
  cv_results <- read.table(paste(training_prefix, ".predictor_cv", sep = ""), header = TRUE)
  gamma_opt <- cv_results$Gamma
  lambda_opt <- cv_results$Lambda
  
  X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Y <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE))
  
  # standardize the input and output for training
  X_means <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  Y_means <- colMeans(Y)
  Y_sds <- apply(Y, 2, sd)
  for(j in 1:ncol(X)) {
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  for(j in 1:ncol(Y)) {
    Y[,j] = (Y[,j] - Y_means[j]) / Y_sds[j]
  }
  
  std_params <- data.frame("Xmeans" = X_means, "Xsds" = X_sds)
  write.table(std_params, paste(training_prefix, ".X_stdparams", sep = ""), row.names = FALSE)
  std_params <- data.frame("Ymeans" = Y_means, "Ysds" = Y_sds)
  write.table(std_params, paste(training_prefix, ".Y_stdparams", sep = ""), row.names = FALSE)
  rm(std_params)
  
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 256, 
                activation = 'relu',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
                input_shape = c(ncol(X))) %>%
    layer_dense(units = 128, 
                activation = 'relu',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
    
    layer_dense(units = ncol(Y),
                activation = 'linear',
                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
  
  model %>% 
    compile(loss = 'mean_squared_error',
            optimizer = optimizer_adam(),
            metrics = 'mean_squared_error')
  
  history <- model %>% 
    fit(X, Y,
        epochs = num_epochs,
        batch_size = batchsize,
        verbose = 0) # verbose = 0  ensures it is silent
  
  model %>% save_model_hdf5(paste(training_prefix, ".predictor.hdf5", sep = ""))
}

CLOUDClassify <- function(training_prefix, testing_prefix) {
  X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
  
  # standardize the input for testing
  std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
  X_means <- c(std_params$Xmeans)
  X_sds <- c(std_params$Xsds)
  for(j in 1:ncol(X)) {
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  rm(std_params)
  
  model <- load_model_hdf5(paste(training_prefix, ".classifier.hdf5", sep = ""))
  
  Yest_num <- model %>% predict_classes(X) 
  
  Yest <- data.frame("Class" = ifelse(Yest_num == 0, "cons", 
                                     ifelse(Yest_num == 1, "neoparent",
                                            ifelse(Yest_num == 2, "neochild",
                                                  ifelse(Yest_num == 3, "sub", "spec")))))
  
  write.table(Yest, paste(testing_prefix, ".classifications", sep = ""), row.names = FALSE)
}

CLOUDPredict <- function(training_prefix, testing_prefix) {
  X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
  
  # standardize the input for testing
  std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
  X_means <- c(std_params$Xmeans)
  X_sds <- c(std_params$Xsds)
  for(j in 1:ncol(X)) {
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  rm(std_params)
  
  model <- load_model_hdf5(paste(training_prefix, ".predictor.hdf5", sep = ""))
  
  std_params <- read.table(paste(training_prefix, ".Y_stdparams", sep = ""), header = TRUE)
  Y_means <- c(std_params$Ymeans)
  Y_sds <- c(std_params$Ysds)
  rm(std_params)
  
  Yest_std <- model %>% predict(X)
  Yest <- Yest_std
  
  # un-standardize the responses
  for(j in 1:ncol(Yest)) {
    Yest[, j] = Yest_std[,j] * Y_sds[j] + Y_means[j]
  }
  
  Y <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE))
  colnames(Yest) <- colnames(Y)
  rm(Y)
  
  write.table(Yest, paste(testing_prefix, ".predictions", sep = ""), row.names = FALSE)
  
  
}