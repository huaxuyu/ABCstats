lambdaOpt = function(data_seq, group_vector, L1 = -3, L2 = 3){
  # All lambda values to be tested
  lambda_seq = seq(L1, L2, 0.01)

  # A variable to indicate if data are transformed or not
  note = "Transformed"

  # Generate a vector to store the product of p values from all groups
  p_seq = c()

  # Acquire all group names
  gv = unique(group_vector)

  # Generate a list to store quantitative data from each group
  ds = c()

  # Separate data into individual groups
  for (i in 1:length(gv)) {
    ds[[i]] = data_seq[group_vector == gv[i]]

    # Check identical values
    if (all(ds[[i]] == ds[[i]][1])) {
      result = list("data_trans" = data_seq,
                    "lambda" = "No ABC due to identical values within one or more groups",
                    "note" = "No ABC due to identical values within one or more groups")
      return(result)
    }
  }

  # Optimize lambda value by calculating the pooled data normality
  for (i in 1:length(lambda_seq)) {
    p_values = c()
    if(lambda_seq[i] == 0){
      for (j in 1:length(gv)) {
        p_values[j] = shapiro.test(log( ds[[j]] ))$p.value
      }
    } else{
      for (j in 1:length(gv)) {
        group_data = ( ds[[j]]^lambda_seq[i] - 1 ) / lambda_seq[i]
        p_values[j] = shapiro.test(group_data)$p.value
      }
    }
    p_seq[i] = sum(log10(p_values))
  }

  # Calculate the pooled data normality for original data
  for (i in 1:length(gv)) {
    p_values[i] = shapiro.test( ds[[i]] )$p.value
  }

  # Find the best lambda
  p_seq = append(p_seq, sum(log10(p_values)))
  index = match(max(p_seq),p_seq)
  lambda_final = round(lambda_seq[index],2)

  # Perform data transformation using the best lambda
  if(is.na(lambda_final)){
    note = "Original"
  } else if(lambda_final == 0) {
    data_seq = log(data_seq)
  } else {data_seq = (data_seq^lambda_final-1)/lambda_final}

  # Return the results in a list
  result = list("data_trans" = data_seq,
                "lambda" = lambda_final,
                "note" = note)
  return(result)
}

# Imputation
impute = function(df, impt="default"){
  if (impt=="default") {
    df[is.na(df)] = 0
    for (i in 1:nrow(df)){
      data_seq = as.numeric(df[i,])
      # Gap-filling by adding a constant small value
      zero_number = sum(data_seq == 0)
      if(zero_number != 0){
        data_seq = data_seq + max(data_seq)/1000
      }
      df[i,] = data_seq
    }
  } else if (impt=="knn") {
    message("Imputation with KNN algorithm...")
    message("This process may take minites to finish.")
    df = t(VIM::kNN(t(df)))[1:nrow(df),]
  } else if (impt=="addition"){
    df[is.na(df)] = 0
    for (i in 1:nrow(df)){
      data_seq = as.numeric(df[i,])
      # Gap-filling by adding a constant small value
      zero_number = sum(data_seq == 0)
      if(zero_number != 0){
        data_seq[data_seq == 0] = max(data_seq)/1000
      }
      df[i,] = data_seq
    }
  } else {
    message("Invalid method name. Default method is used.")
    df = impute(df)
  }
  return(df)
}
