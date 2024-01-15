# Suggest how to use transformation to parametrize a parameter that is...

# a.) only positive
pars<- list(
  log_theta = log(2) # the value of theta will be = 2
)
f<- function(pars) {
  getAll(pars)
  theta<- exp(log_theta)
  #
  # The rest of your negative loglikelihood goes here
  #
}

# b.) only negative
pars<- list(
  nlog_theta = log(-1 * (-2)) # the value of theta will be = -2
)
f<- function(pars) {
  getAll(pars)
  theta<- -exp(nlog_theta)
  #
  # The rest of your negative loglikelihood goes here
  #
}

# c.) between 2 and 5
pars<- list(
  working_theta = qlogis((1 / 3) * (4.8 - 2)) # the value of theta will be 4.8
)
f<- function(pars) {
  getAll(pars)
  theta<- 3 * plogis(working_theta) + 2
  #
  # The rest of your negative loglikelihood goes here
  #
}

# d.) an increasing vector
pars<- list(
  working_theta = c(
    -4, # the first element of theta will be -4
    log(2 - (-4)), # the second element of theta will be 2
    log(7 - (-4) - exp(log(2 - (-4) )) )# the third element of theta will be 7
  )
)
pars<- list(
  getAll(pars)
  theta<- 0 * working_theta # Create a vector the same size as working_theta
  theta[1]<- working_theta[1]
  theta[2]<- working_theta[1] + exp(working_theta[2])
  theta[3]<- working_theta[1] + exp(working_theta[2]) + exp(working_theta[3])
  #
  # The rest of your negative loglikelihood goes here
  #
)
