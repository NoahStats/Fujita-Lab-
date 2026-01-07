#https://www.spsanderson.com/steveondata/posts/2025-03-24/
library(foreach)


##### basic syntax  -----------------------------------------------------------
foreach(variable = sequence) %do% {
#code   
}


##### Example1 : Basic Iteration  ---------------------------------------------
#sums squares of numbers from 1 to 5
result <- foreach(i = 1:5) %do% {
  i^2
}
result
#foreach returns a list while traditional default for-loop doesn't collect results 


##### Example 2: Combining results  -------------------------------------------
#Sum the squares of numbers from 1 to 5 
total = foreach(i = 1:5, .combine = '+') %do% {
  i^2
}
total


#### Example 3: Multiple Input Sequences  ------------------------------------

results = foreach(a = 1:3, b = 4:6) %do% { 
  a*b
  }
results 


#### Example 4: Working with Dataframes  -------------------------------------

# sample data frame 
df = data.frame(
  id = 1:3, 
  values = c(10, 20, 30)
)
#calculate a new column based on values 
results = foreach(id = df$id, val = df$values) %do% { 
  data.frame(id = id, value = val, squared = val^2)
}
print(results)

#combine results into a single dataframe
combined_df = rbind(results)
combined_df = do.call(rbind, results) 
#do.call(rbind, ) is the correct way to combine lists of dataframes into single df
print(combined_df)


# Example 5: Parallel Processing  -----------------------------------------
library(doParallel)

#Register parallel backend 
cores = detectCores() -1 # use one less than available cores 
registerDoParallel(cores)

# Perform parallel computation;%dopar% for 'parallel' 
results = foreach(i = 1:10, .combine = 'c') %dopar% {
  #simulate a computation-heavy task
  Sys.sleep(1) #Sleep for 1second
  i^2
}

#Stop the parallel backend 
stopImplicitCluster()

print(results)


# Practice  ---------------------------------------------------------------
#Register parallel backend
library(doParallel)
cores = detectCores()-1
registerDoParallel(cores)

#Perform parallel computation 
results = foreach(i = 1:5, .combine = 'c') %dopar% {
  factorial(i)
}
print(results)

#Stop the parallel backend 
stopImplicitCluster()


# Example 6: Exporting Variables and Packages  ----------------------------
library(foreach)
library(doParallel)

#Register Parallel backend 
registerDoParallel(2)

#Define a function and variable in the main environment 
my_function = function(x){
  return(x^2 + y)
}
y = 10

# Use .export and .packages to make dependencies available 
results = foreach(i = 1:5, 
                  .export = c('my_function', 'y'), #export variables and functions
                  .packages = 'stats') %dopar% { #export packages
                    my_function(i) + mean(c(i, i+1)) # mean() from stats package
                  }

stopImplicitCluster()
print(results)


# Example 7: Handing Errors with .errorhandling  --------------------------
results = foreach(i = c(1,2,0,4,5), 
                  .combine = 'c',
                  .errorhandling = 'remove') %do% {
                    10/i #Will cause division by zero error for i = 0
                  }
print(results)


# Example 8 Converting For Loop into foreach() ----------------------------
result = foreach(i = 1:5, .combine = 'c') %do% {
  i ^3 
}


# Example 9: Comparing Sequential and Parallel foreach() ------------------

library(foreach)
library(doParallel)
library(tictoc)
# Function to calculate prime numbers up to n 
is_prime = function(n) {
  if(n <= 1) return(FALSE)
  if(n <= 3) return (TRUE)
  if(n %% 2  ==0 || n%%3 == 0) return(FALSE) 
  i = 5 
  while(i*i<= n) {
    if(n%% i == 0 || n%% (i+2) == 0) return(FALSE)
    i = i+6
  } 
  return(TRUE)
              
}
#large numbers to check for primality
numbers = 100000000 + 1:8

# Sequential execution 
tic('Sequential')
seq_result = foreach(num = numbers, .combine = 'c') %do% {
  is_prime(num)
}
toc()

#Parallel execution 
registerDoParallel(4)
tic('Parallel')
par_result = foreach(num = numbers, .combine = 'c') %dopar% {
  is_prime(num)
}
toc()
stopImplicitCluster()

