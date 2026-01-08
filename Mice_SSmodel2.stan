data{
  int<lower = 1> n; 
  int<lower = 1> T; 
  vector[T] y[n]; 
  // no spike probability 


//start from system equations and then separate Population level parameters and individual level parameters
// Also make sure to separate mean, sd and noise. 

parameters{

    // system equations 
    vector[T] mu_t;
    real sd<lower = 0> sd_mu_t; 

    // Population level parameters 
    real level;  //mean
    real A; 
    real B;
    real phi;   
    real<lower= 0> sd_level;  //sd
    real<lower =0> sd_trend_it;
    real<lower = 0> sd_A;
    real<lower =0> sd_B; 
    real<lower = > sd_phi; 

    // Individual level parameters
    vector[n] level_i; 
    vector[n] A_i; 
    vector[n] B_i; 
    vector[n] phi_i; 

    // noise parameters 
    real<lower = 0> sd_y; 
}



model{
// system equations 
for(t in 2:T){
  mu_t[t] ~ normal(mu_t[t-1], sd_mu_t);
}

//Strucutral Equations (Hierarchical Structures)
for(i in 1:n){
    level_i[i] ~ normal(level, sd_level); 
    mu_it[i,t] ~ normal(mu_t[t], sd_trend_it);
    A_i[i] ~ normal(A, sd_A);
    P_i[i] ~ normal(P, sd_P); 
    phi_i[i] ~ normal(phi, sd_phi); 
}


//Observation Equations

for(i in 1:n){
    
    real season[i] = A_i[i] * sin(2*pi()/P_i[i] + phi_i[i]);
        
        for(t in 1:T){
            mu0 = level_i[i] + mu_it[i,t] + season[i];
            y[i][t] ~ normal(mu0, sd_y);

    }
}
}







generated quantities{
    matrix[n,T] y_rep; 

    for(i in 1:n){
        for(t in 1:t){
            real season[i] = A_i[i] * sin(2*pi()/P_i[i] + phi_[i]);
            real mu0 = level_i[i] + mu_it[i,t] + season[i]; 
            y_rep[i,t] = normal_rng(mu0, sd_y);
        }
    }
}