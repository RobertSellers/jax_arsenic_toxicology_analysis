data {
    int<lower=0> N; 
    real y[N];
    real<lower=0> x[N];
    int<lower=0> idc[N]; 
    int<lower=0> J;
    int<lower=0> idr[N]; 
    int<lower=0> K;
    real pb[J]; 
    real pc; 
    real pd; 
    real pe[J];
} 
parameters { 
    real<lower=0> 
    sigmasq_y;
    real slope[J]; 
    real lasy; 
    real uasy; 
    real<lower=0> ed[J];
    real rslope[K]; 
    real rlasy[K]; 
    real ruasy[K]; 
    real<lower=0> sigmasq_slope; 
    real<lower=0> sigmasq_lasy; 
    real<lower=0> sigmasq_uasy;
    real rnslope; 
    real rnlasy; 
    real rnuasy;
} 
transformed parameters {
    real<lower=0> sigma_y;
    real mu[N];
    real<lower=0> sigma_slope; 
    real<lower=0> sigma_lasy; 
    real<lower=0> sigma_uasy; 
    sigma_slope <- sqrt(sigmasq_slope); 
    sigma_lasy <- sqrt(sigmasq_lasy); 
    sigma_uasy <- sqrt(sigmasq_uasy); 
    sigma_y <- sqrt(sigmasq_y); 
    for(i in 1:N){ 
        mu[i] <- ( lasy  + rlasy[idr[i]] )  + ( ( uasy  + ruasy[idr[i]] ) - ( lasy  + rlasy[idr[i]] ) ) / (1 + exp(-exp( ( slope[idc[i]]  + rslope[idr[i]] ) ) * (log(x[i]/ ed[idc[i]] ))))^exp( 0 );
    }
} 
model {
    slope ~ normal(pb, 1) ; 
    lasy ~ normal(pc, 100) ; 
    uasy ~ normal(pd, 100) ; 
    ed ~ normal(pe, 2) ;
    sigmasq_y ~ inv_gamma(0.001, 0.001) ; 
    y ~ normal(mu, sigma_y);
    rslope ~ normal(0, sigma_slope); 
    rlasy ~ normal(0, sigma_lasy); 
    ruasy ~ normal(0, sigma_uasy);
    sigmasq_slope ~ inv_gamma(0.001, 0.001) ; 
    sigmasq_lasy ~ inv_gamma(0.001, 0.001) ; 
    sigmasq_uasy ~ inv_gamma(0.001, 0.001) ;
    rnslope ~ normal(0, sigma_slope); 
    rnlasy ~ normal(0, sigma_lasy); 
    rnuasy ~ normal(0, sigma_uasy);
}
generated quantities {
    real residuals[N];
    real log_lik[N];
    real pslope[J]; 
    real plasy; 
    real puasy; 
    real ped[J]; 
    real passym;
    for (i in 1:N){
        residuals[i] <- y[i] - mu[i];log_lik[i] <- normal_log(y[i], mu[i], sigma_y);
    }
    for (j in 1:J) pslope[j]  <-  slope[j]  + rnslope ; 
    plasy <- lasy  + rnlasy ; 
    puasy  <-  uasy  + rnuasy ; 
    for (j in 1:J) ped[j]  <-  ed[j]  ; passym  <-  0  ;
}
