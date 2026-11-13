# following Jie's code, this is a "class" describing the fused ellipse structure
# for a given semi-major axis and aspect ration
FusedEllipse <- function(a = 1, ratio = 2.0) {
  
  b <- a / ratio
  e2 <- 1 - (b / a)^2
  
  # Integrand for the elliptic integral of the second kind
  f <- function(u) sqrt(1 - e2 * sin(u)^2)
  
  # Complete elliptic integral, quarter circumference
  E_quarter <- integrate(f, 0, pi / 2,
                         rel.tol = 1e-10)$value
  
  # Incomplete elliptic integral
  E <- function(u) {
    integrate(f, 0, u, rel.tol = 1e-10)$value
  }
  
  # Invert the incomplete elliptic integral
  inverse_E <- function(s) {
    
    if (s <= 0) return(0)
    if (s >= 1) return(pi / 2)
    
    uniroot(
      function(u) E(u) - s * E_quarter,
      interval = c(0, pi / 2),
      tol = 1e-10
    )$root
  }
  
  # solve for thetaprime given theta 
  theta_prime <- function(theta) {
    
    t <- theta %% (2 * pi)
    
    upper <- t <= pi
    s <- ifelse(upper, t, 2 * pi - t) / pi
    
    u <- vapply(s, inverse_E, numeric(1))
    
    ifelse(upper, u, pi - u)
  }
  
  # following Jie, additional functionality for nucleoid information
  nucleoid <- function(theta) {
    
    tp <- theta_prime(theta)
    
    cbind(
      x = -a * sin(tp),
      y = b * cos(tp)
    )
  }
  
  list(
    a = a,
    b = b,
    e2 = e2,
    E_quarter = E_quarter,
    theta_prime = theta_prime,
    nucleoid = nucleoid
  )
}

# two fused mitochondria with different a/b
ellipse.1 <- FusedEllipse(a = 1, ratio = 1)
ellipse.2 <- FusedEllipse(a = 1, ratio = 4)

# theta values for plotting
theta <- seq(0, pi, length.out = 100)

# convert to thetaprime for both ellipses
tp.1 <- ellipse.1$theta_prime(theta)
tp.2 <- ellipse.2$theta_prime(theta)

# dataframe with the results
df = rbind(data.frame(theta = theta, thetaprime = tp.1, e = 1),
           data.frame(theta = theta, thetaprime = tp.2, e = 4))
      
# plot the corresponding regions, converting to degrees     
sf = 3
png("jie-model.png", width=600*sf, height=200*sf, res=72*sf)
ggarrange(
ggplot(df, aes(x = theta*180/pi, y = thetaprime*180/pi, color=factor(e))) + 
  geom_line() + theme_minimal() + theme(legend.position = "none") + 
  labs(x = "θ / degrees", y = "φ / degrees") ,
ggplot(df, aes(x = theta*180/pi, y = sin(thetaprime), color=factor(e))) + 
  geom_line() + theme_minimal() + theme(legend.position = "none") + 
  labs(x = "θ / degrees", y = "d/a") ,
nrow = 1)
dev.off()
