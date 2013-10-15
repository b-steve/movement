DATA_SECTION
  init_int n
  init_matrix x_obs(1,n,1,2)

PARAMETER_SECTION
  init_bounded_number kappa(0,20)
  init_bounded_number a(0,100)
  init_bounded_number b(0,50)
  init_bounded_number sigma(0,50)
  random_effects_matrix x(1,n,1,2)
  number x_diff
  number y_diff
  number x_error
  number y_error
  number dist
  number bear
  number alpha
  number beta
  objective_function_value f

PROCEDURE_SECTION
  f = 0.0;
  for (int i = 1; i <= n; i++){
    // Calculating x and y movement
    if (i == 1){
      x_diff = x(i, 1);
      y_diff = x(i, 2);
    } else {
      x_diff = x(i, 1) - x(i - 1, 1);
      y_diff = x(i, 2) - x(i - 1, 2);
    }
    // Calculating distance moved
    dist = sqrt(square(x_diff) + square(y_diff));
    // Subtracting log density of gamma distribution
    alpha = square(a)/square(b);
    beta = a/square(b);
    f -= alpha*log(beta) + (alpha - 1)*log(dist) - (beta*dist) - gammln(alpha);
    // Calculating bearing
    df1b2variable rat;
    rat = x_diff/dist;
    bear = asin(rat);
    dvariable bess;
    bess = bessi0(kappa);
    // Subtractiong log density of von-Mises distribution
    f -= kappa*cos(bear) - log(2*M_PI*bess);
    // Contribution from measurement error
    x_error = x_obs(i, 1) - x(i, 1);
    y_error = x_obs(i, 2) - x(i, 2);
    f -= -0.5*log(2*M_PI) - log(sigma) - square(x_error)/(2*square(sigma));
    f -= -0.5*log(2*M_PI) - log(sigma) - square(y_error)/(2*square(sigma));
  }

GLOBALS_SECTION
  #include <bessel.cpp>
