DATA_SECTION
  init_int n
  init_matrix x_obs(1,n,1,2)

PARAMETER_SECTION
  init_bounded_number kappa(0.1,650,-1)
  init_bounded_number a(0,100,-1)
  init_bounded_number b(0,50,-1)
  init_bounded_number sigma(0,50,1)
  random_effects_matrix x(1,n,1,2,1)
  objective_function_value f

PROCEDURE_SECTION
  f = 0.0;
  ll_1(kappa, a, b, sigma, x(1, 1), x(1, 2));
  cout << "1: " << "kappa: " << kappa << ", a: " << a << ", b: " << b << ", sigma: " << sigma << ", f: " << f << endl;
  ll_2(kappa, a, b, sigma, x(2, 1), x(2, 2), x(1, 1), x(1, 2));
  //cout << "2: " << "kappa: " << kappa << ", a: " << a << ", b: " << b << ", sigma: " << sigma << ", f: " << f << endl;
  for (int i = 3; i <= n; i++){
    ll_rest(i, kappa, a, b, sigma, x(i, 1), x(i, 2), x(i - 1, 1), x(i - 1, 2), x(i - 2, 1), x(i - 2, 2));
    //cout << i << ": " << "kappa: " << kappa << ", a: " << a << ", b: " << b << ", sigma: " << sigma << ", f: "<< f << endl;
  }
  

SEPARABLE_FUNCTION void ll_1(const dvariable& kappa, const dvariable& a, const dvariable& b, const dvariable& sigma, const dvariable& x_1, const dvariable& y_1)
  dvariable x_diff = x_1;
  dvariable y_diff = y_1;
  // Calculating distance moved
  dvariable dist = sqrt(square(x_diff) + square(y_diff) + 1e-150);
  // Subtracting log density of gamma distribution
  dvariable alpha = square(a)/square(b);
  dvariable beta = a/square(b);
  f -= alpha*log(beta) + (alpha - 1)*log(dist) - (beta*dist) - gammln(alpha);
  // Calculating bearing
  dvariable bear = atan(x_diff/y_diff);
  if (value(y_diff) < 0){
    bear += M_PI;
  } else if (value(x_diff) < 0){
    bear += 2*M_PI;
  }
  // Subtracting log density of von-Mises distribution
  // Working out bessel function
  dvariable y, bessi0;
  dvariable ax = fabs(kappa);	     
  if (value(ax) < 3.75){
    y = kappa*(1/3.75);
    y *= y;
    bessi0 = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    bessi0 = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  // The actual von-Mises density
  f -= kappa*cos(bear) - log(2*M_PI*bessi0);
  // Contribution from measurement error
  dvariable x_error = x_obs(1, 1) - x_1;
  dvariable y_error = x_obs(1, 2) - y_1;
  f -= -0.5*log(2*M_PI) - log(sigma) - square(x_error)/(2*square(sigma));
  f -= -0.5*log(2*M_PI) - log(sigma) - square(y_error)/(2*square(sigma));

SEPARABLE_FUNCTION void ll_2(const dvariable& kappa, const dvariable& a, const dvariable& b, const dvariable& sigma, const dvariable& x_1, const dvariable& y_1, const dvariable& x_2, const dvariable& y_2)
  dvariable x_diff = x_1 - x_2;
  dvariable y_diff = y_1 - x_2;   
  dvariable old_x_diff = x_2;
  dvariable old_y_diff = y_2;
  // Calculating distance moved
  dvariable dist = sqrt(square(x_diff) + square(y_diff) + 1e-150);
  // Subtracting log density of gamma distribution
  dvariable alpha = square(a)/square(b);
  dvariable beta = a/square(b);
  f -= alpha*log(beta) + (alpha - 1)*log(dist) - (beta*dist) - gammln(alpha);
  // Calculating bearing
  dvariable bear = atan(x_diff/y_diff);
  if (value(y_diff) < 0){
    bear += M_PI;
  } else if (value(x_diff) < 0){
    bear += 2*M_PI;
  }
  // Calculating old bearing
  dvariable old_bear = atan(old_x_diff/old_y_diff);   
  if (value(old_y_diff) < 0){
    old_bear += M_PI;
  } else if (value(old_x_diff) < 0){
    old_bear += 2*M_PI;
  }
  // Subtracting log density of von-Mises distribution
  // Working out bessel function
  dvariable y, bessi0;
  dvariable ax = fabs(kappa);	  
  if (value(ax) < 3.75){
    y = kappa*(1/3.75);
    y *= y;
    bessi0 = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    bessi0 = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  // The actual von-Mises density
  f -= kappa*cos(bear - old_bear) - log(2*M_PI*bessi0);
  // Contribution from measurement error
  dvariable x_error = x_obs(2, 1) - x_1;
  dvariable y_error = x_obs(2, 2) - y_1;
  f -= -0.5*log(2*M_PI) - log(sigma) - square(x_error)/(2*square(sigma));
  f -= -0.5*log(2*M_PI) - log(sigma) - square(y_error)/(2*square(sigma));

SEPARABLE_FUNCTION void ll_rest(int i, const dvariable& kappa, const dvariable& a, const dvariable& b, const dvariable& sigma, const dvariable& x_1, const dvariable& y_1, const dvariable& x_2, const dvariable& y_2, const dvariable& x_3, const dvariable& y_3)
  dvariable x_diff = x_1 - x_2;
  dvariable y_diff = y_1 - x_2;   
  dvariable old_x_diff = x_2 - x_3;
  dvariable old_y_diff = y_2 - y_3;
  // Calculating distance moved
  dvariable dist = sqrt(square(x_diff) + square(y_diff) + 1e-150);
  // Subtracting log density of gamma distribution
  dvariable alpha = square(a)/square(b);
  dvariable beta = a/square(b);
  f -= alpha*log(beta) + (alpha - 1)*log(dist) - (beta*dist) - gammln(alpha);
  // Calculating bearing
  dvariable bear = atan(x_diff/y_diff);
  if (value(y_diff) < 0){
    bear += M_PI;
  } else if (value(x_diff) < 0){
    bear += 2*M_PI;
  }
  // Calculating old bearing
  dvariable old_bear = atan(old_x_diff/old_y_diff);   
  if (value(old_y_diff) < 0){
    old_bear += M_PI;
  } else if (value(old_x_diff) < 0){
    old_bear += 2*M_PI;
  }
  // Subtracting log density of von-Mises distribution
  // Working out bessel function
  dvariable y, bessi0;
  dvariable ax = fabs(kappa);	  
  if (value(ax) < 3.75){
    y = kappa*(1/3.75);
    y *= y;
    bessi0 = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    bessi0 = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  // The actual von-Mises density
  f -= kappa*cos(bear - old_bear) - log(2*M_PI*bessi0);
  // Contribution from measurement error
  dvariable x_error = x_obs(i, 1) - x_1;
  dvariable y_error = x_obs(i, 2) - y_1;
  f -= -0.5*log(2*M_PI) - log(sigma) - square(x_error)/(2*square(sigma));
  f -= -0.5*log(2*M_PI) - log(sigma) - square(y_error)/(2*square(sigma));

TOP_OF_MAIN_SECTION
  arrmblsize=150000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(4600404);


