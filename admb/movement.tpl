DATA_SECTION
  init_int n
  init_matrix x_obs(1,n,1,2)

PARAMETER_SECTION
  init_bounded_number kappa(0,650,1)
  init_bounded_number a(0,1000,1)
  init_bounded_number b(0,1000,1)
  init_bounded_number sigma(0,1000,1)
  random_effects_matrix x(1,n,1,2,1)
  number x_diff
  number y_diff
  number sign_x_diff
  number sign_y_diff
  number x_dummy
  number y_dummy
  number both_dummy
  number x_error
  number y_error
  number dist
  number bear
  number old_bear
  number alpha
  number beta
  number bessi0
  number ax
  number y
  objective_function_value f

PROCEDURE_SECTION
  f = 0.0;
  for (int i = 1; i <= n; i++){
    // Calculating x and y movement
    if (i == 1){
      x_diff = x(i, 1);
      y_diff = x(i, 2);
      old_bear = 0;
    } else {
      x_diff = x(i, 1) - x(i - 1, 1);
      y_diff = x(i, 2) - x(i - 1, 2);
      old_bear = bear;
    }
    sign_x_diff = x_diff/fabs(x_diff);
    sign_y_diff = y_diff/fabs(y_diff);
    // Calculating distance moved
    dist = sqrt(square(x_diff) + square(y_diff));
    // Subtracting log density of gamma distribution
    alpha = square(a)/square(b);
    beta = a/square(b);
    f -= alpha*log(beta) + (alpha - 1)*log(dist) - (beta*dist) - gammln(alpha);
    // Calculating dummy variables
    // x_dummy is 1 if x_diff is positive and y_diff is negative
    x_dummy = (sign_x_diff - sign_y_diff)*(sign_x_diff + square(sign_y_diff))/4;
    // y_dummy is 1 if x_diff is negative and y_diff is positive
    y_dummy = (sign_y_diff - sign_x_diff)*(sign_y_diff + square(sign_x_diff))/4;
    // both_dummy is 1 if both x_diff and y_diff are negative
    both_dummy = (-sign_x_diff - sign_y_diff)*(square(sign_x_diff) - sign_y_diff)/4;
    // Calculating bearing
    bear = asin(x_diff/dist);
    // Adjustments depending on the bearing's quadrant
    bear += x_dummy*2*(0.5*M_PI - bear);
    bear += y_dummy*2.0*M_PI;
    bear += both_dummy*(M_PI - 2.0*bear);
    // Subtracting log density of von-Mises distribution
    // Working out bessel function
    ax = fabs(kappa);	  
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
    x_error = x_obs(i, 1) - x(i, 1);
    y_error = x_obs(i, 2) - x(i, 2);
    f -= -0.5*log(2*M_PI) - log(sigma) - square(x_error)/(2*square(sigma));
    f -= -0.5*log(2*M_PI) - log(sigma) - square(y_error)/(2*square(sigma));
  }

TOP_OF_MAIN_SECTION
  arrmblsize=150000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(4600404);


