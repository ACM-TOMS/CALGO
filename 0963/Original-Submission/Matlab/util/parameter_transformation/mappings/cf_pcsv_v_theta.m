function phi = cf_pcsv_v_theta(decode, omega, theta, y_t, tau)
  [mu, A, lambda_0, kappa, theta_param, sigma, rho] = ...
    decode_pcsv_param(theta, decode);
  phi = cf_pcsv_v(mu, A, lambda_0, kappa, theta_param, sigma, rho, ...
                omega, y_t, tau);
end
