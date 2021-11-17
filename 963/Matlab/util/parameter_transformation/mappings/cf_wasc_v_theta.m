function phi = cf_wasc_v_theta(decode, omega, theta, y_t, tau)
  [mu, Sigma_0, M, Q, rho, beta] = decode_wasc_param(theta, decode);
  phi = cf_wasc_v(mu, Sigma_0, M, Q, rho, beta, omega, y_t, tau);
end
