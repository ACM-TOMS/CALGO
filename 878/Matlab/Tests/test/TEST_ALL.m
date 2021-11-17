% TEST_ALL Script to run all tests for Vauto package

disp 'RUNNING PRIMARY AND SECONDARY TESTS FOR VAUTO PACKAGE'
% TEST ROUTINES NEEDED FOR COMPLETE-DATA LOG-LIKELIHOOD:
test_der2array_mds_set
test_find_CGW
test_vyw quiet
test_omega_building
test_omega_forward
test_omega quiet
test_varma_llc
test_parmatprod quiet
% Derivative routines:
test_find_CGW_deriv
test_vyw_deriv quiet
test_trisolve_deriv
test_chol_deriv quiet
test_omega_deriv quiet
test_varma_llc_deriv
%
% TEST MISSING VALUE ROUTINES:
test_lambda_multiply quiet
test_find_Sm
test_find_lambda_om
test_find_Vhat
test_product_deriv
test_omega_building miss
test_atba_c
test_profile
test_profile_ar
test_C_multiply
test_omega_ar
test_var_ll
test_var_ll_jac
test_varma_llm
test_varma_jac
test_varma_llm_deriv
%
% MISC TESTS
test_varma_sim quiet
%
disp('ALL TESTS SUCCEEDED')
