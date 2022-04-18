% TEST_PRIMARY  Script to run primary tests of Vauto package

disp 'TESTING PRIMARY VAUTO FUNCTIONS VAR_LL, VARMA_LLC, VARMA_LLM AND VARMA_SIM'
test_varma_llc
test_varma_llc_deriv
test_var_ll
test_var_ll_jac
test_varma_llm
test_varma_jac
test_varma_llm_deriv
test_varma_sim quiet
disp('ALL PRIMARY TESTS SUCCEEDED')
