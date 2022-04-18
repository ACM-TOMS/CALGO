# TTEST Suite

## Introduction: Why *TTESTs*
The toolbox is written with the following ides in mind.
- A test should be easy to write, to prevent bugs in the test suits
- A test suite should be easy to write, to foster the writing of tests
- A test run should not stop just because one test fails, allowing to fix multiple bugs in a single test run
- When a test fails, as much information as possible should be provided
- Tests should be portable between Matlab version
  (So far this is most likely not fulfiled and the *TTESTs* may only work well with R2016b onwards)
- The test suite shall not replace Matlab features, and work with most Matlab test styles 
  (So far this is most likely not fulfiled and the *TTESTs* may only work well with script based tests)

In particular, the *TTESTs* are designed to mimic the UI of *googletests*.


## Type of tests
There are  mostly two types of tests: `EXPECT_` and `ASSERT_` tests.
In the following only `EXPECT_` tests are described. The syntax for `ASSERT_` tests is always the same.

       EXPECT_FAIL()                           always fails
       EXPECT_SUCCEED()                        always succeeds

       EXPECT_FALSE( val )                     tests whether val is (implicitely) false
       EXPECT_TRUE( val )                      tests whether val is (implicitely) true

       EXPECT_EQ( lhs, rhs )                   tests isequal(lhs,rhs)
       EXPECT_NE( lhs, rhs )                   tests ~isequal(lhs,rhs)
       EXPECT_LE/EXPECT_LT/EXPECT_GE/EXPECT_GT( lhs, rhs )         tests lhs <= / < / >= / > rhs
       EXPECT_STREQ/EXPECT_STRNE               tests whether two strings are equal or unequal using strcmp
       EXPECT_STRCASEEQ/EXPECT_STRCASENE       tests whether two strings are equal or unequal (up to case) using strcmpi
                                               All four functions accept char arrays (e.g. 'text') and strings (e.g. "text")
                                               Note that strings are only supported (in a sensible way) from R2018a onwards

       EXPECT_PRED/EXPECT_NPRED( pred, arg1, arg2, ..., argn )      Tests whether pred( arg1, arg2, ..., argn ) is true / false
                                               These functions differ from the GTEST macros, since the number of arguments does not has to be provided in the functionname
       EXPECT_NEAR( lhs, rhs, bd )             tests norm( lhs - rhs ) <= bd
       EXPECT_DOUBLE_EQ( lhs, rhs, [factor] )  tests norm( lhs - rhs ) <= factor*eps(lhs), i.e. the two values only differ up to factor many ULP's. By default, factor = 4
       EXPECT_FLOAT_EQ/EXPECT_SINGLE_EQ( lhs, rhs, [factor] )      the same as _DOUBLE_EQ but the test is made with single-precision ULP's

       EXPECT_WARNING( 'cmd' )                 tests whether the cmd omits at least one warning
       EXPECT_WARNING( 'cmd', warnid1, ..., warnidn )      tests whether the cmd (given as string or char array) throws at least one warning with id's warnid1, ..., warnidn
                                                               If no warning is from warnid1, ..., warnidn is thrown, the test fails
                                                               Note that: Currently this function cannot distinguish well between warnings and errors.

       EXPECT_ERROR( 'cmd' )                                   The same was for EXPECT_WARNING, but with errors
       EXPECT_ERROR( 'cmd', warnid1, ..., warnidn )        

       EXPECT_NO_THROW( 'cmd' )                         tests whether the cmd throws no warning nor error.
       EXPECT_NO_THROW( 'cmd', warnid1, ..., warnidn )  tests whether the cmd throws no warning nor error or any of warnid1, ..., warndidn
       

## Behaviour in failure:

- If an `EXPECT_`/`ASSERT_` test fails because the condition is not met (but the test makes syntactily sense), then 
  - the behaviour depends on the variables `TTEST_ERROR_EXPECT`/`TTEST_ERROR_ASSERT`. By default `TTEST_ERROR_EXPECT='exp'`, `TTEST_ERROR_ASSERT='ass'`
    - If the variable is `'exp'`, then `TTEST_ERROR_EXPECT_WORKER()` is called, which only prints out a summary of the failed test.
       In particular, `assert()` is not called, and thus, Matlab does not recognize this as a failure of the test
    - If the variable is `'ass'`, then `TTEST_ERROR_ASSERT()` is called, which calls `assert(false)`, which is catched and rethrown (to have nicer Matlab output)
    - In all other cases, the variable `TTEST_ERROR_EXPECT`/`TTEST_ERROR_ASSERT` is passed to feval, which allows to define a custom error handler
  - the global variable `TTEST_ERRORFLAG` is set to true

- If a test fails, because the test is impossible e.g. `EXPECT_GE( @sum, @prod )`, `EXPECT_EQ(2)`, then (in most cases) an error is thrown which is not caught.
  This behaviour may be changed in future releases
  
## Compatibility
The interface may be subject to be changed until v1.0
