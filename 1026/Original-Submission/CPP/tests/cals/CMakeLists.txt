set(tests test_cals)
foreach (test ${tests})
    create_test(test.cals.${test} ${test})
endforeach ()
