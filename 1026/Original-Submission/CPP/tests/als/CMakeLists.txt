set(tests test_als)
foreach (test ${tests})
    create_test(test.als.${test} ${test})
endforeach ()