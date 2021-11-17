function TTEST_ERROR_ASSERT_WORKER
% TTEST_ERROR_ASSERT_WORKER()
% Default function called when an ASSERT_ fails
% just calles assert( false )
%
% See also: TTEST_ERROR_EXPECT_WORKER
    assert( false );
end