function test_permuttest
%TEST_PERMUTTEST  Unit tests for permuttest.m.

%   © 2018-2023 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate test data
rng(42);
x = randn(30,20);
y = randn(30,20);
tail = {'both','right','left'};

for i = 1:numel(tail)

    % Perform t-test
    [~,p1,ci1,stats1] = ttest(x,y,'tail',tail{i});

    % Perform permutation test
    [t2,p2,ci2,stats2] = permuttest(x,y,'tail',tail{i},'correct',0,...
        'verbose',0);

    % Assert that results are the same or similar
    assert(max(abs(p1-p2))<0.02)
    assert(max(abs(ci1(:)-ci2(:)))<0.02)
    assert(any(stats1.tstat==t2))
    assert(any(stats1.df==stats2.df))
    assert(any(stats1.sd==stats2.sd))

end

disp('All unit tests for permuttest.m passed.')