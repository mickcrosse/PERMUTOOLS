function test_permuvartest2
%TEST_PERMUVARTEST2  Unit tests for permuvartest2.m.

%   © 2018-2023 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate test data
rng(42);
x = randn(30,20);
y = randn(30,20);
tail = {'both','right','left'};

for i = 1:numel(tail)

    % Perform F-test
    [~,p1,ci1,stats1] = vartest2(x,y,'tail',tail{i});

    % Perform permutation test
    [f2,p2,ci2,stats2] = permuvartest2(x,y,'tail',tail{i},'correct',0);

    % Assert that results are the same or similar
    assert(max(abs(p1-p2))<0.3)
    assert(max(abs(ci1(:)-ci2(:)))<5)
    assert(any(stats1.fstat==f2))
    assert(any(stats1.df1==stats2.df1))
    assert(any(stats1.df2==stats2.df2))

end

disp('All unit tests for permuvartest2.m passed.')