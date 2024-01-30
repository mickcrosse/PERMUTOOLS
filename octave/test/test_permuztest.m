function test_permuztest
%TEST_PERMUZTEST  Unit tests for permuztest.m.

%   © 2018-2023 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate test data
rng(42);
x = randn(30,20);
m = 0;
sigma = 1;
tail = {'both','right','left'};

for i = 1:numel(tail)

    % Perform Z-test
    [~,p1,ci1,zval1] = ztest(x,m,sigma,'tail',tail{i});

    % Perform permutation test
    [zval2,p2,ci2] = permuztest(x,m,sigma,'tail',tail{i},'correct',0,...
        'verbose',0);

    % Assert that results are the same or similar
    assert(max(abs(p1-p2))<0.1)
    assert(max(abs(ci1(:)-ci2(:)))<0.2)
    assert(any(zval1==zval2))

end

disp('All unit tests for permuztest.m passed.')