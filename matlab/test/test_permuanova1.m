function test_permuanova1
%TEST_PERMUANOVA1  Unit tests for permuanova1.m.

%   © 2018-2023 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate test data
rng(42);
x = randn(30,5);
group = 1:5;

% Perform one-way ANOVA
[p1,tbl1,stats1] = anova1(x,group,'off');

% Perform permutation one-way ANOVA
[~,p2,~,stats2,tbl2] = permuanova1(x,group);

% Assert that results are the same or similar
assert(abs(p1-p2)<0.005)
assert(round(tbl1{6},10)==round(tbl2{6},10))
assert(round(tbl1{7},10)==round(tbl2{7},10))
assert(round(tbl1{8},10)==round(tbl2{8},10))
assert(tbl1{10}==tbl2{10})
assert(tbl1{11}==tbl2{11})
assert(tbl1{12}==tbl2{12})
assert(round(tbl1{14},10)==round(tbl2{14},10))
assert(round(tbl1{15},10)==round(tbl2{15},10))
assert(any(ismember(stats1.gnames,stats2.gnames)))
assert(any(stats1.n==stats2.n))
assert(any(stats1.means==stats2.means))
assert(stats1.df==stats2.df)
assert(round(stats1.s,10)==round(stats2.s,10))

disp('All unit tests for permuanova1.m passed.')