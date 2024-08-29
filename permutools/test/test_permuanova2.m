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
[p1,Table1,stats1] = anova1(x,group,'off');

% Perform permutation one-way ANOVA
[~,p2,~,stats2,Table2] = permuanova1(x,group);

% Assert that results are the same or similar
assert(abs(p1-p2)<0.005)
assert(round(Table1{6},10)==round(Table2{6},10))
assert(round(Table1{7},10)==round(Table2{7},10))
assert(round(Table1{8},10)==round(Table2{8},10))
assert(Table1{10}==Table2{10})
assert(Table1{11}==Table2{11})
assert(Table1{12}==Table2{12})
assert(round(Table1{14},10)==round(Table2{14},10))
assert(round(Table1{15},10)==round(Table2{15},10))
assert(any(ismember(stats1.gnames,stats2.gnames)))
assert(any(stats1.n==stats2.n))
assert(any(stats1.means==stats2.means))
assert(stats1.df==stats2.df)
assert(round(stats1.s,10)==round(stats2.s,10))

disp('All unit tests for permuanova1.m passed.')